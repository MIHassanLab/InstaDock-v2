import warnings
warnings.filterwarnings("ignore")

import os
import sys
import torch
import pickle
import hashlib
import ID_utils
import numpy as np
import pandas as pd
from rdkit import Chem
from pathlib import Path
from chemprop import models
from sklearn.svm import SVC, SVR
from rdkit.Chem import Descriptors
from lightning import pytorch as pl
from lightgbm import LGBMClassifier
from chemprop import data, featurizers
from PyQt5.QtCore import QThread, pyqtSignal
from rdkit.Chem import Crippen, rdMolDescriptors
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.tree import DecisionTreeRegressor, DecisionTreeClassifier
from sklearn.neighbors import KNeighborsRegressor, KNeighborsClassifier
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor, AdaBoostRegressor, RandomForestClassifier, GradientBoostingClassifier, AdaBoostClassifier

def get_models_path():
    if sys.platform.startswith("win"):
        base_dir = os.getenv("LOCALAPPDATA", str(Path.home()))
    else:
        base_dir = str(Path.home() / ".local" / "share")
    return Path(base_dir) / "IDv2Models"

def compute_file_hash(file_path):
    hash_sha256 = hashlib.sha256()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_sha256.update(chunk)
    return hash_sha256.hexdigest()

def load_reference_hashes(ref_file_path):
    ref_hashes = {}
    if not os.path.exists(ref_file_path):
        ID_utils.displayError(f"Reference hash file not found: {ref_file_path}")
        return None
    with open(ref_file_path, "r") as f:
        for line in f:
            if "::" in line:
                rel_path, file_hash = line.strip().split("::")
                ref_hashes[rel_path] = file_hash
    return ref_hashes

def verify_folder_integrity(models_path, reference_hashes):
    models_path = Path(models_path)
    local_files = {}
    for root, _, files in os.walk(models_path):
        for file in files:
            file_path = Path(root) / file
            rel_path = str(file_path.relative_to(models_path))
            local_files[rel_path] = compute_file_hash(file_path)
    
    for ref_file, ref_hash in reference_hashes.items():
        if ref_file not in local_files:
            return False
        if local_files[ref_file] != ref_hash:
            return False
    
    for local_file in local_files:
        if local_file not in reference_hashes:
            return False
    return True

class ModelLoader:
    @staticmethod
    def load_catboost_model(pkl_path):
        with open(pkl_path, 'rb') as file:
            model = pickle.load(file)
        return model

    @staticmethod
    def load_dpmnn_model(model_path):
        # Assuming models.MPNN.load_from_checkpoint exists in your module
        mpnn = models.MPNN.load_from_checkpoint(model_path, map_location=torch.device('cpu'))
        return mpnn

    @staticmethod
    def load_meta_model(pkl_path):
        with open(pkl_path, 'rb') as file:
            model = pickle.load(file)
        return model


class predictADMET(QThread):
    finishedSignal = pyqtSignal(object)
    progressBarSignal = pyqtSignal(int)
    toastSignal = pyqtSignal(list)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # flag to ensure we only toast once per run
        self._toast_emitted = False

    def smiles(self, smiles):
        self.data = smiles

    meta_models_regression = {
        'Ridge': Ridge(),
        'RandomForest': RandomForestRegressor(),
        'SVR': SVR(),
        'KNeighbors': KNeighborsRegressor(),
        'ExtraTrees': ExtraTreesRegressor(),
        'AdaBoost': AdaBoostRegressor(),
        'DecisionTree': DecisionTreeRegressor()
    }

    meta_models_classification = {
        'LogisticRegression': LogisticRegression(),
        'RandomForest': RandomForestClassifier(),
        'SVM': SVC(),
        'GradientBoosting': GradientBoostingClassifier(),
        'KNeighbors': KNeighborsClassifier(),
        'AdaBoost': AdaBoostClassifier(),
        'DecisionTree': DecisionTreeClassifier(),
        'LightGBM': LGBMClassifier()
    }

    # --- Global options ---
    OPTIONS = ["ob", "pgp_inhibitor", "pgp_substrate", "hia", "caco2_reg", "bbb", "fu", "ppb",
               "vd", "logd", "logp", "logs", "pka", "pkb", "cyp2d6_inhibitor", "cyp2d6_substrate",
               "cyp3a4_inhibitor", "cyp3a4_substrate", "oatp1b1", "oatp1b3", "oct2", "t0.5", "cl",
               "bcrp", "ames", "carcinogenicity", "dili", "micronucleus_tox", "bp", "mp"]

    # Map options to proper names for DataFrame columns (if needed)
    OPTIONS_PROPER_NAMES = {opt: opt.upper() for opt in OPTIONS}  # customize as needed

    # --- Data Preprocessing Function ---
    def preprocess_data(self, smile_input):
        # Accept SMILES string or list
        smile = smile_input if isinstance(smile_input, str) else smile_input[0]

        # Validate SMILES first
        mol = Chem.MolFromSmiles(smile)
        if mol is None:
            if not self._toast_emitted:
                self.toastSignal.emit([2000, f"Invalid SMILES found"])
                self._toast_emitted = True
            return None, None

        # Sanitize after validation
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            if not self._toast_emitted:
                self.toastSignal.emit([2000, f"Sanitization failed!"])
                self._toast_emitted = True
            return None, None

        # Create featurizer with sanitized molecule
        featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
        
        # Calculate RDKit descriptors
        data_rdkit = Descriptors.CalcMolDescriptors(mol, missingVal=1.0)
        data_rdkit = list(data_rdkit.values())

        # Create test data from sanitized SMILES
        test_data = [data.MoleculeDatapoint.from_smi(Chem.MolToSmiles(mol))]
        test_dataset = data.MoleculeDataset(test_data, featurizer=featurizer)
        test_loader = data.build_dataloader(test_dataset, shuffle=False)

        return data_rdkit, test_loader

    # --- Prediction Function ---
    def predict_with_models(self, catboost_model, mpnn_model, meta_model, smile_input, regression_flag):
        """
        Makes ensemble predictions using CatBoost, MPNN, and a meta model.
        Returns:
            final_predictions: numpy array of predictions.
        """
        processed_rdkit, processed_mpnn = self.preprocess_data(smile_input)
        if processed_rdkit is None or processed_mpnn is None:
            return None

        # Predict with CatBoost model
        catboost_preds = catboost_model.predict(processed_rdkit)

        # Predict with MPNN using PyTorch Lightning
        with torch.inference_mode():
            trainer = pl.Trainer(
                enable_progress_bar=False,
                accelerator="cpu",
                devices=1,
                logger=False
            )
            mpnn_preds = trainer.predict(mpnn_model, processed_mpnn)

        # Process predictions based on task type
        mpnn_preds = np.concatenate(mpnn_preds, axis=0)
        if regression_flag:
            mpnn_preds = mpnn_preds.reshape(-1)
        else:
            # For classification, apply thresholding at 0.5
            mpnn_preds = (mpnn_preds > 0.5).astype(int)

        # Combine features for meta-model input
        combined_features = np.column_stack([catboost_preds, mpnn_preds])
        final_predictions = meta_model.predict(combined_features).flatten()

        return final_predictions

    def get_mol(self, smiles):
        """Convert SMILES to RDKit molecule without computing coordinates."""
        return Chem.MolFromSmiles(smiles)

    def calc_HBA(self, smiles):
        mol = self.get_mol(smiles)
        return rdMolDescriptors.CalcNumHBA(mol) if mol else None

    def calc_HBD(self, smiles):
        mol = self.get_mol(smiles)
        return rdMolDescriptors.CalcNumHBD(mol) if mol else None

    def calc_MW(self, smiles):
        mol = self.get_mol(smiles)
        return Descriptors.MolWt(mol) if mol else None

    def calc_logp(self, smiles):
        mol = self.get_mol(smiles)
        return Descriptors.MolLogP(mol) if mol else None

    def calc_rotatable_bonds(self, smiles):
        mol = self.get_mol(smiles)
        return rdMolDescriptors.CalcNumRotatableBonds(mol) if mol else None

    def calc_aromatic_rings(self, smiles):
        mol = self.get_mol(smiles)
        return rdMolDescriptors.CalcNumAromaticRings(mol) if mol else None

    def calc_fraction_csp3(self, smiles):
        mol = self.get_mol(smiles)
        return rdMolDescriptors.CalcFractionCSP3(mol) if mol else None

    def calc_nitro_groups(self, smiles):
        """Count nitro groups (SMARTS pattern)."""
        mol = self.get_mol(smiles)
        if not mol:
            return None
        nitro_smarts = Chem.MolFromSmarts('[NX3](=O)=O')
        return len(mol.GetSubstructMatches(nitro_smarts))

    def calc_halogen_groups(self, smiles):
        """Count halogen atoms (F, Cl, Br, I)."""
        mol = self.get_mol(smiles)
        if not mol:
            return None
        return sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() in {'F', 'Cl', 'Br', 'I'})

    def calc_aniline_groups(self, smiles):
        """Count aniline groups (SMARTS pattern)."""
        mol = self.get_mol(smiles)
        if not mol:
            return None
        aniline_smarts = Chem.MolFromSmarts('[NX3;H2,H1]c1ccccc1')
        return len(mol.GetSubstructMatches(aniline_smarts))

    def calc_benzene_rings(self, smiles):
        """Count benzene rings using SMARTS pattern."""
        mol = self.get_mol(smiles)
        if not mol:
            return None
        benzene_smarts = Chem.MolFromSmarts('c1ccccc1')
        return len(mol.GetSubstructMatches(benzene_smarts))

    def calc_ghose_rule(self, smiles):
        mol = self.get_mol(smiles)
        if mol:
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            mr = Crippen.MolMR(mol)
            return "Pass" if (160 <= mw <= 480 and -0.4 <= logp <= 5.6 and 40 <= mr <= 130) else "Fail"
        return "N/A"

    def calc_veber_rule(self, smiles):
        mol = self.get_mol(smiles)
        if mol:
            rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
            tpsa = rdMolDescriptors.CalcTPSA(mol)
            return "Pass" if (rot_bonds <= 10 and tpsa <= 140) else "Fail"
        return "N/A"

    def calc_pfizer_rule(self, smiles):
        mol = self.get_mol(smiles)
        if mol:
            logp = Descriptors.MolLogP(mol)
            tpsa = rdMolDescriptors.CalcTPSA(mol)
            return "Fail" if (logp > 3 and tpsa < 75) else "Pass"
        return "N/A"

    def calc_gsk_rule(self, smiles):
        mol = self.get_mol(smiles)
        if mol:
            logp = Descriptors.MolLogP(mol)
            tpsa = rdMolDescriptors.CalcTPSA(mol)
            return "Pass" if (logp <= 4 and tpsa >= 75) else "Fail"
        return "N/A"

    def calc_golden_triangle(self, smiles):
        mol = self.get_mol(smiles)
        if mol:
            logp = Descriptors.MolLogP(mol)
            tpsa = rdMolDescriptors.CalcTPSA(mol)
            return "Pass" if (-2 <= logp <= 5 and 40 <= tpsa <= 130) else "Fail"
        return "N/A"

    def calc_egan_rule(self, smiles):
        mol = self.get_mol(smiles)
        if mol:
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            tpsa = rdMolDescriptors.CalcTPSA(mol)
            return "Pass" if (mw <= 400 and logp <= 2 and tpsa <= 80) else "Fail"
        return "N/A"

    def calc_lipinski_rule(self, smiles):
        mol = self.get_mol(smiles)
        if mol:
            hbd = rdMolDescriptors.CalcNumHBD(mol)
            hba = rdMolDescriptors.CalcNumHBA(mol)
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            return "Pass" if (mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10) else "Fail"
        return "N/A"

    def calc_tpsa(self, smiles):
        mol = self.get_mol(smiles)
        return rdMolDescriptors.CalcTPSA(mol) if mol else None

    # --- Main Prediction Pipeline ---
    def run_prediction_pipeline(self, df, model_base_path):
        """
        Runs the entire prediction pipeline.
        
        Args:
            data_path (str): Path to the CSV file containing input data (SMILES).
            model_base_path (str): Base directory where model subfolders (by option) reside.
        
        Returns:
            results_dict (dict): Dictionary where keys are option names and values are lists of predictions.
        """
        # Assuming the CSV has a column named 'SMILES'
        smiles_list = df.iloc[:, 0]
        
        # Names of regression tasks
        regression_names = ["caco2_reg", "bbb", "fu", "ppb", "vd", "cl", "bp", "logd", "logp", "logs", "mp", "pka", "pkb"]
        
        results_dict = {"SMILES": smiles_list.tolist()}
        
        # Loop over each option and its proper name
        for i, option in enumerate(self.OPTIONS):
            regression_flag = option in regression_names
            
            folder_path = os.path.join(model_base_path, option)
            catboost_model_path = os.path.join(folder_path, 'Catboost.pkl')
            dmpnn_model_path = os.path.join(folder_path, 'NN.ckpt')
            meta_model_path = os.path.join(folder_path, 'Meta_model.pkl')
            
            # Load models using ModelLoader
            catboost_model = ModelLoader.load_catboost_model(catboost_model_path)
            mpnn_model = ModelLoader.load_dpmnn_model(dmpnn_model_path)
            meta_model = ModelLoader.load_meta_model(meta_model_path)
            
            # Apply prediction for each SMILES in the dataset
            option_predictions = []
            for smile in smiles_list:
                # Wrap the smile in a list to match expected input format
                prediction = self.predict_with_models(catboost_model, mpnn_model, meta_model, [smile], regression_flag)
                if prediction is None:
                    option_predictions.append(None)
                else:
                    # Assuming a single prediction value is returned per SMILES
                    option_predictions.append(prediction[0])
            
            # Save predictions in the dictionary using the proper name label
            proper_name = self.OPTIONS_PROPER_NAMES.get(option, option)
            results_dict[proper_name] = option_predictions
            progress = int((i + 1) / len(self.OPTIONS) * 100)
            self.progressBarSignal.emit(progress)

        # Calculate additional molecular properties
        moreproperties = {
            "HBA": self.calc_HBA,
            "HBD": self.calc_HBD,
            "MW": self.calc_MW,
            "LOGP": self.calc_logp,
            "Ghose rule": self.calc_ghose_rule,
            "GSK Rule": self.calc_gsk_rule,
            "Golden Triangle": self.calc_golden_triangle,
            "Lipinski Rule": self.calc_lipinski_rule,
            "Egan's rule": self.calc_egan_rule,
            "Pfizer Rule": self.calc_pfizer_rule,
            "Veber": self.calc_veber_rule,
            "Rotatable Bonds": self.calc_rotatable_bonds,
            "Number of aromatic rings": self.calc_aromatic_rings,
            "Nitro groups": self.calc_nitro_groups,
            "Halogen groups": self.calc_halogen_groups,
            "Aniline groups": self.calc_aniline_groups,
            "Benzene rings": self.calc_benzene_rings,
            "FractionCSP3": self.calc_fraction_csp3,
            "TPSA": self.calc_tpsa,
        }

        for property_name, function in moreproperties.items():
            valueList = []
            for smile in smiles_list:
                valueList.append(function(smile))
            results_dict[property_name] = valueList        

        return results_dict

    def run(self):    
        df = pd.DataFrame(self.data)
        self.toastSignal.emit([2000, "Verifying AI models..."])
        models_path = get_models_path()
        model_base_path = str(models_path)
        self.toastSignal.emit([2000, "Running prediction pipeline..."])
        results = self.run_prediction_pipeline(df, model_base_path)
        self.finishedSignal.emit(results)
