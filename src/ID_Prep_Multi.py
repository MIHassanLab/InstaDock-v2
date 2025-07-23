import os
import re
import sys
import shutil
import ID_utils
import tempfile
import pyKVFinder
import IDv2Resource
from glob import glob
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from PyQt5.QtCore import QThread, pyqtSignal, QProcess, QEventLoop


LIGAND_TAGS = [' UNL ', ' UNK ', ' LIG ', ' <0> ', ' DRG ', ' INH ', ' NAG ', ' SO4 ', ' ATP ', ' ADP ', ' AMP ', ' HEM ', ' FMN ', ' FAD ', ' NAD ', ' GDP ', ' GTP ', ' SAM ']

STD_RESIDUES = [' ALA ', ' ARG ', ' ASN ', ' ASP ', ' CYS ', ' GLN ', ' GLU ', ' GLY ', ' HIS ', ' ILE ', ' LEU ', ' LYS ', ' MET ', ' PHE ', ' PRO ', ' SER ', ' THR ', ' TRP ', ' TYR ', ' VAL ']

STD_RESIDUES2 = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']

key = b'XaQ5nlmxWQyPf9Q9bMjki5M7yEiXgOmTDRpvtTRWbik='

def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    if hasattr(sys, '_MEIPASS'):
        # When bundled by PyInstaller, _MEIPASS contains the temporary extraction path.
        return os.path.join(sys._MEIPASS, relative_path)
    return os.path.join(os.path.abspath("."), relative_path)


def inspectProtein(targetWD, method):
    # 1. Find all PDBs…
    pdbfiles = glob(os.path.join(targetWD, "**", "*.pdb"), recursive=True)
    # 2. Exclude any .pdb that has a matching .pdbqt
    filtered = []
    for pdb in pdbfiles:
        rel_parts = os.path.normpath(pdb).split(os.sep)
        if "Backup" in rel_parts:
            continue
        if os.path.isfile(pdb[:-4] + ".pdbqt"):
            continue
        filtered.append(pdb)
    pdbfiles = filtered

    if not pdbfiles:
        return "No PDB files found (or all have .pdbqt)!", "Error"

    all_msgs, all_titles = [], []
    for pdb in pdbfiles:
        filename = os.path.basename(pdb)
        # GPU-specific processing
        if method.upper() == 'GPU':
            proc = QProcess()
            proc.setProcessChannelMode(QProcess.MergedChannels)
            exe = resource_path("mk_prepare_receptor.exe")
            proc.start(exe, ["--read_pdb", pdb])
            proc.waitForFinished(-1)
            out = proc.readAllStandardOutput().data().decode()
            if 'error' in out.lower() or 'failed' in out.lower():
                all_msgs.append(f"{filename} cannot be processed using GPU methods. Please try again with CPU methods.")
                all_titles.append("GPU Error")
                continue  # skip chemical checks if GPU failed

        # read and analyze
        with open(pdb, 'r') as f:
            lines = f.readlines()
        content = "".join(lines)

        ligand_found = any(tag in content for tag in LIGAND_TAGS)
        std_res_found = any(res in content for res in STD_RESIDUES)
        has_cocrystal = ligand_found and std_res_found

        atomnos, chains = [], set()
        hetcount, anomalies = 0, []
        altloc, ins, hcount = False, False, 0

        for i, L in enumerate(lines):
            rec = L[0:6].strip()
            if rec in ("ATOM","HETATM"):
                occ, tmp = L[54:60].strip(), L[60:66].strip()
                if occ and not _is_float(occ):
                    anomalies.append(f"Line {i+1}: invalid occupancy '{occ}'")
                if tmp and not _is_float(tmp):
                    anomalies.append(f"Line {i+1}: invalid temp '{tmp}'")
                rnm = L[17:20].strip()
                if len(rnm)!=3:
                    anomalies.append(f"Line {i+1}: unusual residue name '{rnm}'")
                a=L[16].strip()
                if a and a!="A": altloc=True
                c=L[26].strip()
                if c: ins=True

                if rec=="ATOM":
                    num=L[22:26].strip()
                    if num:
                        if num.isdigit(): atomnos.append(int(num))
                        else: anomalies.append(f"Line {i+1}: invalid res num '{num}'")
                    cid=L[21].strip()
                    if cid: chains.add(cid)
                    else: anomalies.append(f"Line {i+1}: missing chain ID")
                    name=L[12:16].strip()
                    if name.startswith("H"): hcount+=1
                else:
                    hetcount+=1

        seq = sorted(set(atomnos))
        gaps = [str(n) for n in range(seq[0], seq[-1]+1) if n not in seq] if seq else []

        issues = []
        if has_cocrystal: issues.append("Found co-crystallized ligand")
        if gaps: issues.append(f"Missing residues: {len(gaps)} gaps ({','.join(gaps)})")
        if len(chains)>1: issues.append(f"Multiple chains: {','.join(sorted(chains))}")
        if hetcount: issues.append(f"HETATM records: {hetcount}")
        if altloc: issues.append("Alternate location indicators")
        if ins: issues.append("Insertion codes detected")
        if hcount: issues.append(f"Hydrogen atoms: {hcount}")

        if not issues and not anomalies:
            all_msgs.append(f"{filename} is ready for docking!")
            all_titles.append("OK")
        else:
            lines = [f"PDB analysis for {filename}:"]
            if issues:
                lines.append("Issues:")
                for item in issues:
                    lines.append(f"  - {item}")
            if anomalies:
                lines.append("Anomalies:")
                for item in anomalies:
                    lines.append(f"  - {item}")
            # join into one multiline string
            all_msgs.append("\n".join(lines))
            all_titles.append("Warning")

    overall = "Warning" if any(t != "OK" for t in all_titles) else "OK"

    # 3. Now handle the .pdbqt files (if you still want to report them):
    pdbqtfiles = glob(os.path.join(targetWD, "**", "*.pdbqt"), recursive=True)
    for pdbqt in pdbqtfiles:
        # skip Backup folder here too, if desired
        rel_parts = os.path.normpath(pdbqt).split(os.sep)
        if "Backup" in rel_parts:
            continue
        all_msgs.append(f"{os.path.basename(pdbqt)} is ready for docking!")
        all_titles.append("OK")

    return "\n\n----------------\n\n".join(all_msgs), overall

def _is_float(s):
    try:
        float(s)
        return True
    except:
        return False

class prepareTargetM1(QThread):
    requestResidueDialog = pyqtSignal()
    residuesReady = pyqtSignal(str)
    toastSignal = pyqtSignal(list)
    finishedSignal = pyqtSignal(list)
    
    # Meeko always to preserve the outputs 
    def setWD(self, wd):
        self.targetWD = wd

    def dockMethod(self, method):
        self.dockingMethod = method

    def fixPDB(self, checkFix):
        self.checkFix = checkFix

    def checkVals(self, verifylist):
        self.verifylist = verifylist

    def findProtein(self, folder):
        files_dict = {}  # {base_name: {'pdb': path, ...}}
        
        for path in glob(os.path.join(folder, '*')):
            if not os.path.isfile(path):
                continue  # Skip directories
            
            base_name = os.path.basename(path)  # Get filename without path
            base, ext = os.path.splitext(base_name)
            ext = ext.lower()
            
            if ext in ['.pdb', '.gpf', '.pdbqt']:
                # Initialize the base entry if not exists
                if base not in files_dict:
                    files_dict[base] = {'pdb': None, 'gpf': None, 'pdbqt': None}
                
                # Perform validations (size and content checks)
                if ext in ['.pdb', '.pdbqt']:
                    proteinSize = os.path.getsize(path) / 1024
                    if proteinSize > 20:
                        with open(path, "r") as f:
                            content = f.read()
                            if any(lig in content for lig in LIGAND_TAGS) and any(res in content for res in STD_RESIDUES):
                                self.toastSignal.emit([2000, "Found Co-crystalized ligand in the protein, Continuing..."])
                
                # Assign the path to the correct file type
                files_dict[base][ext[1:]] = path  # e.g., '.pdb' → 'pdb'
                    
        return files_dict

    def run(self):
        os.chdir(self.targetWD)
        if self.checkFix==1:
            self.toastSignal.emit([2000, "Fixing and conversion in progress..."])

        self.toastSignal.emit([2000, 
                    "Assuming that appropriate option has been selected for bounding box calculations."
                ])

        target_folders = list()
        main_dir_name = os.path.basename(os.path.normpath(self.targetWD))
        if main_dir_name != "Backup":
            target_folders.append(self.targetWD)
        
        target_folders.extend(
            entry.path 
            for entry in os.scandir(self.targetWD) 
            if entry.is_dir() and entry.name.lower() != "Backup"
        )

        valid_targets = {}

        for folder in target_folders:
            folder_files = self.findProtein(folder)

            for base, files in folder_files.items():
                # Validate based on docking method
                if self.dockingMethod == 'GPU':
                    # GPU requires at least pdb and gpf (pdbqt is optional)
                    if files['pdb'] and (files['gpf'] or files['pdbqt']):
                        valid_targets[base] = files
                    elif files['pdb'] and not files['gpf'] and not files['pdbqt']:
                        valid_targets[base] = files  # Still valid with minimal requirements
                else:
                    # Non-GPU requires at least pdb (pdbqt is optional)
                    if files['pdb']:
                        valid_targets[base] = files
                    elif (files['pdbqt'] and not files['pdb']) and self.verifylist[1]!=1:
                        valid_targets[base] = files  # Still valid with minimal requirements
                    elif (files['pdbqt'] and not files['pdb']) and self.verifylist[0] == 1:
                        error_msg = (
                            f"❌ ERROR: {base} named target protein had a PDBQT file without the corresponding PDB file. This is required for the selected site prediction method. Please verify your files and try again."
                        )
                        self.finishedSignal.emit([error_msg])
                        return

        if not valid_targets:
            self.finishedSignal.emit(["❌ ERROR: No valid targets found. Ensure each target has at least a PDB file."])
            return
        
        if valid_targets:
            # Handle site selection based on verifylist
            if any(x == 1 for x in self.verifylist):
                if self.verifylist[0] == 1:
                    self.requestResidueDialog.emit()
                    # 2) block *this* thread in its own QEventLoop until we get the text
                    loop = QEventLoop()
                    def _on_text(t):
                        self._residue_text = t
                        loop.quit()
                    # connect the one‐time slot
                    self.residuesReady.connect(_on_text)
                    loop.exec_()  

                    self.toastSignal.emit([2000, "Calculating bounding box dimensions\nusing user selected residues..."])
                    self.finishedSignal.emit(["Site", valid_targets, self._residue_text])
                else:
                    self.toastSignal.emit([
                        2000, 
                        "Calculating bounding box dimensions using IDv2's binding site predictor..."
                    ])
                    self.finishedSignal.emit(["Predict", valid_targets, "PREDICT"])
            else:
                self.toastSignal.emit([
                    2000, 
                    "Using default blind box dimensions..."
                ])
                self.finishedSignal.emit(["V1", valid_targets, "BLIND"])

class prepareTargetM2(QThread):
    finishedSignal = pyqtSignal(list)
    toastSignal = pyqtSignal(list)

    def totalFile(self, file):
        self.mode = file[0]
        self.proteins = file[1]
        self.tag = file[2]

    def IDv1Calc(self, proteinLocation):
        # Your existing file parsing to compute center and size.
        xval, yval, zval = [], [], []
        x, y, z = [], [], []
        with open(proteinLocation, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    xval.append(line[30:39])
                    yval.append(line[38:46])
                    zval.append(line[46:54])
        for i in range(len(xval)):
            x.append(float(xval[i]))
            y.append(float(yval[i]))
            z.append(float(zval[i]))

        cx = float(int((sum(x)/len(x))*1000)) / 1000.0
        lx = float(int(max(x) - min(x)) + 10)
        cy = float(int((sum(y)/len(y))*1000)) / 1000.0
        ly = float(int(max(y) - min(y)) + 10)
        cz = float(int((sum(z)/len(z))*1000)) / 1000.0
        lz = float(int(max(z) - min(z)) + 10)

        return [cx, cy, cz], [lx, ly, lz]
    
    def extract_coordinates(self, filename):
        try:
            coords = {}
            with open(filename, 'r') as file:
                for line in file:
                    # Only consider lines with an equals sign
                    if '=' in line:
                        key, value = line.split('=', 1)
                        coords[key.strip()] = value.strip()

            # Retrieve center and size values from the dictionary
            cx = coords.get('center_x')
            cy = coords.get('center_y')
            cz = coords.get('center_z')
            sx = coords.get('size_x')
            sy = coords.get('size_y')
            sz = coords.get('size_z')
            return [cx, cy, cz, sx, sy, sz]
        except:
            return ["❌ No suitable configuration file was found. Please review your files and try again."]

    def run(self):
        if self.mode == 'Site':
            proteinSites = dict()
            specificResidues = self.tag.replace(' ', '').split(',')
            resDict = ID_utils.separateMultiResults(specificResidues)
            # Get results as {"Hemoglobin": {"A": (256, 320), "B": (315, 330)}, "Cytochrome": {"B": (295, 300)}

            result = ID_utils.getMaxMinMulti(resDict)
            # Gives {'Hemoglobin': ((256, 'A'), (330, 'B')), 'Cytochrome': ((295, 'B'), (300, 'B'))}

            for target_name, files in self.proteins.items():
                pdb = files.get('pdb')
                pdbqt = files.get('pdbqt')

                for protein, mm in result.items():
                    if target_name.lower() == protein.lower():
                        minRes, maxRes = mm  # Each is a tuple: (residue, chain)
                        
                        # Choose pdbqt if available and pdb is not, otherwise use pdb
                        if pdbqt is not None and pdb is None:
                            first_coords = ID_utils.get_residue_coords(pdbqt, minRes[1], minRes[0])
                            last_coords  = ID_utils.get_residue_coords(pdbqt, maxRes[1], maxRes[0])
                        else:
                            first_coords = ID_utils.get_residue_coords(pdb, minRes[1], minRes[0])
                            last_coords  = ID_utils.get_residue_coords(pdb, maxRes[1], maxRes[0])
                        
                        center, size = ID_utils.calculate_bounding_box(first_coords, last_coords)

                proteinSites[target_name.lower()] = [files, center, size]

            self.finishedSignal.emit([proteinSites])

        elif self.mode == 'Predict':
            proteinSites = dict()
            pdbs = [files.get('pdb') in target_name, files in self.proteins.items()]
            for pdb in pdbs:
                pdb_path = pdb
                try:
                    kv_res = pyKVFinder.run_workflow(pdb_path)
                except Exception:
                    self.finishedSignal.emit(["❌ No cavities found by IDv2!"])
                    return

                vols = kv_res.volume

                cid = max(vols, key=vols.get)
                unfiltered = kv_res.residues.get(cid, [])
                if not unfiltered:
                    self.finishedSignal.emit(["❌ No residues for largest cavity."])
                    return
                
                res_list = [entry for entry in unfiltered if entry[2] in STD_RESIDUES2]
                f_res, l_res = res_list[0], res_list[-1]
                c1 = ID_utils.get_residue_coords(pdb_path, f_res[1], int(f_res[0]))
                c2 = ID_utils.get_residue_coords(pdb_path, l_res[1], int(l_res[0]))
                center, size = ID_utils.calculate_bounding_box(c1, c2)

                proteinSites[protein.split('.')[0].lower()] = [files, center, size]

            self.finishedSignal.emit([proteinSites])
            return
                
        else:
            proteinSites = dict()
            for target_name, files in self.proteins.items():
                pdb = files.get('pdb')
                pdbqt = files.get('pdbqt')
                # Fallback using default (blind) box dimensions.
                if pdbqt and not pdb:
                    center, size = self.IDv1Calc(pdbqt)
                else:
                    center, size = self.IDv1Calc(pdb)

                proteinSites[target_name.lower()] = [files, center, size]
            
            self.finishedSignal.emit([proteinSites])

class prepareTargetM3(QThread):
    finishedSignal = pyqtSignal(str)
    toastSignal = pyqtSignal(list)

    def setTargetWD(self, wd):
        self.targetWD = wd

    def allThings(self, text):
        self.resDict = text[0]
        #Dictionary contains target: [filesDictionary, center, size]

    def fixPDB(self, checkFix):
        self.checkFix = checkFix

    def dockMethod(self, method):
        self.method = method

    def fixPdb(self, proteinFile, targetWD):
        try:
            # Prepare a temporary filename for the fixed version.
            base_name = os.path.splitext(os.path.basename(proteinFile))[0]
            temp_fixed_file = os.path.join(os.path.dirname(proteinFile), f"{base_name}_temp_fixed.pdb")
            
            # Fix the PDB file using PDBFixer.
            fixer = PDBFixer(filename=proteinFile)
            fixer.findMissingResidues()
            fixer.findNonstandardResidues()
            fixer.replaceNonstandardResidues()
            fixer.removeHeterogens(True)
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            
            # Write the fixed structure to a temporary file.
            with open(temp_fixed_file, 'w') as f:
                PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
            
            # Back up the original file into the main Backup directory.
            backup_dir = os.path.join(targetWD, "Backup")
            if not os.path.exists(backup_dir):
                os.makedirs(backup_dir)

            # Create the backup file path and move the original there.
            backup_file = os.path.join(backup_dir, f"{base_name}_original.pdb")
            shutil.move(proteinFile, backup_file)
            
            # Replace the original file with the fixed temporary file.
            shutil.move(temp_fixed_file, proteinFile)
            
            return "Done"
        except Exception:
            return f"❌ Unable to fix the PDB, please review your files and try again."
        
    def IDv1Calc(self, proteinLocation):
        # Your existing file parsing to compute center and size.
        xval, yval, zval = [], [], []
        x, y, z = [], [], []
        with open(proteinLocation, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    xval.append(line[30:39])
                    yval.append(line[38:46])
                    zval.append(line[46:54])
        for i in range(len(xval)):
            x.append(float(xval[i]))
            y.append(float(yval[i]))
            z.append(float(zval[i]))

        cx = float(int((sum(x)/len(x))*1000)) / 1000.0
        lx = float(int(max(x) - min(x)) + 10)
        cy = float(int((sum(y)/len(y))*1000)) / 1000.0
        ly = float(int(max(y) - min(y)) + 10)
        cz = float(int((sum(z)/len(z))*1000)) / 1000.0
        lz = float(int(max(z) - min(z)) + 10)

        return [cx, cy, cz], [lx, ly, lz]
    
    def find_file_in_dir(self, filename, search_dir):
        if filename is None:
            return None
        basename = os.path.basename(filename)
        pattern = os.path.join(search_dir, basename)
        return next(iter(glob(pattern)), None)
    
    def convertPDB(self, pdbFile, center, size):
        pdbqtFile = pdbFile.replace(".pdb", ".pdbqt")
        gpfFile = pdbFile.replace(".pdb", ".gpf")
        if self.method == "GPU":
            ar = [
                "--read_pdb", pdbqtFile,
                "-g", gpfFile, "--write_pdbqt", pdbqtFile,
                "--box_center", str(center[0]), str(center[1]), str(center[2]),
                "--box_size", str(int(size[0])), str(int(size[1])), str(int(size[2]))
            ]
        else:
            ar = [
                "--read_pdb", pdbFile,
                "-p", pdbqtFile,
                "--box_center", str(center[0]), str(center[1]), str(center[2]),
                "--box_size", str(int(size[0])), str(int(size[1])), str(int(size[2]))
            ]

        # Validate executable path
        primary_exe = resource_path("mk_prepare_receptor.exe")
        if not os.path.exists(primary_exe):
            return 0

        # Run mk_prepare_receptor
        self.process.start(primary_exe, ar)
        self.process.waitForFinished(-1)
        output = self.process.readAllStandardOutput().data().decode()
        name = (
            os.path.splitext(os.path.basename(pdbFile))[0] 
            if pdbFile 
            else os.path.splitext(os.path.basename(pdbqtFile))[0]
        )

        if "error" in output.lower() or "failed" in output.lower():
            if self.method == "GPU":
                return "❌ Could not prepare the receptor molecule for the selected docking method"
            # Fallback to OpenBabel
            inputform = '-ipdb'
            alt_exe = "obabel"
            ar = [inputform, pdbFile, "-opdbqt", "-O", pdbqtFile, "-xr", "-xc", "-xn", "--partialcharge", "gasteiger"]
            
            self.process.start(alt_exe, ar)
            self.process.waitForFinished(-1)
            output2 = self.process.readAllStandardOutput().data().decode()
            
            if "error" in output2.lower():
                return 0
            else:
                return 1
        else:
            return 1
        
    def run(self):
        totalCount = 0
        numTargets = len(self.resDict.keys())
        for pdbName, (data_dict, center, size) in self.resDict.items():
            pdb_pattern = data_dict.get("pdb")
            tPDB = next(iter(glob(pdb_pattern)), None) if pdb_pattern else None
            pdbLoc = os.path.abspath(tPDB) if tPDB is not None else None
            pdbqt_pattern = data_dict.get("pdbqt")
            tPDBQT = next(iter(glob(pdbqt_pattern)), None) if pdbqt_pattern else None
            pdbqtLoc = os.path.abspath(tPDBQT) if tPDBQT is not None else None
            gpf_pattern = data_dict.get("gpf")
            tGPF = next(iter(glob(gpf_pattern)), None) if gpf_pattern else None
            gpfLoc = os.path.abspath(tGPF) if tGPF is not None else None

            # Use targetWD as the base directory for protein_dir
            protein_dir = os.path.join(self.targetWD, pdbName)  # FIXED HERE

            if pdbLoc and (not pdbqtLoc and not gpfLoc):
                with open(pdbLoc, 'r') as f:
                    lines = f.readlines()

                fixedPDB = []
                for line in lines:
                    if not line.startswith(('ATOM  ', 'HETATM')):
                        continue
                    resname = line[17:20].strip()
                    if resname == 'HOH':
                        continue
                    if resname not in STD_RESIDUES2:
                        continue
                    
                    fixedPDB.append(line)

                # Write back out
                with open(pdbLoc, 'w') as f:
                    f.writelines(fixedPDB)

                def find_altlocs(pdb_path):
                    """
                    Scan a PDB file and report any altLoc codes found.
                    """
                    altlocs = set()
                    with open(pdb_path) as f:
                        for line in f:
                            if line.startswith(("ATOM  ", "HETATM")):
                                code = line[16]              # column 17 in 0-based string indexing
                                if code not in (" ", ""):
                                    altlocs.add(code)
                    return sorted(altlocs)


                def strip_altlocs(pdb_path, output_path, keep="A"):
                    """
                    Write a new PDB with only the desired altLoc record (default “A”),
                    or all atoms that have no altLoc.
                    """
                    with open(pdb_path) as src, open(output_path, "w") as dst:
                        for line in src:
                            if line.startswith(("ATOM  ", "HETATM")):
                                code = line[16]
                                # keep if either no altLoc or exactly the one we want
                                if code in (" ", "") or code == keep:
                                    # blank out column 17 so downstream tools see no alternate locs
                                    dst.write(line[:16] + " " + line[17:])
                            else:
                                dst.write(line)

                pdb_path = pdbLoc
                altlocs = find_altlocs(pdb_path)
                if altlocs:
                    dirn, base = os.path.split(pdb_path)
                    fd, tmp_path = tempfile.mkstemp(prefix=base, dir=dirn, text=True)
                    os.close(fd)
                    strip_altlocs(pdb_path, tmp_path, keep=altlocs[0])
                    os.replace(tmp_path, pdb_path)

                self.process = QProcess()
                self.process.setProcessChannelMode(QProcess.MergedChannels)
                primary_exe = resource_path("mk_prepare_receptor.exe")
                minimal_args = ["--read_pdb", pdbLoc]
                self.process.start(primary_exe, minimal_args)
                self.process.waitForFinished(-1)
                output = self.process.readAllStandardOutput().data().decode()

                if not os.path.exists(protein_dir):
                    os.makedirs(protein_dir)
                    shutil.move(pdbLoc, protein_dir)

                pdbLoc = self.find_file_in_dir(pdbLoc, protein_dir)

                if ("error" in output.lower() or "failed" in output.lower()):
                    if self.method == "GPU":
                        self.toastSignal.emit([2000, "Fixing and conversion in progress..."])
                        testTxt = self.fixPdb(pdbLoc, self.targetWD)
                    elif self.checkFix==1:
                        self.toastSignal.emit([2000, "Fixing and conversion in progress..."])
                        testTxt = self.fixPdb(pdbLoc, self.targetWD)
                    else:
                        self.toastSignal.emit([2000, f"Conversion in progress..."])
                        testTxt = "DoneWithoutFix"
                else:
                    self.toastSignal.emit([2000, "Conversion in progress..."])
                    testTxt = "DoneWithoutFix"

                if testTxt == "Done":
                    testConv = self.convertPDB(pdbLoc, center, size)
                    totalCount += testConv

                elif testTxt == "DoneWithoutFix":
                    testConv = self.convertPDB(pdbLoc, center, size)
                    totalCount += testConv

            elif pdbqtLoc and (not pdbLoc and not gpfLoc):
                if not os.path.exists(protein_dir):
                    os.makedirs(protein_dir)
                    shutil.move(pdbqtLoc, protein_dir)
                totalCount += 1

            elif pdbLoc and pdbqtLoc and not gpfLoc:
                confTest = glob(os.path.join(os.path.dirname(pdbLoc), "*conf.txt"))
                confTest = confTest[0] if confTest else None
                if not os.path.exists(protein_dir):
                    os.makedirs(protein_dir)
                    shutil.move(pdbLoc, protein_dir)
                    shutil.move(pdbqtLoc, protein_dir)
                    if confTest:
                        shutil.move(confTest, protein_dir)
                totalCount += 1
            elif pdbLoc and gpfLoc and pdbqtLoc:
                if not os.path.exists(protein_dir):
                    os.makedirs(protein_dir)
                    shutil.move(pdbLoc, protein_dir)
                    shutil.move(gpfLoc, protein_dir)
                    shutil.move(pdbqtLoc, protein_dir)
                totalCount += 1
            elif pdbLoc and gpfLoc:
                if not os.path.exists(protein_dir):
                    os.makedirs(protein_dir)
                    shutil.move(pdbLoc, protein_dir)
                    shutil.move(gpfLoc, protein_dir)
                totalCount += 1

        if totalCount == numTargets:
            self.finishedSignal.emit("✅ Finished the target conversion process")
        else:           
            self.finishedSignal.emit("❌ Failed the target conversion process")

class prepareLigand(QThread):
    finishedSignal = pyqtSignal(str)
    
    def setLigandWD(self, wd):
        self.wd = wd

    def extractLibs(self, libraries):
        libmulti = None
        libsingle = None
        cwd = os.getcwd()
        ligCount = 0
        if len(libraries) > 1:
            libmulti = os.path.join(cwd, "Ligand_libraries_found")
            lib = libmulti
            if not os.path.exists(libmulti):
                os.makedirs(libmulti)
            for library in libraries:
                if library.endswith(".mol2"):
                    f = open(library, mode='r', encoding='utf-8')
                    text = f.read()
                    f.close()
                    d = "@<TRIPOS>MOLECULE"
                    splitted = [d+e for e in text.split(d) if e]
                    for i in splitted:
                        ligCount += 1
                        test = i.splitlines()
                        name = test[1]
                        outputfile =  name+".mol2"
                        u = open(outputfile,'wb')
                        file = "\n".join(test).decode('utf-8')
                        u.write(file)
                        u.close()
                elif library.endswith(".sdf"):
                    f = open(library, mode='r', encoding='utf-8')
                    text = f.read()
                    arr = text.split("$$$$")[:-1]
                    for file in arr:
                        ligCount += 1
                        if file == arr[0]:
                            name = file.splitlines()[0]
                        elif file == arr[len(arr)-1]:
                            name = file.splitlines()[1]
                        else:
                            name = file.splitlines()[1]
                        u = open(os.path.join(cwd, name+".sdf"),'wb')
                        file = "".join(file[1:len(file)-1])
                        u.write(file.encode('utf-8'))
                        u.close()
                original_library_path = os.path.join(cwd, library)
                target_path = os.path.join(libmulti, os.path.basename(library))
                shutil.move(original_library_path, target_path)

        elif len(libraries) == 1:
            library = libraries[0]
            libsingle = os.path.join(cwd, "Ligand_library_found")
            lib = libsingle
            if not os.path.exists(libsingle):
                os.makedirs(libsingle)
            if libraries[0].endswith(".mol2"):
                with open(libraries[0], mode='r', encoding='utf-8') as f:
                    text = f.read()

                d = "@<TRIPOS>MOLECULE"
                splitted = [d+e for e in text.split(d) if e]
                for i in splitted:
                    ligCount += 1
                    test = i.splitlines()
                    name = test[1]
                    outputfile =  name+".mol2"
                    with open(outputfile,'wb') as u:
                        file = "\n".join(test).decode('utf-8')
                        u.write(file)

            elif libraries[0].endswith(".sdf"):
                with open(libraries[0], mode='r', encoding='utf-8') as f:
                    text = f.read()
                arr = text.split("$$$$")[:-1]
                for file in arr:
                    ligCount += 1
                    if file == arr[0]:
                        name = file.splitlines()[0]
                    elif file == arr[len(arr)-1]:
                        name = file.splitlines()[1]
                    else:
                        name = file.splitlines()[1]
                    with open(os.path.join(cwd, name+".sdf"),'wb') as u:
                        file = "".join(file[1:len(file)-1])
                        u.write(file.encode('utf-8'))

            original_library_path = os.path.join(cwd, library)
            target_path = os.path.join(libsingle, os.path.basename(library))

            shutil.move(original_library_path, target_path)

        return lib
    
    def findLigs(self):       
        ligandFiles = list()
        libraryFiles = list()
        potentialLigands = glob('*.pdb')+ glob('*.mol2') + glob('*.sdf') + glob('*.mol') + glob('*.pdbqt')
        for ligFile in potentialLigands:
            with open(ligFile, "r") as f:
                fileContent = f.read()
            if any(res in fileContent for res in STD_RESIDUES) or "BOX" in fileContent:
                pass
            else:
                ext = os.path.splitext(ligFile)[1].lower()
                if ext == '.sdf':
                    # SDF multi‐molecule delimiter is $$$$
                    if fileContent.count('$$$$') > 1:
                        libraryFiles.append(ligFile)
                    else:
                        ligandFiles.append(ligFile)

                elif ext == '.mol2':
                    # MOL2 multi‐molecule marker is @<TRIPOS>MOLECULE
                    if fileContent.count('@<TRIPOS>MOLECULE') > 1:
                        libraryFiles.append(ligFile)
                    else:
                        ligandFiles.append(ligFile)

                else:
                    ligandFiles.append(ligFile)

        return ligandFiles, libraryFiles
    
    def findPDBQTLigs(self):       
        ligandFiles = list()       
        potentialLigands = glob('*.pdbqt')
        for ligFile in potentialLigands:
            ligandFiles.append(ligFile)

        return ligandFiles

    def convertLig(self, ligFile):
        # Ensure the ligand file is an absolute path.
        ligFile = os.path.abspath(ligFile)
        
        if ligFile.endswith(".pdbqt"):
            pass
        else:       
            # Build an absolute path to the executable.
            cd = "obabel"
            format = ligFile.split('.')[1]
            inputform = "-i"+format
            tem = re.sub('\..*', '.pdbqt', ligFile)
            lig = tem
            ar = [inputform, ligFile, "-opdbqt", "-O", lig, "--partialcharge gasteiger", "-h"]
            
            self.process = QProcess()
            self.process.execute(cd, ar)
            self.process.waitForFinished(-1)

    def run(self):
        os.chdir(self.wd)
        ligandFiles, libraryFiles = self.findLigs()
        if libraryFiles:
            self.extractLibs(libraryFiles)
        ligandFiles, libraryFiles = self.findLigs()

        if ligandFiles:
            for ligandFile in ligandFiles:
                self.convertLig(ligandFile)
        else:            
            self.finishedSignal.emit("❌ No ligands found! Please check your files and try again.")
            return

        ligandFiles = self.findPDBQTLigs()

        if len(ligandFiles) > 1:
            self.finishedSignal.emit(f"✅ {len(ligandFiles)} candidates ready for docking analysis!")
        else:
            self.finishedSignal.emit(f"✅ Ligand {ligandFiles[0]} is ready for docking analysis.")
        
class prepareConfigM1(QThread):
    requestResidueDialog = pyqtSignal()
    residuesReady = pyqtSignal(str)
    toastSignal = pyqtSignal(list)
    finishedSignal = pyqtSignal(list)
    
    def setTargetWD(self, wd):
        self.targetWD = wd
    
    def dockMethod(self, text):
        self.method = text

    def checkVals(self, verifylist):
        self.verifylist = verifylist

    def findProtein(self, folder):
        files_dict = {}  # {base_name: {'pdb': path, ...}}
        
        for path in glob(os.path.join(folder, '*')):
            if not os.path.isfile(path):
                continue  # Skip directories
            
            base_name = os.path.basename(path)  # Get filename without path
            base, ext = os.path.splitext(base_name)
            ext = ext.lower()
            
            if ext in ['.pdb', '.pdbqt']:
                # Initialize the base entry if not exists
                if base not in files_dict:
                    files_dict[base] = {'pdb': None, 'pdbqt': None}
                
                # Perform validations (size and content checks)
                if ext in ['.pdb', '.pdbqt']:
                    proteinSize = os.path.getsize(path) / 1024
                    if proteinSize > 20:
                        with open(path, "r") as f:
                            content = f.read()
                            if any(lig in content for lig in LIGAND_TAGS) and any(res in content for res in STD_RESIDUES):
                                self.toastSignal.emit([2000, "Found Co-crystalized ligand in the protein, Continuing..."])
                
                # Assign the path to the correct file type
                files_dict[base][ext[1:]] = path  # e.g., '.pdb' → 'pdb'
                    
        return files_dict
    
    def run(self):
        os.chdir(self.targetWD)
        if len(glob("conf.txt")) != 0:
            self.finishedSignal.emit(["⚠️ IDv2 detected a configuration file. If you want to create a new configuration file, please remove the existing one and try again."])
            return
        
        if self.method == "AutoDock-GPU":
            self.finishedSignal.emit(
                ["⚠️ AutoDock-GPU was selected as the docking method; therefore, no configuration file will be generated."]
            )
        else:
            self.toastSignal.emit([2000, "Assuming that the receptor files have already been prepared."])
            target_folders = list()
            main_dir_name = os.path.basename(os.path.normpath(self.targetWD)).lower()
            if main_dir_name != "backup":
                target_folders.append(self.targetWD)
            
            target_folders.extend(
                entry.path 
                for entry in os.scandir(self.targetWD) 
                if entry.is_dir() and entry.name.lower() != "backup"
            )

            valid_targets = {}

            for folder in target_folders:
                folder_files = self.findProtein(folder)
                for base, files in folder_files.items():
                    if files['pdb']:
                            valid_targets[base] = files
                    elif (files['pdbqt'] and not files['pdb']) and self.verifylist[1]!=1:
                        valid_targets[base] = files  # Still valid with minimal requirements
                    elif (files['pdbqt'] and not files['pdb']) and self.verifylist[0] == 1:
                        error_msg = (
                            f"❌ ERROR: {base} named target protein had a PDBQT file without the corresponding PDB file. This is required for the selected site prediction method. Please verify your files and try again."
                        )
                        
                        self.finishedSignal.emit([error_msg])
                        return
        
        if not valid_targets:            
            self.finishedSignal.emit(["❌ ERROR: No valid targets found. Ensure each target has at least a PDB file."])
            return
        
        if valid_targets:
            # Handle site selection based on verifylist
            if any(x == 1 for x in self.verifylist):
                if self.verifylist[0] == 1:
                    self.requestResidueDialog.emit()
                    # 2) block *this* thread in its own QEventLoop until we get the text
                    loop = QEventLoop()
                    def _on_text(t):
                        self._residue_text = t
                        loop.quit()
                    # connect the one‐time slot
                    self.residuesReady.connect(_on_text)
                    loop.exec_()  
                    self.toastSignal.emit([2000, "Calculating bounding box dimensions\nusing user selected residues..."])
                    self.finishedSignal.emit(["Site", valid_targets, self._residue_text])
                else:
                    self.toastSignal.emit([
                        2000,
                        "Calculating bounding box dimensions using IDv2's binding site predictor..."])
                    
                    self.finishedSignal.emit(["Predict", valid_targets, "PREDICT"])
            else:
                self.toastSignal.emit([
                    2000,
                    "Using default blind box dimensions..."
                ])  
                self.finishedSignal.emit(["V1", valid_targets, "BLIND"])

class prepareConfigM2(QThread):
    finishedSignal = pyqtSignal(list)

    def totalFile(self, file):
        self.mode = file[0]
        self.proteins = file[1]
        self.tag = file[2]

    def IDv1Calc(self, proteinLocation):
        # Your existing file parsing to compute center and size.
        xval, yval, zval = [], [], []
        x, y, z = [], [], []
        with open(proteinLocation, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    xval.append(line[30:39])
                    yval.append(line[38:46])
                    zval.append(line[46:54])
        for i in range(len(xval)):
            x.append(float(xval[i]))
            y.append(float(yval[i]))
            z.append(float(zval[i]))

        cx = float(int((sum(x)/len(x))*1000)) / 1000.0
        lx = float(int(max(x) - min(x)) + 10)
        cy = float(int((sum(y)/len(y))*1000)) / 1000.0
        ly = float(int(max(y) - min(y)) + 10)
        cz = float(int((sum(z)/len(z))*1000)) / 1000.0
        lz = float(int(max(z) - min(z)) + 10)

        return [cx, cy, cz], [lx, ly, lz]
    
    def extract_coordinates(self, filename):
        try:
            coords = {}
            with open(filename, 'r') as file:
                for line in file:
                    # Only consider lines with an equals sign
                    if '=' in line:
                        key, value = line.split('=', 1)
                        coords[key.strip()] = value.strip()

            # Retrieve center and size values from the dictionary
            cx = coords.get('center_x')
            cy = coords.get('center_y')
            cz = coords.get('center_z')
            sx = coords.get('size_x')
            sy = coords.get('size_y')
            sz = coords.get('size_z')
            return [cx, cy, cz, sx, sy, sz]
        except:
            return ["❌ No suitable configuration file was found. Please review your files and try again."]

    def run(self):
        if self.mode == 'Site':
            proteinSites = dict()
            specificResidues = self.tag.replace(' ', '').split(',')
            resDict = ID_utils.separateMultiResults(specificResidues)
            # Get results as {"Hemoglobin": {"A": (256, 320), "B": (315, 330)}, "Cytochrome": {"B": (295, 300)}

            result = ID_utils.getMaxMinMulti(resDict)
            # Gives {'Hemoglobin': ((256, 'A'), (330, 'B')), 'Cytochrome': ((295, 'B'), (300, 'B'))}

            for target_name, files in self.proteins.items():
                pdb = files.get('pdb')
                pdbqt = files.get('pdbqt')

                for protein, mm in result.items():
                    if target_name.lower() == protein.lower():
                        minRes, maxRes = mm  # Each is a tuple: (residue, chain)
                        
                        # Choose pdbqt if available and pdb is not, otherwise use pdb
                        if pdbqt is not None and pdb is None:
                            first_coords = ID_utils.get_residue_coords(pdbqt, minRes[1], minRes[0])
                            last_coords  = ID_utils.get_residue_coords(pdbqt, maxRes[1], maxRes[0])
                        else:
                            first_coords = ID_utils.get_residue_coords(pdb, minRes[1], minRes[0])
                            last_coords  = ID_utils.get_residue_coords(pdb, maxRes[1], maxRes[0])
                        
                        center, size = ID_utils.calculate_bounding_box(first_coords, last_coords)

                proteinSites[target_name.lower()] = [files, center, size]
            
            self.finishedSignal.emit([proteinSites])

        elif self.mode == 'Predict':
            proteinSites = dict()
            pdbs = [files.get('pdb') in target_name, files in self.proteins.items()]
            for pdb in pdbs:
                pdb_path = pdb
                try:
                    kv_res = pyKVFinder.run_workflow(pdb_path)
                except Exception:
                    self.finishedSignal.emit(["❌ No cavities found by IDv2!"])
                    return

                vols = kv_res.volume

                cid = max(vols, key=vols.get)
                unfiltered = kv_res.residues.get(cid, [])
                if not unfiltered:
                    self.finishedSignal.emit(["❌ No residues for largest cavity."])
                    return
                
                res_list = [entry for entry in unfiltered if entry[2] in STD_RESIDUES2]
                f_res, l_res = res_list[0], res_list[-1]
                c1 = ID_utils.get_residue_coords(pdb_path, f_res[1], int(f_res[0]))
                c2 = ID_utils.get_residue_coords(pdb_path, l_res[1], int(l_res[0]))
                center, size = ID_utils.calculate_bounding_box(c1, c2)

                proteinSites[protein.split('.')[0].lower()] = [files, center, size]

            self.finishedSignal.emit([proteinSites])
            return
        
        else:
            proteinSites = dict()
            for target_name, files in self.proteins.items():
                pdb = files.get('pdb')
                pdbqt = files.get('pdbqt')
                # Fallback using default (blind) box dimensions.
                if pdbqt and not pdb:
                    center, size = self.IDv1Calc(pdbqt)
                else:
                    center, size = self.IDv1Calc(pdb)

                proteinSites[target_name.lower()] = [files, center, size]
            
            self.finishedSignal.emit([proteinSites])

class prepareConfigM3(QThread):
    finishedSignal = pyqtSignal(str)

    def setTargetWD(self, wd):
        self.targetWD = wd

    def allThings(self, text):
        self.resDict = text[0]
    
    def dockMethod(self, method):
        self.method = method

    def createConfig(self, proteinName, pdbqtFile, center, size):
        # Write the configuration file in the required format.
        with open(os.path.join(proteinName, f"{proteinName.lower()}_conf.txt"), 'w') as f:
            f.write(f"receptor = {pdbqtFile}\n\n")
            f.write(f"center_x = {center[0]:.2f}\n")
            f.write(f"center_y = {center[1]:.2f}\n")
            f.write(f"center_z = {center[2]:.2f}\n\n")
            f.write(f"size_x = {int(size[0])}\n")
            f.write(f"size_y = {int(size[1])}\n")
            f.write(f"size_z = {int(size[2])}\n")

    def run(self):
        totalCount = 0
        numTargets = len(self.resDict.keys())
        for pdbName, (data_dict, center, size) in self.resDict.items():
            pdbqt_pattern = data_dict.get("pdbqt")
            tPDBQT = next(iter(glob(pdbqt_pattern)), None) if pdbqt_pattern else None
            pdbqtLoc = os.path.abspath(tPDBQT) if tPDBQT is not None else None

            protein_dir = os.path.join(self.targetWD, pdbName)

            if pdbqtLoc:
                if not os.path.exists(protein_dir):
                    os.makedirs(protein_dir)
                name = os.path.splitext(os.path.basename(pdbqtLoc))[0]
                self.createConfig(name, pdbqtLoc, center, size)
                totalCount += 1

        if totalCount == numTargets:
            self.finishedSignal.emit(f"✅ {totalCount} configuration files successfully prepared!")

class convertFiles(QThread):
    toastSignal = pyqtSignal(list)
    finishedSignal = pyqtSignal(str)

    def setTargetWD(self, wd):
        self.targetWD = wd
        
    def setLigandWD(self, wd):
        self.ligandWD = wd

    def run(self):
        # --- Step 1: Validate directories ---
        if not self.targetWD:
            self.finishedSignal.emit("❌ Target directory not set.")
            return
        if not self.ligandWD:
            self.finishedSignal.emit("❌ Ligand directory not set.")
            return

        obabel_exe = "obabel"

        # --- Step 2: Convert receptor files ---
        try:
            os.chdir(self.targetWD)
        except Exception:
            self.finishedSignal.emit(f"❌ Cannot access target directory '{self.targetWD}'")
            return

        receptor_patterns = ['*.pdb', '*.gpf', '*.pdbqt']
        receptor_files = []
        for pat in receptor_patterns:
            receptor_files.extend(glob(pat))

        if not receptor_files:
            pass
        else:
            # Group by base name
            groups = {}
            for path in receptor_files:
                if not os.path.isfile(path):
                    continue
                base, ext = os.path.splitext(path)
                ext = ext.lower().lstrip('.')
                groups.setdefault(base, set()).add(ext)

            for base, exts in groups.items():
                if 'pdbqt' in exts:
                    continue  # already converted

                # choose input file: prefer .pdb over .gpf
                if 'pdb' in exts:
                    src = f"{base}.pdb"
                elif 'gpf' in exts:
                    src = f"{base}.gpf"
                else:
                    continue
                out = f"{base}.pdbqt"

                # Size & content checks
                try:
                    size_kb = os.path.getsize(src) / 1024
                    with open(src, 'r') as f:
                        content = f.read()
                except Exception:
                    continue

                # Notify co‑crystallized ligand
                if any(lig in content for lig in LIGAND_TAGS) and any(res in content for res in STD_RESIDUES):
                    self.toastSignal.emit([2000, f"🔍 Found co‑crystallized ligand in '{src}'."])

                # Build OpenBabel args for protein
                args = [
                    f"-i{ext}", src,
                    "-opdbqt", "-O", out,
                    "-xr", "-xn", "-xc", "-xp",
                    "--partialcharge", "gasteiger"
                ]
                QProcess.execute(obabel_exe, args)

        # --- Step 3: Convert ligand files ---
        try:
            os.chdir(self.ligandWD)
        except Exception:
            self.finishedSignal.emit(f"❌ Cannot access ligand directory '{self.ligandWD}'")
            return

        ligand_patterns = ['*.mol2', '*.sdf', '*.mol', '*.pdb', '*.pdbqt']
        ligand_files = []
        for pat in ligand_patterns:
            ligand_files.extend(glob(pat))

        if not ligand_files:
            pass
        else:
            for src in ligand_files:
                base, ext = os.path.splitext(src)
                ext = ext.lower().lstrip('.')
                out = f"{base}.pdbqt"

                # Skip if already converted
                if os.path.exists(out):
                    continue

                # Read & size check
                try:
                    size_kb = os.path.getsize(src) / 1024
                    with open(src, 'r') as f:
                        content = f.read()
                except Exception:
                    continue

                # --- library‐file checks for .sdf and .mol2 ---
                if ext == 'sdf':
                    count = content.count('$$$$')
                    if count > 1:
                        continue
                elif ext == 'mol2':
                    count = len(re.findall(r'@<TRIPOS>MOLECULE', content))
                    if count > 1:
                        continue

                # existing library logic
                is_library = (size_kb > 20 and
                              (any(res in content for res in STD_RESIDUES) or "BOX" in content))
                if is_library:
                    continue  # skip true library files

                # skip too‐small non‐PDB ligands
                if size_kb <= 20 and ext != 'pdb':
                    continue

                # Build OpenBabel args for ligand
                args = [
                    f"-i{ext}", src,
                    "-opdbqt", "-O", out,
                    "--partialcharge", "gasteiger",
                    "--DelNonPolarH"
                ]
                code = QProcess.execute(obabel_exe, args)

        # --- Step 4: Done ---
        self.finishedSignal.emit("✅ Conversion complete: all missing .pdbqt files generated.")


