import os
import re
import sys
import shutil
import ID_utils
import tempfile
import pyKVFinder
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

def inspectProtein(wd, method):
    """
    Inspects PDB files in the given directory and returns a plain-text
    message (and title) summarizing any issues.

    If method='GPU', attempts GPU processing and reports if it fails.

    Returns:
      msg (str): multiline text summary suitable for QMessageBox.setDetailedText
      title (str): dialog title
    """

    # 1. Find all PDB files in the directory.
    pattern = os.path.join(wd, "*.pdb")

    # Find all matching files
    pdb_files = glob(pattern)

    # Pick the one with the maximum file size
    pdbfile = max(pdb_files, key=os.path.getsize)

    if not pdbfile:
        return "No PDB file found!", "Error"

    # 2. Select a target protein
    target_protein = None
    has_cocrystal_ligand = False
    with open(pdbfile, 'r') as f:
        content = f.read()
    ligand_found = any(tag in content for tag in LIGAND_TAGS)
    std_res_found = any(res in content for res in STD_RESIDUES)
    if ligand_found and std_res_found:
        target_protein = pdbfile
        has_cocrystal_ligand = True
    elif ligand_found and not std_res_found:
        return "No PDB file found!", "Error"
    else:
        target_protein = pdbfile
        
    filename = os.path.basename(target_protein)

    # GPU-specific processing
    if method.upper() == 'GPU':
        process = QProcess()
        process.setProcessChannelMode(QProcess.MergedChannels)
        primary_exe = resource_path("mk_prepare_receptor.exe")
        args = ["--read_pdb", target_protein]  # add GPU flag if needed
        process.start(primary_exe, args)
        process.waitForFinished(-1)
        output = process.readAllStandardOutput().data().decode()

        # Check for GPU processing errors
        if 'error' in output.lower() or 'failed' in output.lower():
            msg = f"{filename} cannot be processed using GPU methods. Please try again with CPU methods."
            return msg, "GPU Error"
        # Otherwise, continue with analysis on PDB file

    # 3. Read lines
    with open(target_protein, 'r') as f:
        lines = f.readlines()

    atomnos = []
    chains = set()
    hetcount = 0
    anomalies = []
    altloc_found = False
    insertion_found = False
    hydrogen_count = 0

    for i, line in enumerate(lines):
        record = line[0:6].strip()
        if record in ("ATOM", "HETATM"):
            # occupancy and temp factor
            occ = line[54:60].strip()
            temp = line[60:66].strip()
            if occ and not _is_float(occ):
                anomalies.append(f"Line {i+1}: invalid occupancy '{occ}'")
            if temp and not _is_float(temp):
                anomalies.append(f"Line {i+1}: invalid temperature factor '{temp}'")

            resname = line[17:20].strip()
            if len(resname) != 3:
                anomalies.append(f"Line {i+1}: unusual residue name '{resname}'")

            alt = line[16].strip()
            if alt and alt != "A":
                altloc_found = True

            ins = line[26].strip()
            if ins:
                insertion_found = True

            if record == "ATOM":
                num = line[22:26].strip()
                if num:
                    if num.isdigit():
                        atomnos.append(int(num))
                    else:
                        anomalies.append(f"Line {i+1}: invalid residue number '{num}'")
                cid = line[21].strip()
                if cid:
                    chains.add(cid)
                else:
                    anomalies.append(f"Line {i+1}: missing chain ID")

                name = line[12:16].strip()
                if name.startswith("H"):
                    hydrogen_count += 1
            else:
                hetcount += 1

    # check sequence gaps
    missing = []
    if atomnos:
        seq = sorted(set(atomnos))
        missing = [str(n) for n in range(seq[0], seq[-1]+1) if n not in seq]

    # build issues
    issues = []
    if has_cocrystal_ligand:
        issues.append("Found co-crystallized ligand")
    if missing:
        issues.append(f"Missing residues: {len(missing)} gaps ({', '.join(missing)})")
    if len(chains) > 1:
        issues.append(f"Multiple chains: {', '.join(sorted(chains))}")
    if hetcount:
        issues.append(f"HETATM records: {hetcount}")
    if altloc_found:
        issues.append("Alternate location indicators detected")
    if insertion_found:
        issues.append("Insertion codes detected")
    if hydrogen_count:
        issues.append(f"Hydrogen atoms: {hydrogen_count}")

    # format output
    lines_out = []
    if not issues and not anomalies:
        lines_out.append(f"{filename} is ready for docking.")
        title = "OK"
    else:
        title = "Warning"
        lines_out.append(f"PDB analysis for {filename}:")
        if issues:
            lines_out.append("Issues:")
            for item in issues:
                lines_out.append(f"  - {item}")
        if anomalies:
            lines_out.append("Anomalies:")
            for item in anomalies:
                lines_out.append(f"  - {item}")

    msg = "\n".join(lines_out)
    return msg, title

def _is_float(s):
    try:
        float(s)
        return True
    except:
        return False

class prepareTargetS1(QThread):
    requestResidueDialog = pyqtSignal()
    residuesReady = pyqtSignal(str)
    toastSignal = pyqtSignal(list)
    finishedSignal = pyqtSignal(list)

    def setWD(self, wd):
        self.wd = wd

    def checkVals(self, verifylist):
        self.verifylist = verifylist

    def dockMethod(self, text):
        self.dockingMethod = text

    def findProtein(self):
        targets = {}
        for potentialTargetProtein in glob('*'):
            base = os.path.splitext(potentialTargetProtein)[0]
            ext = os.path.splitext(potentialTargetProtein)[1].lower()
            proteinSize = os.path.getsize(potentialTargetProtein) / 1024
            if ext in ['.pdb', '.gpf', '.pdbqt']:
                if base not in targets:
                    targets[base] = {'pdb': None, 'gpf': None, 'pdbqt': None}
                if ext == '.pdb':
                    if proteinSize > 20:
                        with open(potentialTargetProtein, "r") as f:
                            fileContent = f.read()
                            if any(lig in fileContent for lig in LIGAND_TAGS) and any(res in fileContent for res in STD_RESIDUES):
                                self.toastSignal.emit([2000, "Found Co-crystalized ligand in the protein, Continuing..."])
                                targets[base]['pdb'] = potentialTargetProtein
                            elif not any(lig in fileContent for lig in LIGAND_TAGS) and any(res in fileContent for res in STD_RESIDUES):
                                targets[base]['pdb'] = potentialTargetProtein
                elif ext == '.gpf':
                    targets[base]['gpf'] = potentialTargetProtein
                elif ext == '.pdbqt':
                    if proteinSize > 20:
                        with open(potentialTargetProtein, "r") as f:
                            fileContent = f.read()
                            if any(lig in fileContent for lig in LIGAND_TAGS) and any(res in fileContent for res in STD_RESIDUES):
                                self.toastSignal.emit([2000, "Found Co-crystalized ligand in the protein, Continuing..."])
                                targets[base]['pdbqt'] = potentialTargetProtein
                            elif not any(lig in fileContent for lig in LIGAND_TAGS) and any(res in fileContent for res in STD_RESIDUES):
                                targets[base]['pdbqt'] = potentialTargetProtein

        # Filter: must have either pdb or pdbqt
        filtered = {k: v for k, v in targets.items() if v['pdb'] or v['pdbqt']}
        
        # Choose the one with the most populated fields
        def count_populated(entry):
            return sum(1 for v in entry.values() if v is not None)
        
        if filtered:
            best_key = max(filtered, key=lambda k: count_populated(filtered[k]))
            return {best_key: filtered[best_key]}
        else:
            return {None}
    
    def run(self):
        os.chdir(self.wd)
        self.toastSignal.emit([2000, "Assuming that appropriate option has been selected for bounding box calculations."])

        targets = self.findProtein()
        proteinBaseName, completeProtein = next(iter(targets.items()))

        if completeProtein['gpf'] and not completeProtein['pdb']:
            self.finishedSignal.emit([ f"❌ The target protein '{proteinBaseName}' was recognized as a GPF file but is missing a corresponding PDB file. Please verify your files and try again."])

        if all(completeProtein[attr] is not None for attr in ['pdb', 'gpf', 'pdbqt']):
            if self.dockingMethod == 'GPU':               
                self.finishedSignal.emit(['✅'])
            else:
                self.finishedSignal.emit([f"❌ The target protein '{proteinBaseName}' has corresponding PDB, GPF and PDBQT file but the docking mode is set to using CPU. Please verify your selection and try again."])
        elif (completeProtein['pdb'] and completeProtein['gpf'] and not completeProtein['pdbqt']) and self.dockingMethod == 'GPU':
            self.finishedSignal.emit(['✅'])
        elif (not completeProtein['pdb'] and not completeProtein['gpf'] and completeProtein['pdbqt']) and self.dockingMethod == 'CPU':
            if len(glob("conf.txt")) != 0:
                self.toastSignal.emit([2000, "Calculating bounding box dimensions\nusing user prepared configuration file..."])
                self.finishedSignal.emit(["Site", completeProtein, "CONF"])
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
                    self.finishedSignal.emit(["Site", completeProtein, self._residue_text])

                else:
                    self.finishedSignal.emit([f"❌ The target protein '{proteinBaseName}' was recognized as a PDBQT file but is missing a corresponding PDB file. Please verify your files and try again."])
            
            else:
                self.toastSignal.emit([2000, "Using default blind box dimensions..."])
                self.finishedSignal.emit(["V1", completeProtein, "BLIND"])

        elif (completeProtein['pdb'] and not completeProtein['pdbqt'] and not completeProtein['gpf']) or (completeProtein['pdb'] and completeProtein['pdbqt'] and not completeProtein['gpf']):
            if len(glob("conf.txt")) != 0:
                self.toastSignal.emit([2000, "Calculating bounding box dimensions\nusing user prepared configuration file..."])
                self.finishedSignal.emit(["Site", completeProtein, "CONF"])

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
                    self.finishedSignal.emit(["Site", completeProtein, self._residue_text])
                else:
                    self.toastSignal.emit([2000, "Calculating bounding box dimensions using IDv2's binding site predictor..."])
                    self.finishedSignal.emit(["Predict", completeProtein, "PREDICT"])
            
            else:
                self.toastSignal.emit([2000, "Using default blind box dimensions..."])
                self.finishedSignal.emit(["V1", completeProtein, "BLIND"]    )
            
        else:
            self.finishedSignal.emit(["❌ Please verify your files and try again."])

class prepareTargetS2(QThread):
    finishedSignal = pyqtSignal(list)

    def totalFile(self, file):
        self.mode = file[0]
        self.protein = file[1]
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
        lx = float(int(max(x) - min(x)) + 15)
        cy = float(int((sum(y)/len(y))*1000)) / 1000.0
        ly = float(int(max(y) - min(y)) + 15)
        cz = float(int((sum(z)/len(z))*1000)) / 1000.0
        lz = float(int(max(z) - min(z)) + 15)

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
            if self.tag == "CONF":
                conf = glob("conf.txt")[0]
                config = self.extract_coordinates(conf)
                if len(config) > 1:
                    center = [config[0], config[1], config[2]]
                    size = [config[3], config[4], config[5]]
                    self.finishedSignal.emit([center, size, self.protein])
                else:
                    self.finishedSignal.emit(config)

            else:
                specificResidues = self.tag.replace(' ', '').split(',')
                resDict = ID_utils.separateSingleResults(specificResidues)
                #Get results as {'A': ['256', '299', '302'], 'B': ['295', '315']} 
                maxRes, minRes = ID_utils.getMaxMinSingle(resDict)
                # Gives (256, 'A'), (310, 'B')
                if self.protein['pdbqt'] and not self.protein['pdb']:
                    first_coords = ID_utils.get_residue_coords(self.protein['pdbqt'], minRes[1], minRes[0])
                    last_coords = ID_utils.get_residue_coords(self.protein['pdbqt'], maxRes[1], maxRes[0])
                else:
                    first_coords = ID_utils.get_residue_coords(self.protein['pdb'], minRes[1], minRes[0])
                    last_coords = ID_utils.get_residue_coords(self.protein['pdb'], maxRes[1], maxRes[0])

                center, size = ID_utils.calculate_bounding_box(first_coords, last_coords)
                self.finishedSignal.emit([center, size, self.protein])

        elif self.mode == 'Predict':
            pdb_path = self.protein.get('pdb')
            try:
                kv_res = pyKVFinder.run_workflow(pdb_path)
            except Exception:
                self.finishedSignal.emit(["❌ No cavities found by IDv2!"])
                return

            vols = kv_res.volume
            cid  = max(vols, key=vols.get)
            unfiltered = kv_res.residues.get(cid, [])
            if not unfiltered:
                self.finishedSignal.emit(["❌ No residues for largest cavity."])
                return
            
            res_list = [entry for entry in unfiltered if entry[2] in STD_RESIDUES2]
            f_res, l_res = res_list[0], res_list[-1]
            c1 = ID_utils.get_residue_coords(self.protein['pdb'], f_res[1], int(f_res[0]))
            c2 = ID_utils.get_residue_coords(self.protein['pdb'], l_res[1], int(l_res[0]))

            center, size = ID_utils.calculate_bounding_box(c1, c2)
            self.finishedSignal.emit([center, size, self.protein])
            return
    
        else:
            # Fallback using default (blind) box dimensions.
            if self.protein['pdbqt'] and not self.protein['pdb']:
                self.proteinLocation = self.protein['pdbqt'] #PDBQT File
                center, size = self.IDv1Calc(self.proteinLocation)
            else:
                self.proteinLocation = self.protein['pdb'] #PDB File
                center, size = self.IDv1Calc(self.proteinLocation)
            
            self.finishedSignal.emit([center, size, self.protein])

class prepareTargetS3(QThread):
    toastSignal = pyqtSignal(list)
    finishedSignal = pyqtSignal(str)

    def allThings(self, text):
        self.center = text[0]
        self.size = text[1]
        self.proteinFile = text[2]

    def fixPDB(self, checkFix):
        self.checkFix = checkFix

    def dockMethod(self, method):
        self.method = method

    def fixPdb(self, proteinFile):
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
            
            # Back up the original file into a backup folder.
            backup_dir = os.path.join(os.path.dirname(proteinFile), "Backup")
            if not os.path.exists(backup_dir):
                os.makedirs(backup_dir)
            # Here we can rename the original file by appending a suffix.
            backup_file = os.path.join(backup_dir, f"{base_name}_original.pdb")
            shutil.move(proteinFile, backup_file)
            
            # Rename the temporary fixed file to the original filename.
            os.rename(temp_fixed_file, proteinFile)
            
            return "Done"
        except Exception:
            return "❌ Unable to fix the PDB, please review your files and try again."
        
    def convertPDB(self, name, proteinFile):        
        if self.method == "GPU":
            gpfFile = name + ".gpf"
            pdbqtFile = name + ".pdbqt"
            ar = [
                "--read_pdb", os.path.abspath(proteinFile),
                "-g", gpfFile, "--write_pdbqt", pdbqtFile,
                "--box_center", str(self.center[0]), str(self.center[1]), str(self.center[2]),
                "--box_size", str(int(self.size[0])), str(int(self.size[1])), str(int(self.size[2]))
            ]
        else:
            pdbqtFile = name + ".pdbqt"
            ar = [
                "--read_pdb", os.path.abspath(proteinFile),
                "-p", pdbqtFile,
                "--box_center", str(self.center[0]), str(self.center[1]), str(self.center[2]),
                "--box_size", str(int(self.size[0])), str(int(self.size[1])), str(int(self.size[2]))
            ]

        # Validate executable path
        primary_exe = resource_path("mk_prepare_receptor.exe")
        if not os.path.exists(primary_exe):
            return f"❌ Executable not found: {primary_exe}"

        # Run mk_prepare_receptor
        self.process.start(primary_exe, ar)
        self.process.waitForFinished(-1)
        output = self.process.readAllStandardOutput().data().decode()

        if "error" in output.lower() or "failed" in output.lower():
            if self.method == "GPU":
                return "❌ Could not prepare the receptor molecule for the selected docking method"
            # Fallback to OpenBabel
            inputform = '-ipdb'
            alt_exe = "obabel"
            ar = [inputform, proteinFile, "-opdbqt", "-O", pdbqtFile, "-xr", "-xc", "-xn", "--partialcharge", "gasteiger"]
            
            self.process.start(alt_exe, ar)
            self.process.waitForFinished(-1)
            output2 = self.process.readAllStandardOutput().data().decode()
            
            if "error" in output2.lower():
                return f"❌ Failed with OpenBabel: {output2}"
            else:
                return "Done"
        else:
            return "Done"
    
    def run(self):
        with open(self.proteinFile['pdb'], 'r') as f:
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
        with open(self.proteinFile['pdb'], 'w') as f:
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
                        if code in (" ", "") or code == keep:
                            dst.write(line[:16] + " " + line[17:])
                    else:
                        dst.write(line)

        pdb_path = self.proteinFile['pdb']
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
        minimal_args = ["--read_pdb", os.path.abspath(self.proteinFile['pdb'])]
        self.process.start(primary_exe, minimal_args)
        self.process.waitForFinished(-1)
        output = self.process.readAllStandardOutput().data().decode()

        if (("error" in output.lower() or "failed" in output.lower()) and (self.method == "GPU" or self.checkFix == 1)):
            self.toastSignal.emit([2000, "Fixing and conversion in progress..."])
            testTxt = self.fixPdb(self.proteinFile['pdb'])
        else:
            self.toastSignal.emit([2000, "Conversion in progress..."])
            testTxt = "DoneWithoutFix"
        
        if testTxt == "Done":
            name = self.proteinFile['pdb'].split('.')[0]
            proteinFixed = self.proteinFile['pdb']
            testConv = self.convertPDB(name, proteinFixed)    
            self.finishedSignal.emit(testConv)

        elif testTxt == "DoneWithoutFix":
            name = self.proteinFile['pdb'].split('.')[0]
            testConv = self.convertPDB(name, self.proteinFile['pdb'])               
            self.finishedSignal.emit(testConv)

        else:               
            self.finishedSignal.emit("❌ Failed the target conversion process")

class prepareLigand(QThread):
    toastSignal = pyqtSignal(list)
    finishedSignal = pyqtSignal(str)
    #Meeko always to preserve the outputs 
    def setWD(self, wd):
        self.wd = wd

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
    
    def moveLibs(self, libraries):
        cwd = os.getcwd()
        if len(libraries) > 1:
            libmulti = "Ligand_libraries_found"
            lib = libmulti
            if not os.path.exists(libmulti):
                os.makedirs(libmulti)
                for library in libraries:
                    shutil.move(os.path.join(cwd, library), libmulti)
        elif len(libraries) == 1:
            libsingle = "Ligand_library_found"
            lib = libsingle
            if not os.path.exists(libsingle):
                os.makedirs(libsingle)
                shutil.move(os.path.join(cwd, libraries), libsingle)            

        self.toastSignal.emit([2000, f"Libraries detected - transferring to:\n {os.path.abspath(lib)}"])
                
    def convertLig(self, ligFile):
        if ligFile.endswith(".pdbqt"):
            pass
        else:
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
            self.moveLibs(libraryFiles)
        ligandFiles, libraryFiles = self.findLigs()
        if libraryFiles:
            ID_utils.displayError(f"Library File(s) found while preparing ligands, please remove/relocate them and try again.")
        else:
            if ligandFiles:
                for ligandFile in ligandFiles:
                    self.convertLig(ligandFile)
                if len(ligandFiles) > 1:                    
                    self.finishedSignal.emit(f"✅ Ligands ready - {len(ligandFiles)} candidates for docking analysis!")
                else:                    
                    self.finishedSignal.emit(f"✅ Ligand {ligandFiles[0]} is ready for docking analysis.")
            else:               
                ID_utils.displayError("❌ No ligands found! Please check your files and try again.")
        
class prepareConfigS1(QThread):
    requestResidueDialog = pyqtSignal()
    residuesReady = pyqtSignal(str)
    toastSignal = pyqtSignal(list)
    finishedSignal = pyqtSignal(list)
    
    def setWD(self, wd):
        self.wd = wd
    
    def dockMethod(self, text):
        self.method = text

    def checkVals(self, verifylist):
        self.verifylist = verifylist

    def findProtein(self):
        targets = {}
        for potentialTargetProtein in glob('*'):
            base = os.path.splitext(potentialTargetProtein)[0]
            ext = os.path.splitext(potentialTargetProtein)[1].lower()
            proteinSize = os.path.getsize(potentialTargetProtein) / 1024
            if ext in ['.pdb', '.pdbqt']:
                if base not in targets:
                    targets[base] = {'pdb': None, 'pdbqt': None}
                if ext == '.pdb':
                    if proteinSize > 20:
                        with open(potentialTargetProtein, "r") as f:
                            fileContent = f.read()
                            if any(lig in fileContent for lig in LIGAND_TAGS) and any(res in fileContent for res in STD_RESIDUES):
                                self.toastSignal.emit([2000, "Found Co-crystalized ligand in the protein, Continuing..."])
                                targets[base]['pdb'] = potentialTargetProtein
                            elif not any(lig in fileContent for lig in LIGAND_TAGS) and any(res in fileContent for res in STD_RESIDUES):
                                targets[base]['pdb'] = potentialTargetProtein
                elif ext == '.pdbqt':
                    if proteinSize > 20:
                        with open(potentialTargetProtein, "r") as f:
                            fileContent = f.read()
                            if any(lig in fileContent for lig in LIGAND_TAGS) and any(res in fileContent for res in STD_RESIDUES):
                                self.toastSignal.emit([2000, "Found Co-crystalized ligand in the protein, Continuing..."])
                                targets[base]['pdbqt'] = potentialTargetProtein
                            elif not any(lig in fileContent for lig in LIGAND_TAGS) and any(res in fileContent for res in STD_RESIDUES):
                                targets[base]['pdbqt'] = potentialTargetProtein

        # Filter: must have either pdb or pdbqt
        filtered = {k: v for k, v in targets.items() if v['pdb'] or v['pdbqt']}

        def count_populated(entry):
            return sum(1 for v in entry.values() if v is not None)
        
        if filtered:
            best_key = max(filtered, key=lambda k: count_populated(filtered[k]))
            return {best_key: filtered[best_key]}
        else:
            return {None}            

    def run(self):
        os.chdir(self.wd)
        if len(glob("conf.txt")) != 0:
            self.finishedSignal.emit(["⚠️ IDv2 detected a configuration file. If you want to create a new configuration file, please remove the existing one and try again."])
            return
            
        if self.method == "AutoDock-GPU":         
            self.finishedSignal.emit(
                ["⚠️ AutoDock-GPU was selected as the docking method; therefore, a configuration file will not be generated."])
            return
        
        else:
            self.toastSignal.emit([2000, "Assuming that the target has already been prepared..."])
            targets = self.findProtein()
            proteinBaseName, completeProtein = next(iter(targets.items()))

            pdbTargetProtein = [completeProtein['pdb']]
            pdbqtTargetProtein = [completeProtein['pdbqt']]

            if pdbqtTargetProtein and not pdbTargetProtein:
                self.toastSignal.emit("Detected a compatible target PDBQT file. Continuing...")
                self.toastSignal.emit([2000, "Assuming that appropriate option has been selected for bounding box calculations."])
                
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
                        self.finishedSignal.emit(["Site", completeProtein, self._residue_text])

                    else:
                        self.finishedSignal.emit([f"❌ The target protein '{proteinBaseName}' was recognized as a PDBQT file but is missing a corresponding PDB file. Please verify your files and try again."])
                
                else:
                    self.toastSignal.emit([2000, "Using default blind box dimensions..."])
                    self.finishedSignal.emit(["V1", completeProtein, "BLIND"])

            elif pdbTargetProtein and pdbqtTargetProtein:
                self.toastSignal.emit([2000, "Assuming that appropriate option has been selected for bounding box calculations."])
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
                    self.finishedSignal.emit(["Site", completeProtein, self._residue_text])
                elif self.verifylist[1] == 1:
                    self.toastSignal.emit([2000, "Calculating bounding box dimensions using IDv2's binding site predictor..."])
                    self.finishedSignal.emit(["Predict", completeProtein, "PREDICT"])
                else:
                    self.toastSignal.emit([2000, "Using default blind box dimensions..."])
                    self.finishedSignal.emit(["V1", completeProtein, "BLIND"])

            elif pdbTargetProtein and not pdbqtTargetProtein:     
                self.finishedSignal.emit(
                    ["❌ Error: Only the target PDB file was found; the required PDBQT file is missing. Consequently, the configuration file could not be generated. Please review your input files and try again."])
            else:
                self.finishedSignal.emit(["❌ No suitable target protein found! Please check your files and try again."])

class prepareConfigS2(QThread):
    finishedSignal = pyqtSignal(list)

    def totalFile(self, file):
        self.mode = file[0]
        self.protein = file[1]
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
        lx = float(int(max(x) - min(x)) + 15)
        cy = float(int((sum(y)/len(y))*1000)) / 1000.0
        ly = float(int(max(y) - min(y)) + 15)
        cz = float(int((sum(z)/len(z))*1000)) / 1000.0
        lz = float(int(max(z) - min(z)) + 15)

        return [cx, cy, cz], [lx, ly, lz]

    def run(self):
        if self.mode == 'Site':
            specificResidues = self.tag.replace(' ', '').split(',')
            resDict = ID_utils.separateSingleResults(specificResidues)
            #Get results as {'A': ['256', '299', '302'], 'B': ['295', '315']} 
            maxRes, minRes = ID_utils.getMaxMinSingle(resDict)
            # Gives (256, 'A'), (310, 'B')
            if self.protein['pdbqt'] and not self.protein['pdb']:
                first_coords = ID_utils.get_residue_coords(self.protein['pdbqt'], minRes[1], minRes[0])
                last_coords = ID_utils.get_residue_coords(self.protein['pdbqt'], maxRes[1], maxRes[0])
            else:
                first_coords = ID_utils.get_residue_coords(self.protein['pdb'], minRes[1], minRes[0])
                last_coords = ID_utils.get_residue_coords(self.protein['pdb'], maxRes[1], maxRes[0])

            center, size = ID_utils.calculate_bounding_box(first_coords, last_coords)
            self.finishedSignal.emit([center, size, self.protein])

        elif self.mode == 'Predict':
            pdb_path = self.protein['pdb']
            try:
                kv_res = pyKVFinder.run_workflow(pdb_path)
            except Exception:
                self.finishedSignal.emit(["❌ No cavities found by IDv2!"])
                return

            vols = kv_res.volume
            cid  = max(vols, key=vols.get)
            unfiltered = kv_res.residues.get(cid, [])
            if not unfiltered:
                self.finishedSignal.emit(["❌ No residues for largest cavity."])
                return
            
            res_list = [entry for entry in unfiltered if entry[2] in STD_RESIDUES2]
            f_res, l_res = res_list[0], res_list[-1]
            c1 = ID_utils.get_residue_coords(self.protein['pdb'], f_res[1], int(f_res[0]))
            c2 = ID_utils.get_residue_coords(self.protein['pdb'], l_res[1], int(l_res[0]))

            center, size = ID_utils.calculate_bounding_box(c1, c2)
            self.finishedSignal.emit([center, size, self.protein])
            return
    
        else:
            self.proteinLocation = self.protein['pdbqt'] #PDBQT File
            center, size = self.IDv1Calc(self.proteinLocation)
            
            self.finishedSignal.emit([center, size, self.protein])

class prepareConfigS3(QThread):
    finishedSignal = pyqtSignal(str)
    
    def allThings(self, text):
        self.center = text[0]
        self.size = text[1]
        self.proteinFile = text[2]

    def run(self):
        with open("conf.txt", "w") as f:
            f.write("receptor = " + self.proteinFile['pdbqt'] + "\n\n")
            f.write(f"center_x = {self.center[0]:.2f}\n")
            f.write(f"center_y = {self.center[1]:.2f}\n")
            f.write(f"center_z = {self.center[2]:.2f}\n\n")
            f.write(f"size_x = {self.size[0]:.2f}\n")
            f.write(f"size_y = {self.size[1]:.2f}\n")
            f.write(f"size_z = {self.size[2]:.2f}\n")

        self.finishedSignal.emit("✅ Configuration file successfully prepared!")

class convertFiles(QThread):
    toastSignal = pyqtSignal(list)
    finishedSignal = pyqtSignal(str)

    def setWD(self, wd):
        self.wd = wd

    def run(self):
        # 1) Change into the working directory
        try:
            os.chdir(self.wd)
        except Exception:
            self.finishedSignal.emit(
                f"❌ Could not access working directory '{self.wd}'"
            )
            return

        # 2) Gather all files with these extensions (including any existing .pdbqt)
        patterns = ['*.mol2', '*.sdf', '*.pdb', '*.mol', '*.pdbqt']
        raw_files = []
        for pat in patterns:
            raw_files.extend(glob(pat))

        if not raw_files:
            self.finishedSignal.emit(
                "❌ No convertible files found. Please check the directory."
            )
            return

        # 3) Group by base filename (so we can skip any that already have a .pdbqt)
        groups = {}
        for fname in raw_files:
            base, ext = os.path.splitext(fname)
            ext = ext.lower().lstrip('.')
            groups.setdefault(base, set()).add(ext)

        obabel_exe = "obabel"

        # 4) Process each group
        for base, exts in groups.items():
            # Skip if we already have a .pdbqt for this base name
            if 'pdbqt' in exts:
                continue

            # Otherwise convert each source file in turn
            for ext in exts:
                src = f"{base}.{ext}"
                out = f"{base}.pdbqt"

                # Read the file once to decide protein vs. ligand
                try:
                    with open(src, 'r') as f:
                        content = f.read()
                except Exception:
                    continue

                is_protein = any(res in content for res in STD_RESIDUES)
                has_ligand = any(lig in content for lig in LIGAND_TAGS)

                # --- extra library‐file check for ligands ---
                if not is_protein and ext in ('sdf', 'mol2'):
                    if ext == 'sdf':
                        count = content.count('$$$$')
                        if count > 1:
                            continue
                    else:  # mol2
                        count = len(re.findall(r'@<TRIPOS>MOLECULE', content))
                        if count > 1:
                            continue

                # For proteins, always use the protein conversion args;
                # also notify if there was a co‑crystallized ligand.
                if is_protein:
                    if has_ligand:
                        self.toastSignal.emit([2000,
                            f"🔍 Found co-crystallized ligand in '{src}'."
                        ])
                    args = [
                        f"-i{ext}", src,
                        "-opdbqt", "-O", out,
                        "-xr", "-xn", "-xc", "-xp",
                        "--partialcharge", "gasteiger"
                    ]
                else:
                    # Ligand conversion (single‐molecule)
                    args = [
                        f"-i{ext}", src,
                        "-opdbqt", "-O", out,
                        "--partialcharge", "gasteiger",
                        "--DelNonPolarH"
                    ]

                # Run OpenBabel
                exit_code = QProcess.execute(obabel_exe, args)

        # 5) Finished
        self.finishedSignal.emit(
            "✅ Conversion complete: all missing .pdbqt files have been generated."
        )


