import re
import os
import csv
import math
import time
import shutil
import asyncio
import aiohttp
import pyKVFinder
import IDv2Resource
from glob import glob
from statistics import mean
from PyQt5 import QtWidgets
from urllib.parse import quote
from PyQt5.QtCore import QSize
from cryptography.fernet import Fernet
from PyQt5.QtGui import QIcon, QPixmap
from PyQt5.QtGui import QPixmap, QFont
from QSwitchControl import SwitchControl
from pyqttoast import Toast, ToastPreset
from PyQt5.QtWidgets import QInputDialog
from PyQt5.QtCore import QSize, Qt, QThread, pyqtSignal, QEasingCurve
from PyQt5.QtWidgets import (
    QPushButton, QGridLayout, QHBoxLayout, QLabel, QMainWindow, 
    QVBoxLayout, QWidget, QSpacerItem, QSizePolicy, QGroupBox,
    QDesktopWidget, QToolButton, QFrame, QFileDialog, QSpinBox
)
    

toastfont = QFont('Arial', 12, QFont.Weight.Medium)
key = b'XaQ5nlmxWQyPf9Q9bMjki5M7yEiXgOmTDRpvtTRWbik='

'''
SAMPLE single binding site prediction list
[['294', 'A', 'LEU'], ['295', 'A', 'GLY'], ['299', 'A', 'TYR'], ['302', 'A', 'VAL'], ['315', 'A', 'ALA'], ['345', 'A', 'VAL'], ['361', 'A', 'THR'], ['362', 'A', 'GLU'], ['363', 'A', 'TYR'], ['364', 'A', 'MET'], ['367', 'A', 'GLY'], ['407', 'A', 'HIS'], ['409', 'A', 'ASP'], ['413', 'A', 'ARG'], ['414', 'A', 'ASN'], ['416', 'A', 'LEU'], ['426', 'A', 'ALA'], ['427', 'A', 'ASP'], ['428', 'A', 'PHE'], ['429', 'A', 'GLY']] 
'''
LIGAND_TAGS = [' UNL ', ' UNK ', ' LIG ', ' <0> ', ' DRG ', ' INH ', ' NAG ', ' SO4 ', ' ATP ', ' ADP ', ' AMP ', ' HEM ', ' FMN ', ' FAD ', ' NAD ', ' GDP ', ' GTP ', ' SAM ']

STD_RESIDUES = [' ALA ', ' ARG ', ' ASN ', ' ASP ', ' CYS ', ' GLN ', ' GLU ', ' GLY ', ' HIS ', ' ILE ', ' LEU ', ' LYS ', ' MET ', ' PHE ', ' PRO ', ' SER ', ' THR ', ' TRP ', ' TYR ', ' VAL ']


def show_toast(text):
    toastfont = QFont('Arial', 12, QFont.Weight.Medium)
    toast = Toast()
    toast.setDuration(text[0])
    toast.setText(text[1])
    toast.setTextFont(toastfont)
    toast.applyPreset(ToastPreset.INFORMATION)
    toast.setResetDurationOnHover(False)
    toast.setBorderRadius(3) 
    toast.setIcon(QPixmap(':/IDLogo.png'))
    toast.setShowCloseButton(False) 
    toast.setIconSize(QSize(40, 40)) 
    toast.setMaximumOnScreen(1)
    toast.setIconColor(None)
    toast.show()

def separateSingleResults(textList):
    """
    Processes a list of strings in the format "<Residue Number><Chain Letter>"
    or "<Chain Letter><Residue Number>" and returns a dictionary keyed by the
    chain letter with a list of residue numbers.
    """
    bindingDict = {}
    for text in textList:
        m = re.match(r'^([A-Za-z])?(\d+)([A-Za-z])?$', text.strip())
        if not m:
            continue
        letter1, number, letter2 = m.groups()
        chain = letter1 or letter2
        if not chain:
            continue
        bindingDict.setdefault(chain, []).append(number)

    return bindingDict

def separateMultiResults(textList):
    """
    Processes a list of strings in the format:
      "<Residue Number><Chain Letter> [Protein Name]" or "<Chain Letter><Residue Number> [Protein Name]"
    Returns a dict keyed by protein name with values as a list of (chain letter, residue number) tuples.
    """
    bindingDict = {}
    for text in textList:
        # Get protein name from square brackets
        protein = re.search(r'\[([^\]]+)\]', text)
        if not protein:
            continue
        protein = protein.group(1).strip()
        
        # Get residue info from text before '[' and extract optional letters and number.
        prefix = text.split('[')[0].strip()
        m = re.match(r'^([A-Za-z])?(\d+)([A-Za-z])?$', prefix)
        if not m:
            continue
        letter1, number, letter2 = m.groups()
        chain = letter1 or letter2
        if not chain:
            continue
        
        bindingDict.setdefault(protein, []).append((chain, number))

    return bindingDict

def get_max_min(dictionary): #Obsolete
    # Find the max and min values and their corresponding keys
    max_key = max(dictionary, key=dictionary.get)
    min_key = min(dictionary, key=dictionary.get)
    return max_key, dictionary[max_key], min_key, dictionary[min_key]

def getMaxMinSingle(bindingDict):
    all_residues = []

    for chain, residues in bindingDict.items():
        for res in residues:
            all_residues.append((int(res), chain))  # Store (residue, chain) as a tuple

    if all_residues:
        min_res, min_chain = min(all_residues)  # Get the tuple with the lowest residue
        max_res, max_chain = max(all_residues)  # Get the tuple with the highest residue
        return (min_res, min_chain), (max_res, max_chain)

    return None, None

def getMaxMinMulti(bindingDict):
    result = {}
    for protein, entries in bindingDict.items():
        # Find the entry with the minimum residue number.
        min_entry = min(entries, key=lambda x: int(x[1]))
        max_entry = max(entries, key=lambda x: int(x[1]))
        # Format: (min_residue, chain) and (max_residue, chain)
        result[protein] = ((int(min_entry[1]), min_entry[0]), (int(max_entry[1]), max_entry[0]))
    return result

def get_residue_coords(pdb_file, chain_id, residue_number):
    coords = None
    residue_found = False
    
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("ATOM"):  # Only consider ATOM records
                atom_name = line[12:16].strip()
                chain = line[21].strip()
                res_num = int(line[22:26].strip())
                
                if chain == chain_id and res_num == residue_number:
                    residue_found = True
                    if atom_name == "CA" or (atom_name == "C" and coords is None):
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        coords = (x, y, z)
            elif residue_found and line.startswith("TER"):  # Exit early when residue is complete
                break
    
    return coords

def calculate_bounding_box(first_coords, last_coords, buffer=10):
    x_min, y_min, z_min = min(first_coords[0], last_coords[0]), min(first_coords[1], last_coords[1]), min(first_coords[2], last_coords[2])
    x_max, y_max, z_max = max(first_coords[0], last_coords[0]), max(first_coords[1], last_coords[1]), max(first_coords[2], last_coords[2])

    # Center of the bounding box
    center = ((x_min + x_max) / 2, (y_min + y_max) / 2, (z_min + z_max) / 2)

    # Size with buffer if specified
    size = [
        math.ceil(x_max - x_min + 2 * buffer),
        math.ceil(y_max - y_min + 2 * buffer),
        math.ceil(z_max - z_min + 2 * buffer)
    ]

    return center, size

def generateConfig(wd):
    proteinLocation = wd
    receptorName = os.path.basename(proteinLocation)
    xval, yval, zval = [], [], []
    x, y, z = [], [], []
    with open(proteinLocation, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                xtest = line[30:39]
                xval.append(xtest)
                ytest = line[38:46]
                yval.append(ytest)
                ztest = line[46:54]
                zval.append(ztest)
    for i in range(len(xval)):
        x.append(float(xval[i]))
        y.append(float(yval[i]))
        z.append(float(zval[i]))
    
    cx = str(int((mean(x))*1000)/1000.)
    lx = str(int(max(x) - min(x)) + 15)
    cy = str(int((mean(y))*1000)/1000.)
    ly = str(int(max(y) - min(y)) + 15)
    cz = str(int((mean(z))*1000)/1000.)
    lz = str(int(max(z) - min(z)) + 15)

    text = "receptor = "+receptorName+"\n\ncenter_x = "+cx+"\ncenter_y = "+cy+"\ncenter_z = "+cz+"\n\nsize_x = "+lx+"\nsize_y = "+ly+"\nsize_z = "+lz
    with open("conf.txt", 'w') as f:
        f.write(text)
        f.close()

class RateLimiter:
    """
    A simple asynchronous rate limiter that allows up to max_calls per period seconds.
    """
    def __init__(self, max_calls, period):
        self.max_calls = max_calls
        self.period = period
        self.calls = []
        self.lock = asyncio.Lock()

    async def acquire(self):
        async with self.lock:
            now = time.monotonic()
            # Remove calls that are older than the period.
            self.calls = [t for t in self.calls if t > now - self.period]
            if len(self.calls) >= self.max_calls:
                # Calculate time to sleep until the oldest call is outside the period.
                sleep_time = self.calls[0] + self.period - now
                await asyncio.sleep(sleep_time)
            self.calls.append(time.monotonic())

async def download_sdf_pubchem(session, identifier, rate_limiter, downloader, idx):
    """
    Downloads the 3D SDF file for the given compound identifier.
    If the identifier is all digits, it is assumed to be a CID; otherwise, a name.
    The file is saved as '<working_directory>/<identifier>.sdf'.
    Once finished (or failed), the downloader’s progress signal is emitted.
    """
    # Determine the input type and encode if needed.
    if identifier.isdigit():
        input_spec = f"cid/{identifier}"
    else:
        encoded_name = quote(identifier)
        input_spec = f"name/{encoded_name}"
    
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{input_spec}/SDF?record_type=3d"
    show_toast([1000, f"Downloading {identifier} from {url}"])
    

    # Enforce rate limiting.
    await rate_limiter.acquire()

    try:
        async with session.get(url) as response:
            if response.status == 200:
                data = await response.read()
                # Save file to the working directory (downloader.wd).
                output_filename = f"{downloader.wd}/{identifier}.sdf"
                with open(output_filename, "wb") as f:
                    f.write(data)
    except Exception:
        downloader.toastSignal.emit(f"❌ Download failed")
    # Emit progress update after finishing this identifier.
    downloader.progressBarSignal.emit(idx + 1)

class pubchemEngine(QThread):
    toastSignal = pyqtSignal(str)
    finishedSignal = pyqtSignal(str)
    progressBarSignal = pyqtSignal(int)
    progressBarSignalTotal = pyqtSignal(int)

    def setWD(self, wd):
        self.wd = wd
    
    def setIDs(self, ids):
        self.ids = ids

    def setFormat(self, formats):
        self.formats = formats

    def run(self):
        """
        This run() method creates an asyncio event loop and runs asynchronous downloads.
        Depending on the format flag (e.g. self.formats[0]==1), you can handle compound IDs or names.
        """
        # Split the provided IDs/names into a list.
        identifiers = re.split(r'[\n]+', self.ids)
        # Remove any empty strings.
        identifiers = [id.strip() for id in identifiers if id.strip()]
        if not identifiers:
            self.finishedSignal.emit("❌ Please enter appropriate compound IDs or names from PubChem and try again.")
            return

        # Emit total progress for UI (total number of downloads).
        self.progressBarSignalTotal.emit(len(identifiers))

        # Run the asynchronous downloading tasks in a new event loop.
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        try:
            loop.run_until_complete(self.download_all_compounds(identifiers))
            self.finishedSignal.emit(f"✅ Download complete! All files have been successfully saved to {self.wd}.")
        except Exception:
            self.finishedSignal.emit(f"❌ Download failed")
        finally:
            loop.close()

    async def download_all_compounds(self, identifiers):
        """
        This asynchronous method creates the rate limiter and client session, then schedules downloads.
        """
        rate_limiter = RateLimiter(max_calls=5, period=1.0)
        async with aiohttp.ClientSession() as session:
            # Schedule download tasks.
            tasks = []
            for idx, identifier in enumerate(identifiers):
                tasks.append(download_sdf_pubchem(session, identifier, rate_limiter, self, idx))
            await asyncio.gather(*tasks)

async def download_sdf_zinc(session, zinc_id, rate_limiter, downloader, idx):
    """
    Downloads the SDF file from Zinc for the given ZINC id.
    The ZINC id must exactly start with "ZINC" (case sensitive).
    Saves the file as '<working_directory>/<zinc_id>.sdf'.
    Once finished (or failed), the downloader’s progress signal is emitted.
    """
    if not zinc_id.startswith("ZINC"):
        displayError(f"Error: Identifier '{zinc_id}' is not a valid ZINC id (must start with 'ZINC'). Skipping.")
        downloader.progressBarSignal.emit(idx + 1)
        return

    url = f"https://zinc.docking.org/substances/{zinc_id}.sdf"
    show_toast([1000, f"Downloading {zinc_id} from {url}"])

    # Enforce rate limiting.
    await rate_limiter.acquire()

    try:
        async with session.get(url) as response:
            if response.status == 200:
                data = await response.read()
                # Save file to the working directory (downloader.wd).
                output_filename = f"{downloader.wd}/{zinc_id}.sdf"
                with open(output_filename, "wb") as f:
                    f.write(data)
    except Exception:
        downloader.toastSignal.emit(f"❌ Download failed")
    # Emit progress update after finishing this identifier.
    downloader.progressBarSignal.emit(idx + 1)

class zincEngine(QThread):
    toastSignal = pyqtSignal(str)
    finishedSignal = pyqtSignal(str)
    progressBarSignal = pyqtSignal(int)
    progressBarSignalTotal = pyqtSignal(int)

    def setWD(self, wd):
        self.wd = wd
    
    def setIDs(self, ids):
        self.ids = ids

    def run(self):
        """
        This run() method creates an asyncio event loop and runs asynchronous downloads for Zinc.
        """
        # Split the provided Zinc IDs into a list.
        zinc_ids = re.split(r'[,\t\n]+', self.ids)
        # Remove any empty strings.
        zinc_ids = [z.strip() for z in zinc_ids if z.strip()]
        if not zinc_ids:
            self.finishedSignal.emit("❌ Please enter appropriate ZINC IDs and try again.")
            return

        # Emit total progress for UI (total number of downloads).
        self.progressBarSignalTotal.emit(len(zinc_ids))

        # Run the asynchronous downloading tasks in a new event loop.
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        try:
            loop.run_until_complete(self.download_all_zinc(zinc_ids))
            self.finishedSignal.emit(f"✅ Download complete! All files have been successfully saved to {self.wd}.")
        except Exception:
            self.finishedSignal.emit(f"❌ Download failed")
        finally:
            loop.close()

    async def download_all_zinc(self, zinc_ids):
        """
        This asynchronous method creates the rate limiter and client session, then schedules Zinc downloads.
        """
        rate_limiter = RateLimiter(max_calls=5, period=1.0)
        async with aiohttp.ClientSession() as session:
            tasks = []
            for idx, zinc_id in enumerate(zinc_ids):
                tasks.append(download_sdf_zinc(session, zinc_id, rate_limiter, self, idx))
            await asyncio.gather(*tasks)

async def download_imppat_file(session, filetype, imppat_id, rate_limiter, downloader, idx):
    """
    Downloads a file from the IMPPAT images server.
    
    Constructs the URL as:
      https://cb.imsc.res.in/imppat/images/3D/{filetype}/{imppat_id}_3D.{ext}
    where {ext} is the lowercase version of filetype.
    
    Saves the file in downloader.wd and emits a progress update via downloader.progressBarSignal.
    """
    
    ext = filetype.lower()
    url = f"https://cb.imsc.res.in/imppat/images/3D/{filetype}/{imppat_id}_3D.{ext}"
    show_toast([1000, f"Downloading {imppat_id} from {url}"])

    # Enforce rate limiting.
    await rate_limiter.acquire()

    try:
        async with session.get(url) as response:
            if response.status == 200:
                data = await response.read()
                output_filename = f"{downloader.wd}/{imppat_id}_3D.{ext}"
                with open(output_filename, "wb") as f:
                    f.write(data)
    except Exception:
        downloader.toastSignal.emit(f"❌ Download failed")
    # Emit progress update after finishing this download.
    downloader.progressBarSignal.emit(idx + 1)

class imppatEngine(QThread):
    toastSignal = pyqtSignal(str)
    finishedSignal = pyqtSignal(str)
    progressBarSignal = pyqtSignal(int)
    progressBarSignalTotal = pyqtSignal(int)

    def setWD(self, wd):
        self.wd = wd

    def setIDs(self, ids):
        self.ids = ids

    def setFormat(self, formats):
        self.formats = formats

    def run(self):
        """
        Creates an asyncio event loop and runs asynchronous downloads for IMPPAT.
        Expects:
        - self.ids as a string of IMPPAT identifiers.
        - self.formats as a bool list representing [SDF, MOL2, PDB] (e.g., [1, 0, 0] for SDF).
        """
        # Split the provided IMPPAT IDs into a list.
        imppat_ids = re.split(r'[,\t\n]+', self.ids)
        imppat_ids = [imppat_id.strip() for imppat_id in imppat_ids if imppat_id.strip()]
        if not imppat_ids:
            self.finishedSignal.emit("❌ Please enter appropriate IMPPAT IDs and try again.")
            return

        # Determine filetype from the bool list: [SDF, MOL2, PDB].
        filetype = None
        if isinstance(self.formats, list) and len(self.formats) >= 3:
            if self.formats[0]:
                filetype = "SDF"
            elif self.formats[1]:
                filetype = "MOL2"
            elif self.formats[2]:
                filetype = "PDB"
        
        # Emit the total number of downloads for the progress bar.
        self.progressBarSignalTotal.emit(len(imppat_ids))

        # Create a new asyncio event loop for this thread.
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        try:
            loop.run_until_complete(self.download_all_imppat(filetype, imppat_ids))
            self.finishedSignal.emit(f"✅ Download complete! All files have been successfully saved to {self.wd}.")
        except Exception:
            self.finishedSignal.emit(f"❌ Download failed")
        finally:
            loop.close()

    async def download_all_imppat(self, filetype, imppat_ids):
        """
        Creates the rate limiter and aiohttp session, then schedules downloads for all IMPPAT IDs.
        """
        rate_limiter = RateLimiter(max_calls=5, period=1.0)
        async with aiohttp.ClientSession() as session:
            tasks = []
            for idx, imppat_id in enumerate(imppat_ids):
                tasks.append(download_imppat_file(session, filetype, imppat_id, rate_limiter, self, idx))
            await asyncio.gather(*tasks)

class bindingPred(QThread):
    centerSignal = pyqtSignal(list)
    sizeSignal = pyqtSignal(list)

    def getProtein(self, location):
        self.proteinLocation = location

    def run(self):
        predResults = pyKVFinder.run_workflow(self.proteinLocation)

        volRes = predResults.volume
        cavity = max(volRes, key=volRes.get)

        specificResults = predResults.residues
        first_residue = specificResults.get(cavity)[0]
        last_residue = specificResults.get(cavity)[-1]

        first_coords = get_residue_coords(self.proteinLocation, first_residue[1], int(first_residue[0]))
        last_coords = get_residue_coords(self.proteinLocation, last_residue[1], int(last_residue[0]))

        center, size = calculate_bounding_box(first_coords, last_coords)
        
        self.centerSignal.emit(center)
        self.sizeSignal.emit(size)
    
class siteSpecificPred(QThread):
    def getProtein(self, location):
        self.proteinLocation = location

    def residues(self, specificResidues):
        self.specificResidues = specificResidues

    def run(self, specificResidues):
        specificResidues = specificResidues.replace(' ', '').split(',')
        resDict = separateSingleResults(specificResidues)
        max_key, max_value, min_key, min_value = get_max_min(resDict)
        first_coords = get_residue_coords(self.proteinLocation, min_key, min_value)
        last_coords = get_residue_coords(self.proteinLocation, max_key, max_value)

        center, size = calculate_bounding_box(first_coords, last_coords)
        return center, size
    
def displayError(text):
    """
    Displays a standard error message box with custom styling.
    
    :param text: The error message to display.
    """
    error_box = QtWidgets.QMessageBox()
    error_box.setIcon(QtWidgets.QMessageBox.Critical)
    error_box.setText(text)
    error_box.setWindowIcon(QIcon(':/IDLogo.png'))
    error_box.setWindowTitle("Error")
    error_box.setStandardButtons(QtWidgets.QMessageBox.Ok)
    
    # Apply a custom style sheet to match your UI theme.
    error_box.setStyleSheet("""
        QMessageBox {
            background-color: #F5F7FA;  /* light gradient background can be simulated with a solid color */
            font-family: Arial;
            font-size: 14px;
            color: #000000;
        }
        QPushButton {
            font-size: 14px;
            padding: 4px 10px;
            border: 1px solid #A9A9A9;
            border-radius: 3px;
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1,
                                        stop:0 #F2F2F2, stop:0.5 #EBEBEB, stop:1 #CFCFCF);
        }
        QPushButton:hover {
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1,
                                        stop:0 #EAF6FD, stop:0.5 #D9F0FC, stop:1 #A7D9F5);
        }
        QPushButton:pressed {
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1,
                                        stop:0 #E5F4FC, stop:0.5 #C4E5F6, stop:1 #68B3DB);
        }
    """)
    
    error_box.exec_()

def displayWarn(text):
    """
    Displays a standard error message box with custom styling.
    
    :param text: The error message to display.
    """
    error_box = QtWidgets.QMessageBox()
    error_box.setIcon(QtWidgets.QMessageBox.Warning)
    error_box.setText(text)
    error_box.setWindowIcon(QIcon(':/IDLogo.png'))
    error_box.setWindowTitle("Warning")
    error_box.setStandardButtons(QtWidgets.QMessageBox.Ok)
    
    # Apply a custom style sheet to match your UI theme.
    error_box.setStyleSheet("""
        QMessageBox {
            background-color: #F5F7FA;  /* light gradient background can be simulated with a solid color */
            font-family: Arial;
            font-size: 14px;
            color: #000000;
        }
        QPushButton {
            font-size: 14px;
            padding: 4px 10px;
            border: 1px solid #A9A9A9;
            border-radius: 3px;
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1,
                                        stop:0 #F2F2F2, stop:0.5 #EBEBEB, stop:1 #CFCFCF);
        }
        QPushButton:hover {
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1,
                                        stop:0 #EAF6FD, stop:0.5 #D9F0FC, stop:1 #A7D9F5);
        }
        QPushButton:pressed {
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1,
                                        stop:0 #E5F4FC, stop:0.5 #C4E5F6, stop:1 #68B3DB);
        }
    """)
    
    error_box.exec_()

def displayWarning(text):
    warning_box = QtWidgets.QMessageBox()
    warning_box.setIcon(QtWidgets.QMessageBox.Warning)
    warning_box.setWindowIcon(QIcon(':/IDLogo.png'))
    warning_box.setWindowTitle("Warning!")

    # A fixed summary shown as the main text.
    warning_box.setText("Click 'Show Details...' for full report.")
    warning_box.setDetailedText(text)
    warning_box.setStandardButtons(QtWidgets.QMessageBox.Ok)

    # Apply a custom style sheet to match your UI theme.
    warning_box.setStyleSheet("""
        QMessageBox {
            background-color: #F5F7FA;
            font-family: Arial;
            font-size: 14px;
            color: #000000;
        }
        QPushButton {
            font-size: 14px;
            padding: 4px 10px;
            border: 1px solid #A9A9A9;
            border-radius: 3px;
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1,
                                        stop:0 #F2F2F2, stop:0.5 #EBEBEB, stop:1 #CFCFCF);
        }
        QPushButton:hover {
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1,
                                        stop:0 #EAF6FD, stop:0.5 #D9F0FC, stop:1 #A7D9F5);
        }
        QPushButton:pressed {
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1,
                                        stop:0 #E5F4FC, stop:0.5 #C4E5F6, stop:1 #68B3DB);
        }
    """)

    warning_box.exec_()

def displayMess(text):
    sucess_box = QtWidgets.QMessageBox()
    sucess_box.setIcon(QtWidgets.QMessageBox.Information)
    sucess_box.setWindowIcon(QIcon(':/IDLogo.png'))
    sucess_box.setWindowTitle("Success!")

    # A fixed summary shown as the main text.
    sucess_box.setText("Click 'Show Details...' for full report.")
    sucess_box.setDetailedText(text)
    sucess_box.setStandardButtons(QtWidgets.QMessageBox.Ok)

    # Apply a custom style sheet to match your UI theme.
    sucess_box.setStyleSheet("""
        QMessageBox {
            background-color: #F5F7FA;
            font-family: Arial;
            font-size: 14px;
            color: #000000;
        }
        QPushButton {
            font-size: 14px;
            padding: 4px 10px;
            border: 1px solid #A9A9A9;
            border-radius: 3px;
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1,
                                        stop:0 #F2F2F2, stop:0.5 #EBEBEB, stop:1 #CFCFCF);
        }
        QPushButton:hover {
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1,
                                        stop:0 #EAF6FD, stop:0.5 #D9F0FC, stop:1 #A7D9F5);
        }
        QPushButton:pressed {
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1,
                                        stop:0 #E5F4FC, stop:0.5 #C4E5F6, stop:1 #68B3DB);
        }
    """)

    sucess_box.exec_()


def displayMessage(text):
    """
    Displays a standard error message box with custom styling.
    
    :param text: The error message to display.
    """
    sucess_box = QtWidgets.QMessageBox()
    sucess_box.setIcon(QtWidgets.QMessageBox.Information)
    sucess_box.setText(text)
    sucess_box.setWindowTitle("Sucess!")
    sucess_box.setWindowIcon(QIcon(':/IDLogo.png'))
    sucess_box.setStandardButtons(QtWidgets.QMessageBox.Ok)
    
    # Apply a custom style sheet to match your UI theme.
    sucess_box.setStyleSheet("""
        QMessageBox {
            background-color: #F5F7FA;  /* light gradient background can be simulated with a solid color */
            font-family: Arial;
            font-size: 14px;
            color: #000000;
        }
        QPushButton {
            font-size: 14px;
            padding: 4px 10px;
            border: 1px solid #A9A9A9;
            border-radius: 3px;
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1,
                                        stop:0 #F2F2F2, stop:0.5 #EBEBEB, stop:1 #CFCFCF);
        }
        QPushButton:hover {
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1,
                                        stop:0 #EAF6FD, stop:0.5 #D9F0FC, stop:1 #A7D9F5);
        }
        QPushButton:pressed {
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1,
                                        stop:0 #E5F4FC, stop:0.5 #C4E5F6, stop:1 #68B3DB);
        }
    """)
    
    sucess_box.exec_()

def displayComplete(text):
    """
    Displays a standard error message box with custom styling.
    
    :param text: The error message to display.
    """
    complete = QtWidgets.QMessageBox()
    complete.setText(text)
    complete.setWindowTitle("Process Complete!")
    complete.setWindowIcon(QIcon(':/IDLogo.png'))
    complete.setStandardButtons(QtWidgets.QMessageBox.Ok)
    
    # Apply a custom style sheet to match your UI theme.
    complete.setStyleSheet("""
        QMessageBox {
            background-color: #F5F7FA;  /* light gradient background can be simulated with a solid color */
            font-family: Arial;
            font-size: 14px;
            color: #000000;
        }
        QPushButton {
            font-size: 14px;
            padding: 4px 10px;
            border: 1px solid #A9A9A9;
            border-radius: 3px;
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1,
                                        stop:0 #F2F2F2, stop:0.5 #EBEBEB, stop:1 #CFCFCF);
        }
        QPushButton:hover {
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1,
                                        stop:0 #EAF6FD, stop:0.5 #D9F0FC, stop:1 #A7D9F5);
        }
        QPushButton:pressed {
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1,
                                        stop:0 #E5F4FC, stop:0.5 #C4E5F6, stop:1 #68B3DB);
        }
    """)
    
    complete.exec_()


class CustomTitleBar(QWidget):
    def __init__(self, parent):
        super().__init__(parent)
        self.initial_pos = None
        title_bar_layout = QHBoxLayout(self)
        title_bar_layout.setContentsMargins(1, 1, 1, 1)
        title_bar_layout.setSpacing(2)

        self.title = QLabel(f"{self.__class__.__name__}", self)
        self.title.setAlignment(Qt.AlignLeft)
        self.title.setStyleSheet(
            """
        QLabel { font-size: 10pt; margin-left: 5px; }
        """
        )
        if title := parent.windowTitle():
            self.title.setText(title)
        title_bar_layout.addWidget(self.title)

        # Spacer to push buttons to the right
        title_bar_layout.addSpacerItem(QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum))

        # Min button
        self.min_button = QToolButton(self)
        self.min_button.setText("\u005F")        
        self.min_button.setFixedSize(QSize(24, 24))
        self.min_button.clicked.connect(self.window().showMinimized)
        self.min_button.setFocusPolicy(Qt.NoFocus)
        #self.min_button.setStyleSheet('''
        #    QToolButton::menu-indicator { image: none; } 
        #    QToolButton { background-color: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:1, stop:0 rgba(206, 255, 188, 255), 
        #                              stop:0.255682 rgba(212, 255, 179, 255), stop:0.636364 rgba(168, 255, 89, 255), stop:1 rgba(122, 255, 0, 255));
        #                              padding: 2px; border-radius: 6px; border: 1px solid #b5b3b5;
        #        padding: 20px; font-size: 10pt;}
        #    QToolButton::hover { background-color: lightgreen;}
        #                              ''')
        title_bar_layout.addWidget(self.min_button)

        # Close button
        self.close_button = QToolButton(self)
        self.close_button.setText("\u00D7")
        self.close_button.setStyleSheet("color:darkred;")
        self.close_button.setFixedSize(QSize(24, 24))
        self.close_button.clicked.connect(self.window().close)
        self.close_button.setFocusPolicy(Qt.NoFocus)
        #self.close_button.setStyleSheet('''
        #    QToolButton::menu-indicator { image: none; } 
        #    QToolButton { background-color: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:1, stop:0 rgba(255, 188, 188, 255), 
        #                                stop:0.255682 rgba(255, 179, 179, 255), stop:0.636364 rgba(255, 89, 89, 255), stop:1 rgba(255, 0, 0, 255));
        #                                padding: 2px; border-radius: 6px; border: 1px solid #b5b3b5; color: black; font-weight:bold;
        #        padding: 20px; font-size: 10pt;}
        #    QToolButton::hover {background-color: red;}
        #                               ''')
        title_bar_layout.addWidget(self.close_button)

    def changeEvent(self, event):
        super().changeEvent(event)
        event.accept()

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            self.initial_pos = event.pos()
        super().mousePressEvent(event)
        event.accept()

    def mouseMoveEvent(self, event):
        try:
            if self.initial_pos is not None:
                delta = event.pos() - self.initial_pos
                self.window().move(
                    self.window().x() + delta.x(),
                    self.window().y() + delta.y(),
                )
        except Exception:
            pass
        super().mouseMoveEvent(event)
        event.accept()

    def mouseReleaseEvent(self, event):
        self.initial_pos = None
        super().mouseReleaseEvent(event)
        event.accept()

class librarySplitter(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.setupUi(self)
    def center_above_screen(self):
        screen_geometry = QDesktopWidget().screenGeometry()
        window_geometry = self.frameGeometry()
        screen_center = screen_geometry.center()
        screen_center.setX(screen_center.x() + 200)
        screen_center.setY(screen_center.y() - 30)
        window_geometry.moveCenter(screen_center)
        self.move(window_geometry.topLeft())

    def __init__(self):
        super().__init__()
        self.setWindowIcon(QIcon(':/IDLogo.png'))
        self.setWindowTitle("IDv2 Library Splitter")
        self.resize(300, 300)
        self.is_dark_mode = False
        self.center_above_screen()
        self.setStyleSheet("""
        QMainWindow {
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:1, 
            stop:0 #F5F7FA, 
            stop:0.5 #E4E7EB, 
            stop:1 #D9DEE4); /* Light and minimalistic gradient */
            color: #000000; /* Default dark text */
            border: 1px solid #A9A9A9;
        }
        
        QLabel{
            background: transparent;
            font-size: 14px;
            font-family: "Arial";
        }
        
        QGroupBox{
            border:0;
            background: transparent;
        }

        QPushButton {
            font-size: 14px; /* Matches the font size */
            padding: 2px 8px; /* Inner spacing */
            border: 1px solid #A9A9A9; /* Default border */
            border-radius: 3px; /* Rounded corners */
            color: black; /* Default text color */
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, stop:0 #F2F2F2, stop:0.5 #EBEBEB, stop:0.51 #DDDDDD, stop:1 #CFCFCF); /* Subtle 3D gradient */
        }

        QPushButton:hover {
            border: 1px solid #A9A9A9; /* Light border on hover */
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, stop:0 #EAF6FD, stop:0.5 #D9F0FC, stop:0.51 #BEE6FD, stop:1 #A7D9F5); /* Light blue gradient on hover */
        }

        QPushButton:pressed {
            padding: 2px 7px 3px 9px; /* Slight padding shift for pressed effect */
            border: 1px solid #A9A9A9; /* Darker border when pressed */
            border-bottom: 0px; /* Simulates "pressing in" by hiding the bottom border */
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, stop:0 #E5F4FC, stop:0.5 #C4E5F6, stop:0.51 #98D1EF, stop:1 #68B3DB); /* Darker gradient when pressed */
        }
                            
        QPushButton:checked {
            padding: 2px 7px 3px 9px; /* Slight padding shift for pressed effect */
            border: 1px solid #A9A9A9; /* Darker border when pressed */
            border-bottom: 0px; /* Simulates "pressing in" by hiding the bottom border */
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, stop:0 #E5F4FC, stop:0.5 #C4E5F6, stop:0.51 #98D1EF, stop:1 #68B3DB); /* Darker gradient when pressed */
        }

        QPushButton#startButton{
            font-size: 35px;
        }
                            
    """
)
        sizePolicy = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sizePolicy)
        self.setMinimumSize(QSize(300, 300))
        self.setWindowFlags(Qt.FramelessWindowHint)
        central_widget = QWidget()
        central_widget.setObjectName("Container")
        self.title_bar = CustomTitleBar(self)

        central_widget_layout = QVBoxLayout()
        central_widget_layout.setContentsMargins(0, 0, 0, 0)
        central_widget_layout.setAlignment(Qt.AlignTop)

        central_widget_layout.addWidget(self.title_bar)

        work_space_layout = QVBoxLayout()
        work_space_layout.setContentsMargins(11, 11, 11, 11)

        sizePolicy1 = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Minimum)
        sizePolicy1.setHorizontalStretch(0)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy2 = QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)
        sizePolicy2.setHorizontalStretch(0)
        sizePolicy2.setVerticalStretch(0)

        self.introLabel = QLabel(self)
        self.introLabel.setObjectName("introLabel")
        self.introLabel.setTextFormat(Qt.RichText)
        self.introLabel.setText('''
<span style="font-size: 14px; font-family: Arial;"> 
  This submodule is designed to efficiently divide a single library file, such as those formatted in MOL2 or SDF, into multiple distinct ligand files, thereby streamlining data organization and enhancing the ease of subsequent analyses.
</span>
                                ''')
        self.introLabel.setWordWrap(True)
        sizePolicy2.setHeightForWidth(self.introLabel.sizePolicy().hasHeightForWidth())
        self.introLabel.setSizePolicy(sizePolicy2)
        self.introLabel.setAlignment(Qt.AlignJustify)

        work_space_layout.addWidget(self.introLabel)

        self.groupWD = QGroupBox(self)
        self.groupWD.setObjectName("groupWD")
        
        sizePolicy2.setHeightForWidth(self.groupWD.sizePolicy().hasHeightForWidth())
        self.groupWD.setSizePolicy(sizePolicy2)
        self.groupWD.setMinimumSize(QSize(280, 75))

        self.gridLayout1 = QGridLayout(self.groupWD)
        self.gridLayout1.setObjectName("gridLayout1")
        self.gridLayout1.setContentsMargins(0, 0, 0, 0)

        self.browseWDLabel = QLabel(self)
        self.browseWDLabel.setObjectName("browseWDLabel")
        self.browseWDLabel.setText("Browse libary file...")
        self.browseWDLabel.setDisabled(True)
        self.browseWDLabel.setStyleSheet('''QLabel {
        padding: 1px 1px; /* Inner spacing */
        }''')
        sizePolicy1.setHeightForWidth(self.browseWDLabel.sizePolicy().hasHeightForWidth())
        self.browseWDLabel.setSizePolicy(sizePolicy1)
        self.browseWDLabel.setEnabled(False)
        self.browseWDLabel.setMinimumSize(QSize(400, 0))
        self.browseWDLabel.setFrameShape(QFrame.Box)
        self.gridLayout1.addWidget(self.browseWDLabel, 0, 0, 1, 2, Qt.AlignLeft)

        self.browseWDButton = QPushButton("  Browse", self)
        self.browseWDButton.setObjectName("browseWDButton")
        browseIcon = QPixmap(":/Browse.png")
        self.browseWDButton.setIcon(QIcon(browseIcon))
        self.browseWDButton.setIconSize(QSize(30, 30)) 
        sizePolicy1.setHeightForWidth(self.browseWDButton.sizePolicy().hasHeightForWidth())
        self.browseWDButton.setSizePolicy(sizePolicy1)
        self.browseWDButton.clicked.connect(self.browse)
        self.gridLayout1.addWidget(self.browseWDButton, 0, 1, 1, 1, Qt.AlignRight)

        work_space_layout.addWidget(self.groupWD)
        
        self.continueButton = QPushButton("Split Library", self)
        self.continueButton.setObjectName("continueButton")
        self.continueButton.setMinimumSize(QSize(150, 50))
        self.continueButton.clicked.connect(self.splitLib)
        work_space_layout.addWidget(self.continueButton)

        central_widget_layout.addLayout(work_space_layout)
        central_widget.setLayout(central_widget_layout)
        self.setCentralWidget(central_widget)

# Use the generated key to encrypt the text
    def browse(self):
        self.browseWDLabel.setText(QFileDialog.getOpenFileName(None, 'Open library file',"", "Standard Formats (*.sdf *.mol2)")[0])

    def splitLib(self):
        wdLabel = os.path.dirname(self.browseWDLabel.text())
        library = self.browseWDLabel.text()
        libsingle = os.path.join(wdLabel, "Library_found")
        ligs = os.path.join(wdLabel, "Ligands_extracted")
        ligCount = 0
        if not os.path.exists(libsingle):
            os.makedirs(libsingle)

        if not os.path.exists(ligs):
            os.makedirs(ligs)

        try:
            if library.endswith(".mol2"):
                with open(library, mode='r', encoding='utf-8') as f:
                    text = f.read()
                d = "@<TRIPOS>MOLECULE"
                splitted = [d+e for e in text.split(d) if e]
                for i in splitted:
                    ligCount += 1
                    test = i.splitlines()
                    name = test[1]
                    with open(os.path.join(ligs, name+".mol2"),'wb') as u:
                        file = "\n".join(test).decode('utf-8')
                        u.write(file)

            elif library.endswith(".sdf"):
                with open(library, mode='r', encoding='utf-8') as f:
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
                    with open(os.path.join(ligs, name+".sdf"),'wb') as u:
                        file = "".join(file[1:len(file)-1])
                        u.write(file.encode('utf-8'))

            shutil.move(library, libsingle)
            displayMessage(f"✅ Library extracted & moved to: {os.path.abspath(libsingle)}\n\n{ligCount} ligands found & moved to {os.path.abspath(ligs)}.")

        except Exception:
            displayError("❌ Error encountered, please check your working directory and try again.")

class ligandSplitter(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.setupUi(self)
    def center_above_screen(self):
        screen_geometry = QDesktopWidget().screenGeometry()
        window_geometry = self.frameGeometry()
        screen_center = screen_geometry.center()
        screen_center.setX(screen_center.x() + 200)
        screen_center.setY(screen_center.y() - 30)
        window_geometry.moveCenter(screen_center)
        self.move(window_geometry.topLeft())

    def __init__(self):
        super().__init__()
        self.setWindowIcon(QIcon(':/IDLogo.png'))
        self.setWindowTitle("IDv2 Ligand Splitter")
        self.resize(300, 300)
        self.is_dark_mode = False
        self.center_above_screen()
        self.setStyleSheet("""
        QMainWindow {
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:1, 
            stop:0 #F5F7FA, 
            stop:0.5 #E4E7EB, 
            stop:1 #D9DEE4); /* Light and minimalistic gradient */
            color: #000000; /* Default dark text */
            border: 1px solid #A9A9A9;
        }
        
        QLabel{
            background: transparent;
            font-size: 14px;
            font-family: "Arial";
        }
        
        QGroupBox{
            border:0;
            background: transparent;
        }

        QPushButton {
            font-size: 14px; /* Matches the font size */
            padding: 2px 8px; /* Inner spacing */
            border: 1px solid #A9A9A9; /* Default border */
            border-radius: 3px; /* Rounded corners */
            color: black; /* Default text color */
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, stop:0 #F2F2F2, stop:0.5 #EBEBEB, stop:0.51 #DDDDDD, stop:1 #CFCFCF); /* Subtle 3D gradient */
        }

        QPushButton:hover {
            border: 1px solid #A9A9A9; /* Light border on hover */
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, stop:0 #EAF6FD, stop:0.5 #D9F0FC, stop:0.51 #BEE6FD, stop:1 #A7D9F5); /* Light blue gradient on hover */
        }

        QPushButton:pressed {
            padding: 2px 7px 3px 9px; /* Slight padding shift for pressed effect */
            border: 1px solid #A9A9A9; /* Darker border when pressed */
            border-bottom: 0px; /* Simulates "pressing in" by hiding the bottom border */
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, stop:0 #E5F4FC, stop:0.5 #C4E5F6, stop:0.51 #98D1EF, stop:1 #68B3DB); /* Darker gradient when pressed */
        }
                            
        QPushButton:checked {
            padding: 2px 7px 3px 9px; /* Slight padding shift for pressed effect */
            border: 1px solid #A9A9A9; /* Darker border when pressed */
            border-bottom: 0px; /* Simulates "pressing in" by hiding the bottom border */
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, stop:0 #E5F4FC, stop:0.5 #C4E5F6, stop:0.51 #98D1EF, stop:1 #68B3DB); /* Darker gradient when pressed */
        }

        QPushButton#startButton{
            font-size: 35px;
        }
                            
    """
)
        sizePolicy = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sizePolicy)
        self.setMinimumSize(QSize(300, 300))
        self.setWindowFlags(Qt.FramelessWindowHint)
        central_widget = QWidget()
        central_widget.setObjectName("Container")
        self.title_bar = CustomTitleBar(self)

        central_widget_layout = QVBoxLayout()
        central_widget_layout.setContentsMargins(0, 0, 0, 0)
        central_widget_layout.setAlignment(Qt.AlignTop)

        central_widget_layout.addWidget(self.title_bar)

        work_space_layout = QVBoxLayout()
        work_space_layout.setContentsMargins(11, 11, 11, 11)

        sizePolicy1 = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Minimum)
        sizePolicy1.setHorizontalStretch(0)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy2 = QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)
        sizePolicy2.setHorizontalStretch(0)
        sizePolicy2.setVerticalStretch(0)

        self.introLabel = QLabel(self)
        self.introLabel.setObjectName("introLabel")
        self.introLabel.setTextFormat(Qt.RichText)
        self.introLabel.setText('''
<span style="font-size: 14px; font-family: Arial;"> This submodule is specifically designed to efficiently separate docked ligand files, into multiple distinct ligand poses. <br/>It ensures that each pose is extracted accurately for further analysis or visualization. </span>
        ''')
        self.introLabel.setWordWrap(True)
        sizePolicy2.setHeightForWidth(self.introLabel.sizePolicy().hasHeightForWidth())
        self.introLabel.setSizePolicy(sizePolicy2)
        self.introLabel.setAlignment(Qt.AlignJustify)

        work_space_layout.addWidget(self.introLabel)

        self.groupWD = QGroupBox(self)
        self.groupWD.setObjectName("groupWD")
        
        sizePolicy2.setHeightForWidth(self.groupWD.sizePolicy().hasHeightForWidth())
        self.groupWD.setSizePolicy(sizePolicy2)
        self.groupWD.setMinimumSize(QSize(280, 75))

        self.gridLayout1 = QGridLayout(self.groupWD)
        self.gridLayout1.setObjectName("gridLayout1")
        self.gridLayout1.setContentsMargins(0, 0, 0, 0)

        self.browseWDLabel = QLabel(self)
        self.browseWDLabel.setObjectName("browseWDLabel")
        self.browseWDLabel.setText("Browse docked ligand file...")
        self.browseWDLabel.setDisabled(True)
        self.browseWDLabel.setStyleSheet('''QLabel {
        padding: 1px 1px; /* Inner spacing */
        }''')
        sizePolicy1.setHeightForWidth(self.browseWDLabel.sizePolicy().hasHeightForWidth())
        self.browseWDLabel.setSizePolicy(sizePolicy1)
        self.browseWDLabel.setEnabled(False)
        self.browseWDLabel.setMinimumSize(QSize(400, 0))
        self.browseWDLabel.setFrameShape(QFrame.Box)
        self.gridLayout1.addWidget(self.browseWDLabel, 0, 0, 1, 2, Qt.AlignLeft)

        self.browseWDButton = QPushButton("  Browse", self)
        self.browseWDButton.setObjectName("browseWDButton")
        browseIcon = QPixmap(":/Browse.png")
        self.browseWDButton.setIcon(QIcon(browseIcon))
        self.browseWDButton.setIconSize(QSize(30, 30)) 
        sizePolicy1.setHeightForWidth(self.browseWDButton.sizePolicy().hasHeightForWidth())
        self.browseWDButton.setSizePolicy(sizePolicy1)
        self.browseWDButton.clicked.connect(self.browse)
        self.gridLayout1.addWidget(self.browseWDButton, 0, 1, 1, 1, Qt.AlignRight)

        work_space_layout.addWidget(self.groupWD)
        
        self.continueButton = QPushButton("Split Ligand", self)
        self.continueButton.setObjectName("continueButton")
        self.continueButton.setMinimumSize(QSize(150, 50))
        self.continueButton.clicked.connect(self.splitLig)
        work_space_layout.addWidget(self.continueButton)

        central_widget_layout.addLayout(work_space_layout)
        central_widget.setLayout(central_widget_layout)
        self.setCentralWidget(central_widget)

# Use the generated key to encrypt the text
    def browse(self):
        files, _ = QFileDialog.getOpenFileNames(
            self,
            'Select Vina Output Files',
            '',
            'Vina Output Files (*_out.pdbqt)'
        )
        if files:
            self.selected_files = files
            # show summary on label
            if len(files) == 1:
                text = os.path.basename(files[0])
            else:
                text = f"{len(files)} files selected"
            self.browseWDLabel.setText(text)

    def splitLig(self):
        try:
            original_file = self.browseWDLabel.text()
            target_dir = os.path.join(os.path.dirname(original_file), "Splitted files")

            if not os.path.exists(target_dir):
                os.makedirs(target_dir)

            with open(original_file, 'r') as f:
                content = f.read()

            models = content.split("MODEL")[1:]
            j = 1

            for model in models:
                # Generate a new file name without modifying the original file path
                base_name = re.sub('_out\.pdbqt$', '', os.path.basename(original_file))
                new_file_name = f"{base_name}_out_ligand_{j}.pdbqt"
                
                lines = model.split('\n')
                filtered_lines = []
                for line in lines:
                    if re.match(r"\s\d", line) or line.startswith("ENDMDL"):
                        continue
                    filtered_lines.append(line)
                
                final_content = '\n'.join(filtered_lines)
                with open(new_file_name, 'w') as out_file:
                    out_file.write(final_content)
                
                j += 1

            generated_files = glob(os.path.join(os.getcwd(), "*_out_ligand_*.pdbqt"))

            for file in generated_files:
                shutil.move(os.path.abspath(file), target_dir)

            displayMessage(f"✅ Ligand extracted into {len(models)} poses & moved to {target_dir}.")

        except Exception:
            displayError("❌ Error encountered, please check your working directory and try again.")

class topHits(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.setupUi(self)
    def center_above_screen(self):
        screen_geometry = QDesktopWidget().screenGeometry()
        window_geometry = self.frameGeometry()
        screen_center = screen_geometry.center()
        screen_center.setX(screen_center.x() + 200)
        screen_center.setY(screen_center.y() - 30)
        window_geometry.moveCenter(screen_center)
        self.move(window_geometry.topLeft())

    def __init__(self):
        super().__init__()
        self.setWindowIcon(QIcon(':/IDLogo.png'))
        self.setWindowTitle("IDv2 Top Hit Identifier")
        self.resize(400, 400)
        self.is_dark_mode = False
        self.center_above_screen()
        self.setStyleSheet("""
        QMainWindow {
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:1, 
            stop:0 #F5F7FA, 
            stop:0.5 #E4E7EB, 
            stop:1 #D9DEE4); /* Light and minimalistic gradient */
            color: #000000; /* Default dark text */
            border: 1px solid #A9A9A9;
        }
        
        QLabel{
            background: transparent;
            font-size: 14px;
            font-family: "Arial";
        }
        
        QGroupBox{
            border:0;
            background: transparent;
        }

        QPushButton {
            font-size: 14px; /* Matches the font size */
            padding: 2px 8px; /* Inner spacing */
            border: 1px solid #A9A9A9; /* Default border */
            border-radius: 3px; /* Rounded corners */
            color: black; /* Default text color */
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, stop:0 #F2F2F2, stop:0.5 #EBEBEB, stop:0.51 #DDDDDD, stop:1 #CFCFCF); /* Subtle 3D gradient */
        }

        QPushButton:hover {
            border: 1px solid #A9A9A9; /* Light border on hover */
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, stop:0 #EAF6FD, stop:0.5 #D9F0FC, stop:0.51 #BEE6FD, stop:1 #A7D9F5); /* Light blue gradient on hover */
        }

        QPushButton:pressed {
            padding: 2px 7px 3px 9px; /* Slight padding shift for pressed effect */
            border: 1px solid #A9A9A9; /* Darker border when pressed */
            border-bottom: 0px; /* Simulates "pressing in" by hiding the bottom border */
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, stop:0 #E5F4FC, stop:0.5 #C4E5F6, stop:0.51 #98D1EF, stop:1 #68B3DB); /* Darker gradient when pressed */
        }
                            
        QPushButton:checked {
            padding: 2px 7px 3px 9px; /* Slight padding shift for pressed effect */
            border: 1px solid #A9A9A9; /* Darker border when pressed */
            border-bottom: 0px; /* Simulates "pressing in" by hiding the bottom border */
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, stop:0 #E5F4FC, stop:0.5 #C4E5F6, stop:0.51 #98D1EF, stop:1 #68B3DB); /* Darker gradient when pressed */
        }

        QPushButton#startButton{
            font-size: 35px;
        }
                           
        QSpinBox {
            font-family: "Arial";
            font-size: 14px;
            padding: 2px 20px 2px 8px; 
            border: 1px solid #A9A9A9;
            border-radius: 3px;
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:1, 
                        stop:0 #F5F7FA, 
                        stop:0.5 #E4E7EB, 
                        stop:1 #D9DEE4);
        }
                        
    """
)
        sizePolicy = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sizePolicy)
        self.setMinimumSize(QSize(400, 400))
        self.setWindowFlags(Qt.FramelessWindowHint)
        central_widget = QWidget()
        central_widget.setObjectName("Container")
        self.title_bar = CustomTitleBar(self)

        central_widget_layout = QVBoxLayout()
        central_widget_layout.setContentsMargins(0, 0, 0, 0)
        central_widget_layout.setAlignment(Qt.AlignTop)

        central_widget_layout.addWidget(self.title_bar)

        work_space_layout = QVBoxLayout()
        work_space_layout.setContentsMargins(11, 11, 11, 11)

        sizePolicy1 = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Minimum)
        sizePolicy1.setHorizontalStretch(0)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy2 = QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)
        sizePolicy2.setHorizontalStretch(0)
        sizePolicy2.setVerticalStretch(0)
        sizePolicy3 = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy3.setHorizontalStretch(0)
        sizePolicy3.setVerticalStretch(0)

        self.introLabel = QLabel(self)
        self.introLabel.setObjectName("introLabel")
        self.introLabel.setTextFormat(Qt.RichText)
        self.introLabel.setText('''
<span style="font-size: 14px; font-family: Arial;">
The purpose of this submodule is to accurately extract and analyze the binding poses of docked ligands in PDBQT format. <br/>By providing a robust framework for precise pose extraction, it streamlines the workflow for detailed analysis and visualization. 
</span>
''')

        self.introLabel.setWordWrap(True)
        sizePolicy2.setHeightForWidth(self.introLabel.sizePolicy().hasHeightForWidth())
        self.introLabel.setSizePolicy(sizePolicy3)
        self.introLabel.setAlignment(Qt.AlignJustify)

        work_space_layout.addWidget(self.introLabel)

        self.toggleGroup = QGroupBox(self)
        self.toggleGroup.setObjectName("toggleGroup")
        self.horizontalLayout1 = QHBoxLayout(self.toggleGroup)

        self.confGroup = QGroupBox(self)
        self.confGroup.setObjectName("confGroup")
        self.horizontalLayout2 = QHBoxLayout(self.confGroup)

        self.confLabel = QLabel(self.confGroup)
        self.confLabel.setObjectName("confLabel")
        self.confLabel.setText("Single-target")
        sizePolicy1.setHeightForWidth(self.confLabel.sizePolicy().hasHeightForWidth())
        self.confLabel.setSizePolicy(sizePolicy1)
        self.confLabel.setStyleSheet("text-align:left;")
        self.horizontalLayout2.addWidget(self.confLabel)

        self.toggle_1 = SwitchControl(bg_color="#777777", circle_color="#DDD", active_color="#12c6c2", animation_curve=QEasingCurve.InOutCubic, animation_duration=100, checked=True, change_cursor=True)
        self.toggle_1.stateChanged.connect(self.on_switch1_toggled)
        self.horizontalLayout2.addWidget(self.toggle_1)

        self.horizontalLayout1.addWidget(self.confGroup)

        self.siteGroup = QGroupBox(self)
        self.siteGroup.setObjectName("confGroup")
        self.horizontalLayout3 = QHBoxLayout(self.siteGroup)

        self.siteLabel = QLabel(self.siteGroup)
        self.siteLabel.setObjectName("siteLabel")
        self.siteLabel.setText("Multi-target")
        sizePolicy2.setHeightForWidth(self.siteLabel.sizePolicy().hasHeightForWidth())
        self.siteLabel.setSizePolicy(sizePolicy2)
        self.siteLabel.setStyleSheet("text-align:left;")
        self.horizontalLayout3.addWidget(self.siteLabel)

        self.toggle_2 = SwitchControl(bg_color="#777777", circle_color="#DDD", active_color="#12c6c2", animation_curve=QEasingCurve.InOutCubic, animation_duration=100, checked=False, change_cursor=True)
        self.toggle_2.stateChanged.connect(self.on_switch2_toggled)
        self.horizontalLayout3.addWidget(self.toggle_2)

        self.horizontalLayout1.addWidget(self.siteGroup)

        work_space_layout.addWidget(self.toggleGroup)

        self.groupWD = QGroupBox(self)
        self.groupWD.setObjectName("groupWD")
        
        sizePolicy2.setHeightForWidth(self.groupWD.sizePolicy().hasHeightForWidth())
        self.groupWD.setSizePolicy(sizePolicy2)
        self.groupWD.setMinimumSize(QSize(380, 50))

        self.gridLayout1 = QGridLayout(self.groupWD)
        self.gridLayout1.setObjectName("gridLayout1")
        self.gridLayout1.setContentsMargins(0, 0, 0, 0)

        self.browseWDLabel = QLabel(self)
        self.browseWDLabel.setObjectName("browseWDLabel")
        self.browseWDLabel.setText("Browse IDv2 result folder...")
        self.browseWDLabel.setDisabled(True)
        self.browseWDLabel.setStyleSheet('''QLabel {
        padding: 1px 1px; /* Inner spacing */
        }''')
        sizePolicy1.setHeightForWidth(self.browseWDLabel.sizePolicy().hasHeightForWidth())
        self.browseWDLabel.setSizePolicy(sizePolicy1)
        self.browseWDLabel.setEnabled(False)
        self.browseWDLabel.setMinimumWidth(387)
        self.browseWDLabel.setFrameShape(QFrame.Box)
        self.gridLayout1.addWidget(self.browseWDLabel, 0, 0, 1, 2, Qt.AlignLeft)

        self.browseWDButton = QPushButton("  Browse", self)
        self.browseWDButton.setObjectName("browseWDButton")
        browseIcon = QPixmap(":/Browse.png")
        self.browseWDButton.setIcon(QIcon(browseIcon))
        self.browseWDButton.setIconSize(QSize(30, 30)) 
        sizePolicy1.setHeightForWidth(self.browseWDButton.sizePolicy().hasHeightForWidth())
        self.browseWDButton.setSizePolicy(sizePolicy1)
        self.browseWDButton.clicked.connect(self.browse)
        self.gridLayout1.addWidget(self.browseWDButton, 0, 1, 1, 1, Qt.AlignRight)

        work_space_layout.addWidget(self.groupWD)

        self.topHitsSpinBox = QSpinBox()
        self.topHitsSpinBox.setMinimum(1)
        self.topHitsSpinBox.setMaximum(5)
        self.topHitsSpinBox.setValue(1)
        self.topHitsSpinBox.setMinimumHeight(50)
        self.topHitsSpinBox.setObjectName("topHitsSpinBox")
        work_space_layout.addWidget(self.topHitsSpinBox)

        self.continueButton = QPushButton("Get the top hits!", self)
        self.continueButton.setObjectName("continueButton")
        self.continueButton.setMinimumSize(QSize(150, 50))
        self.continueButton.clicked.connect(self.hitIt)
        work_space_layout.addWidget(self.continueButton)

        central_widget_layout.addLayout(work_space_layout)
        central_widget.setLayout(central_widget_layout)
        self.setCentralWidget(central_widget)
    
    def findLigsinFolder(self, directory):
        """Return list of ligand files matching the pattern *_out.pdbqt."""
        return glob(os.path.join(directory, '*_out.pdbqt'))
    
    def browse(self):
        self.browseWDLabel.setText(QFileDialog.getExistingDirectory(None, 'Select the docking result folder'))
        results_dir = self.browseWDLabel.text()
        if self.toggle_1.isChecked():
            docked_dir = os.path.join(results_dir, "Docked Ligands")
            os.chdir(results_dir)
            ligands = self.findLigsinFolder(docked_dir)
        else:
            os.chdir(results_dir)
            subfolders = [entry.name for entry in os.scandir() if entry.is_dir()]
            for target in subfolders:
                if target != "Plots":
                    ligands = self.findLigsinFolder(os.path.join(results_dir, target))

        self.topHitsSpinBox.setMaximum(len(ligands))
        self.topHitsSpinBox.setValue(1)     

    def create_affinity_file_from_pdbqt(self, input_dir, pdbqt_files, multiTag):
        results = []
        for pdbqt_file in pdbqt_files:
            ligand_entry = {}
            if multiTag:
                ligand_entry.update({'Target': os.path.basename(input_dir)})
            ligand_entry.update({'Ligand': os.path.splitext(os.path.basename(pdbqt_file))[0].replace('_out', '')})
            with open(pdbqt_file, 'r') as f:
                content = f.readlines()
            energy_value = None
            # Compute torsional energy using the FIRST occurrence of TORSDOF
            torsional_energy = 'NA'
            try:
                tors_line = next(line for line in content if "TORSDOF" in line)
                # Assume the line format is: "TORSDOF   <number>"
                tors_value = float(tors_line.split()[1])
                torsional_energy = round(tors_value * 0.3113, 3)
            except (StopIteration, IndexError, ValueError):
                torsional_energy = 'NA'
            try:
                energy_line = next(line for line in content 
                                   if line.startswith('REMARK') and 'Binding energy' in line)
                energy_value = float(energy_line.split('=')[1].strip().split()[0])
            except (StopIteration, IndexError, ValueError):
                try:
                    energy_line = next(line for line in content if line.startswith('REMARK VINA RESULT'))
                    parts = energy_line.strip().split()
                    energy_value = float(parts[3])
                except (StopIteration, IndexError, ValueError):
                    show_toast([1000, f"Skipping {pdbqt_file}: No valid energy line found"])
                    continue
            nh_count = sum(1 for line in content if line.startswith(('ATOM', 'HETATM')) and ' H ' not in line)
            ligand_entry.update({
                'Binding Free Energy (kcal/mol)': energy_value,
                'Ligand Efficiency (kcal/mol/non-H atom)': round(-energy_value / nh_count, 3) if nh_count else 0,
                'pKi': round(-math.log10(math.exp(energy_value * 1000 / (1.986 * 298.15))), 3),
                'Torsional Energy (kcal/mol)': torsional_energy
            })           
            results.append(ligand_entry)
        
        return results

# Use the generated key to encrypt the text
    def hitIt(self):
        if self.toggle_1.isChecked():
            self.singleTarget()
        else:
            self.multiTarget()

    def singleTarget(self):
        results_dir = self.browseWDLabel.text()
        num_top_hits = int(self.topHitsSpinBox.cleanText())
        if results_dir == "Browse IDv2 result folder...":
            displayError("Please browse an appropriate IDv2 result folder, and try again...")
        else:
            try:
                docked_dir = os.path.join(results_dir, "Docked Ligands")
                pdbqt_files = self.findLigsinFolder(docked_dir)
                # Include multiTag if needed in your method signature (set to True/False as per your logic)
                affinity_results = self.create_affinity_file_from_pdbqt(results_dir, pdbqt_files, multiTag=False)
                                
                if not affinity_results:
                    displayError("No valid ligand results were found.")
                    return

                # Sort results by binding energy (assuming more negative is better)
                sorted_hits = sorted(affinity_results, key=lambda entry: float(entry['Binding Free Energy (kcal/mol)']))
                num_top_hits = min(num_top_hits, len(sorted_hits))

            # Select the top N ligands
                top_hits = sorted_hits[:num_top_hits]

                # OR, if you need the entire row, you might write out a single-row CSV
                csv_path = os.path.join(str(results_dir), f"Top_{num_top_hits}_Hits_IDv2_Results.csv")
                fieldnames = ['Ligand', 'Binding Free Energy (kcal/mol)', 'pKi', 
                            'Ligand Efficiency (kcal/mol/non-H atom)', 'Torsional Energy (kcal/mol)']
                with open(csv_path, 'w', newline='') as f:
                    writer = csv.DictWriter(f, fieldnames=fieldnames)
                    writer.writeheader()
                    for row in top_hits:
                        writer.writerow({k: v if not isinstance(v, float) else round(v, 3) for k, v in row.items()})

                # Display a message with only the numerical binding energy if desired
                displayMessage(f"Top {num_top_hits} Hits saved to {csv_path}")            
            except Exception:
                displayError("❌ Error encountered, please check your working directory and try again.")

    def multiTarget(self):
        results_dir = self.browseWDLabel.text()
        num_top_hits = int(self.topHitsSpinBox.cleanText())
        if results_dir == "Browse IDv2 result folder...":
            displayError("Please browse an appropriate IDv2 result folder, and try again...")
        else:
            try:
                os.chdir(results_dir)
                combined_top_hits = []  # To collect top hits from every target
                subfolders = [entry.name for entry in os.scandir() if entry.is_dir()]
                
                for target in subfolders:
                    if target == "Plots":
                        continue
                    target_path = os.path.join(results_dir, target)
                    pdbqt_files = self.findLigsinFolder(target_path)
                    # Ensure that the affinity file created includes the 'Target' field; 
                    # you may need to pass an appropriate flag (e.g., multiTag=True)
                    affinity_list = self.create_affinity_file_from_pdbqt(target_path, pdbqt_files, multiTag=True)
                    
                    if not affinity_list:
                        show_toast([1000, f"No valid ligand results found for target '{target}'."])
                        continue

                    # Sort the affinities for this target in ascending order (more negative is better)
                    sorted_hits = sorted(affinity_list, key=lambda entry: float(entry['Binding Free Energy (kcal/mol)']))

                    # Keep only the number of hits specified by the user for this target
                    num_top_hits = min(num_top_hits, len(sorted_hits))

                    # Select the top N ligands
                    top_hits = sorted_hits[:num_top_hits]
                    
                    combined_top_hits.extend(top_hits)
                
                if not combined_top_hits:
                    displayError("No valid ligand results were found across the targets.")
                    return

                # Write the combined top hits to CSV
                csv_path = os.path.join(results_dir, 'Top_Hits_IDv2_Results.csv')
                fieldnames = ['Target', 'Ligand', 'Binding Free Energy (kcal/mol)', 'pKi', 
                            'Ligand Efficiency (kcal/mol/non-H atom)', 'Torsional Energy (kcal/mol)']
                with open(csv_path, 'w', newline='') as f:
                    writer = csv.DictWriter(f, fieldnames=fieldnames)
                    writer.writeheader()
                    for row in combined_top_hits:
                        processed_row = {k: round(v, 3) if isinstance(v, float) else v for k, v in row.items()}
                        writer.writerow(processed_row)
                        
                displayMessage(f"Results saved to {csv_path}")
                
            except Exception:
                displayError("❌ Error encountered, please check your working directory and try again.")

    def on_switch1_toggled(self, state):
        # Turn off switch2 if switch1 is turned on
        if state==2:
            self.toggle_2.start_animation(0)
        else:
            self.toggle_2.start_animation(2)

    def on_switch2_toggled(self, state):
        # Turn off switch1 if switch2 is turned on
        if state==2:
            self.toggle_1.start_animation(0)
        else:
            self.toggle_1.start_animation(2)
