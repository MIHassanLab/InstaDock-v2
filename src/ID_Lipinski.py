import os
import sys
import csv
import ID_utils
import tempfile
import IDv2Resource
from glob import glob 
from rdkit import Chem
from PyQt5 import QtWidgets, QtCore
from pyqttoast import Toast, ToastPreset
from rdkit.Chem import Descriptors, Lipinski
from PyQt5.QtGui import QIcon, QPixmap
from PyQt5.QtCore import QSize, Qt, QProcess
from PyQt5.QtWidgets import (
    QPushButton, QGridLayout, QHBoxLayout, QLabel, QMainWindow, 
    QVBoxLayout, QWidget, QSpacerItem, QSizePolicy, QGroupBox,
    QDesktopWidget, QToolButton, QFrame, QFileDialog, QProgressBar)

LIGAND_TAGS = [' UNL ', ' UNK ', ' LIG ', ' <0> ', ' DRG ', ' INH ', ' NAG ', ' SO4 ', ' ATP ', ' ADP ', ' AMP ', ' HEM ', ' FMN ', ' FAD ', ' NAD ', ' GDP ', ' GTP ', ' SAM ']

STD_RESIDUES = [' ALA ', ' ARG ', ' ASN ', ' ASP ', ' CYS ', ' GLN ', ' GLU ', ' GLY ', ' HIS ', ' ILE ', ' LEU ', ' LYS ', ' MET ', ' PHE ', ' PRO ', ' SER ', ' THR ', ' TRP ', ' TYR ', ' VAL ']

def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    if hasattr(sys, '_MEIPASS'):
        # When bundled by PyInstaller, _MEIPASS contains the temporary extraction path.
        return os.path.join(sys._MEIPASS, relative_path)
    return os.path.join(os.path.abspath("."), relative_path)

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
        title_bar_layout.addWidget(self.min_button)

        # Close button
        self.close_button = QToolButton(self)
        self.close_button.setText("\u00D7")
        self.close_button.setStyleSheet("color:darkred;")
        self.close_button.setFixedSize(QSize(24, 24))
        self.close_button.clicked.connect(self.window().close)
        self.close_button.setFocusPolicy(Qt.NoFocus)
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
        except Exception as e:
            pass
        super().mouseMoveEvent(event)
        event.accept()

    def mouseReleaseEvent(self, event):
        self.initial_pos = None
        super().mouseReleaseEvent(event)
        event.accept()

class ro5(QMainWindow):
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
        self.setWindowTitle("IDv2 Lipinski Calculator")
        self.resize(320, 300)
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
        
        QProgressBar
        { 
            border: 2px solid #37cbe6; 
            border-radius: 5px; 
            text-align: center
        } 
        
        QProgressBar::chunk 
        { 
            background-color: lightblue; 
            width: 2px; 
            margin: 0.5px; 
            border-radius: 1.5px; 
        }
                            
        """
        )
        sizePolicy = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sizePolicy)
        self.setMinimumSize(QSize(320, 300))
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
        sizePolicy3 = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        sizePolicy3.setHorizontalStretch(0)
        sizePolicy3.setVerticalStretch(0)
        sizePolicy4 = QSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.MinimumExpanding)
        sizePolicy4.setHorizontalStretch(0)
        sizePolicy4.setVerticalStretch(0)

        self.introLabel = QLabel(self)
        self.introLabel.setObjectName("introLabel")
        self.introLabel.setTextFormat(Qt.RichText)
        self.introLabel.setText('''
<span style="font-size: 14px; font-family: Arial;"> 
This submodule is designed to provide a fast and accurate filter for Lipinski's Rule of Five, a set of guidelines used to evaluate the drug-likeness of given ligands.
</span>
                                ''')
        self.introLabel.setWordWrap(True)
        sizePolicy2.setHeightForWidth(self.introLabel.sizePolicy().hasHeightForWidth())
        self.introLabel.setSizePolicy(sizePolicy2)
        self.introLabel.setMinimumHeight(150)
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
        self.browseWDLabel.setText("Browse ligand directory...")
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

        self.progress_bar = QProgressBar()
        self.progress_bar.setTextVisible(False)
        self.progress_bar.setMaximumSize(self.width(), 20)
        work_space_layout.addWidget(self.progress_bar)
        
        self.startButton = QPushButton("E V A L U A T E", self)
        self.startButton.setObjectName("continueButton")
        self.startButton.setMinimumSize(QSize(150, 50))
        self.startButton.clicked.connect(self.startPred)
        work_space_layout.addWidget(self.startButton)

        central_widget_layout.addLayout(work_space_layout)
        central_widget.setLayout(central_widget_layout)
        self.setCentralWidget(central_widget)

    def browse(self):
        self.browseWDLabel.setText(QFileDialog.getExistingDirectory(None, 'Select Ligand directory'))

    def findligs(self, path):
        out_ligands = []
        unprocessed_ligands = []
        ligands = []
        cd = "obabel"
        self.process = QProcess()
        for file in glob(path + "/*.pdbqt"):
            if file.endswith("out.pdbqt"):
                out_ligands.append(file)
            unprocessed_ligands.append(file)
        
        if len(out_ligands) == 0 and len(unprocessed_ligands) == 0:
            ligands = None
        else:
            for ligand in unprocessed_ligands:
                with tempfile.NamedTemporaryFile(suffix=".mol2", delete=False) as tmp_mol2:
                    tmp_name = tmp_mol2.name
                try:
                    inputform = "-ipdbqt"
                    args_to_mol2 = [inputform, ligand, "-omol2", "-O", tmp_name]
                    self.process.start(cd, args_to_mol2)
                    self.process.waitForFinished(-1)
                    
                    inputform2 = "-imol2"
                    args_to_smi = [inputform2, tmp_name, "-osmi"]
                    self.process.start(cd, args_to_smi)
                    self.process.waitForFinished(-1)

                    output = self.process.readAllStandardOutput().data().decode().strip()

                    if ligand in out_ligands:
                        smiles = output.split("\n")[0].split("\t")[0]
                    else:
                        smiles = output.split("\t")[0]

                    ligands.append(smiles)

                finally:
                    if os.path.exists(tmp_name):
                        os.remove(tmp_name)

        return ligands    

    def startPred(self):
        if self.browseWDLabel.text() == "":
            ID_utils.displayError("❌ Please select a ligand directory first.")
        else:
            try:
                ligands = self.findligs(self.browseWDLabel.text())
                
                if ligands == None:
                    ID_utils.displayError("❌ No ligands found in the selected directory.")
                else:
                    lenLigand = len(ligands)
                    output_file = os.path.join(self.browseWDLabel.text(), "Lipinski_results.csv")
                    results = []
                    for i, smiles in enumerate(ligands): 
                        molecule = Chem.MolFromSmiles(smiles) 
                        self.update_progress_bar(int(((i + 1) / lenLigand) * 100))
                        if molecule: 
                            molecularWeight = Descriptors.MolWt(molecule) 
                            logP = Descriptors.MolLogP(molecule) 
                            hydrogenDonors = Lipinski.NumHDonors(molecule) 
                            hydrogenAcceptors = Lipinski.NumHAcceptors(molecule) 
                            passesLipinski = not (molecularWeight <= 500 and logP <= 5 and hydrogenDonors <= 5 and hydrogenAcceptors <= 10) # Showing violation rather than passing
                            results.append((smiles, molecularWeight, logP, hydrogenDonors, hydrogenAcceptors, passesLipinski)) 
                        else: 
                            results.append((smiles, None, None, None, None, False)) 
                    
                    with open(output_file, "w", newline="") as f: 
                        writer = csv.writer(f) 
                        writer.writerow(["SMILES", "Molecular Weight", "LogP", "HBD", "HBA", "Lipinski Violation"]) 
                        writer.writerows(results)

                ID_utils.displayMessage(f"✅ Lipinski Filter applied to {lenLigand} ligand(s), saved to {os.path.abspath(output_file)}.")
            except Exception as e:
                ID_utils.displayError("❌ Error encountered, please check your files and try again.")

    def update_progress_bar(self, val):
        self.progress_bar.setValue(val)

if hasattr(QtCore.Qt, 'AA_EnableHighDpiScaling'):
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
if hasattr(QtCore.Qt, 'AA_UseHighDpiPixmaps'):
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)