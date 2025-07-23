import warnings
warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=DeprecationWarning)

import re
import os
import random
import ID_utils
import prolif as plf
from rdkit import Chem
from PyQt5.QtGui import QIcon, QPixmap
from prolif.plotting.network import LigNetwork
from http.server import HTTPServer, SimpleHTTPRequestHandler
from PyQt5 import QtCore, QtWidgets, QtWebEngineWidgets
from PyQt5.QtCore import QSize, Qt, QThread, pyqtSignal, QProcess
from PyQt5.QtWidgets import (QFileDialog, QMainWindow, QDesktopWidget, QHBoxLayout, 
                             QSpacerItem, QSizePolicy, QWidget, QToolButton, QApplication, 
                             QGridLayout, QVBoxLayout, QGroupBox, QLabel, QPushButton, QFrame)


HOST = '0.0.0.0'
PORT = random.randint(1024, 65535)

class CORSRequestHandler(SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Access-Control-Allow-Methods', '*')
        self.send_header('Access-Control-Allow-Headers', '*')
        self.send_header('Cache-Control', 'no-store, no-cache, must-revalidate')
        return super(CORSRequestHandler, self).end_headers()

class HttpDaemon(QThread):
    def run(self):
        self._server = HTTPServer((HOST, PORT), CORSRequestHandler)
        self._server.serve_forever()

    def stop(self):
        self._server.shutdown()
        self._server.socket.close()
        self.wait()

class Worker(QThread):
    finished = pyqtSignal()
    errorOccurred = pyqtSignal(str)
    progress = pyqtSignal(str)

    def __init__(self, receptor_file, ligand_file, results_folder):
        super().__init__()
        self.receptor_file = receptor_file
        self.ligand_file = ligand_file
        self.results_folder = results_folder

    def run(self):
        try:
            # Create necessary directories
            ligand_pdbs_folder = os.path.join(self.results_folder, "ligand_pdbs")
            os.makedirs(ligand_pdbs_folder, exist_ok=True)
            
            # Convert PDBQT to PDB using QProcess
            base_name = os.path.splitext(os.path.basename(self.ligand_file))[0]
            converted_ligand_path = os.path.join(ligand_pdbs_folder, f"{base_name}.pdb")
            
            exe = "obabel"
            arg = ["-ipdbqt", self.ligand_file, "-xh", "-opdb", "-O", converted_ligand_path]
            
            process = QProcess()
            process.start(exe, arg)
            process.waitForFinished(-1)

            # Process molecules
            molPro = Chem.MolFromPDBFile(self.receptor_file, removeHs=False, sanitize=False)
            prot = plf.Molecule(molPro)
            
            molLig = Chem.MolFromPDBFile(converted_ligand_path, removeHs=False, sanitize=False)
            lig = plf.Molecule(molLig)
            
            fp = plf.Fingerprint(parameters={"HBAcceptor": {"DHA_angle": (120, 180)}})
            fp.run_from_iterable([lig], prot, n_jobs=1)
            df = fp.to_dataframe()

            if df.empty:
                raise RuntimeError("❌ No interactions detected between your ligand and protein.")
            
            # Generate visualization
            test = LigNetwork.from_fingerprint(fp, lig, kind="frame", frame=0, display_all=True)
            html_filename = os.path.join(self.results_folder, 
                f"{base_name}_{os.path.basename(self.receptor_file).split('.')[0]}.html")
            test.save(html_filename)
            
            # Modify HTML content
            with open(html_filename, 'r', encoding='utf-8') as f:
                html_content = f.read()
                
            new_line = (
                "ifp = drawGraph('mynetwork', nodes, edges, "
                '{"width": "100%", "height": "80%", "nodes": {"font": {"size": 20}}, '
                '"physics": {"barnesHut": {"avoidOverlap": 0.8, "springConstant": 0.1}}});'
            )
            
            pattern = r"ifp = drawGraph\('mynetwork', nodes, edges, .*?\);"
            modified_content = re.sub(pattern, new_line, html_content)
            
            with open(html_filename, 'w', encoding='utf-8') as f:
                f.write(modified_content)
                
            self.finished.emit()
            
        except Exception as e:
            self.errorOccurred.emit(str(e))

class CustomTitleBar(QWidget):
    def __init__(self, parent):
        super().__init__(parent)
        self.initial_pos = None
        title_bar_layout = QHBoxLayout(self)
        title_bar_layout.setContentsMargins(1, 1, 1, 1)
        title_bar_layout.setSpacing(2)

        # Spacer to push buttons to the right
        title_bar_layout.addSpacerItem(QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum))

        # Minimize button
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

class WebEngineWindow(QMainWindow):
    def __init__(self, html_content, parent=None):
        super().__init__(parent)
        self.setWindowTitle("IDv2 2D Interaction Visualizer")
        self.resize(800, 600)
        
        # Create a QWebEngineView widget
        self.view = QtWebEngineWidgets.QWebEngineView(self)
        self.setCentralWidget(self.view)
        
        # Set local content to access remote URLs if needed.
        settings = QtWebEngineWidgets.QWebEngineSettings.defaultSettings()
        settings.setAttribute(QtWebEngineWidgets.QWebEngineSettings.LocalContentCanAccessRemoteUrls, True)
        
        # Set HTML content directly into the view
        self.view.setHtml(html_content)
        
        # Connect download request signal if needed.
        self.view.page().profile().downloadRequested.connect(parent.on_downloadRequested)

class idTwoDInteract(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowIcon(QIcon(':/IDLogo.png'))
        self.setWindowTitle("IDv2 2D Interaction Visualizer")
        self.resize(300, 400)
        self.is_dark_mode = False
        self.center_above_screen()
        self.setStyleSheet("""
            QMainWindow {
                background: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:1, 
                stop:0 #F5F7FA, 
                stop:0.5 #E4E7EB, 
                stop:1 #D9DEE4);
                color: #000000;
                border: 1px solid #A9A9A9;
            }
            
            QLabel{
                background: transparent;
                font-size: 14px;
                font-family: "Arial";
            }
                               
            QLabel#browseWDLabel{
                background: transparent;
                font-size: 14px;
                border: 1px solid #A9A9A9;
                border-radius: 3px;
            }
            
            QGroupBox{
                border:0;
                background: transparent;
            }
    
            QPushButton {
                font-size: 14px;
                padding: 2px 8px;
                border: 1px solid #A9A9A9;
                border-radius: 3px;
                color: black;
                background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, stop:0 #F2F2F2, stop:0.5 #EBEBEB, stop:0.51 #DDDDDD, stop:1 #CFCFCF);
            }
    
            QPushButton:hover {
                border: 1px solid #A9A9A9;
                background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, stop:0 #EAF6FD, stop:0.5 #D9F0FC, stop:0.51 #BEE6FD, stop:1 #A7D9F5);
            }
    
            QPushButton:pressed, QPushButton:checked {
                padding: 2px 7px 3px 9px;
                border: 1px solid #A9A9A9;
                border-bottom: 0px;
                background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, stop:0 #E5F4FC, stop:0.5 #C4E5F6, stop:0.51 #98D1EF, stop:1 #68B3DB);
            }
                               
            QPushButton#startButton{
                font-size: 20px;
            }
        """)

        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy2 = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Minimum)
        sizePolicy3 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Minimum)
        sizePolicy4 = QSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.MinimumExpanding)
        self.setSizePolicy(sizePolicy)
        self.setMinimumSize(QSize(300, 400))
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
        
        # Intro label
        self.introLabel = QLabel(self)
        self.introLabel.setObjectName("introLabel")
        self.introLabel.setTextFormat(Qt.RichText)
        self.introLabel.setText('''
        <span style="font-size: 13px; font-family: Arial;"> 
        This module provides an interactive 2D view for protein-ligand complexes in PDB/PDBQT formats,  
        enabling users to explore molecular interactions with clarity!<br/>
        Its intuitive design allows both beginners and experts to navigate molecular data easily.  
        <br/>
        Additionally, it offers customizable visualization options to enhance analysis and usability.  
        </span>
        ''')
        self.introLabel.setWordWrap(True)
        self.introLabel.setSizePolicy(sizePolicy)
        self.introLabel.setAlignment(Qt.AlignJustify)
        work_space_layout.addWidget(self.introLabel)
        
        # Target file group
        self.groupTargetWD = QGroupBox(self)
        self.groupTargetWD.setObjectName("groupTargetWD")
        self.groupTargetWD.setSizePolicy(sizePolicy4)
        
        self.gridLayout2 = QGridLayout(self.groupTargetWD)
        self.gridLayout2.setContentsMargins(0, 0, 0, 0)
        
        self.targetLabel = QLabel(self.groupTargetWD)
        self.targetLabel.setObjectName("targetLabel")
        self.targetLabel.setText("Browse target file (PDB)...")
        self.targetLabel.setSizePolicy(sizePolicy4)
        self.targetLabel.setMinimumWidth(300)
        self.targetLabel.setDisabled(True)
        self.targetLabel.setFrameShape(QFrame.Box)
        self.gridLayout2.addWidget(self.targetLabel, 0, 0, 1, 2, Qt.AlignLeft)
        
        self.browseTargetButton = QPushButton("  Browse", self.groupTargetWD)
        self.browseTargetButton.setObjectName("browseTargetButton")
        browseIcon = QPixmap(":/Browse.png")
        self.browseTargetButton.setIcon(QIcon(browseIcon))
        self.browseTargetButton.setIconSize(QSize(30, 30))
        self.browseTargetButton.setSizePolicy(sizePolicy)
        self.browseTargetButton.clicked.connect(self.browseTarget)
        self.gridLayout2.addWidget(self.browseTargetButton, 0, 1, 1, 1, Qt.AlignRight)
        
        work_space_layout.addWidget(self.groupTargetWD)
        
        # Ligand file group – updated for PDBQT files
        self.groupLigandWD = QGroupBox(self)
        self.groupLigandWD.setObjectName("groupLigandWD")
        self.groupLigandWD.setSizePolicy(sizePolicy4)
        
        self.gridLayout3 = QGridLayout(self.groupLigandWD)
        self.gridLayout3.setContentsMargins(0, 0, 0, 0)
        
        self.ligandLabel = QLabel(self.groupLigandWD)
        self.ligandLabel.setObjectName("ligandLabel")
        # Update label text to reflect that we expect a PDBQT file
        self.ligandLabel.setText("Browse ligand file (PDBQT)...")
        self.ligandLabel.setSizePolicy(sizePolicy4)
        self.ligandLabel.setDisabled(True)
        self.ligandLabel.setFrameShape(QFrame.Box)
        self.ligandLabel.setMinimumWidth(300)
        self.gridLayout3.addWidget(self.ligandLabel, 0, 0, 1, 2, Qt.AlignLeft)
        
        self.browseLigandButton = QPushButton("  Browse", self.groupLigandWD)
        self.browseLigandButton.setObjectName("browseLigandButton")
        browseIcon = QPixmap(":/Browse.png")
        self.browseLigandButton.setIcon(QIcon(browseIcon))
        self.browseLigandButton.setIconSize(QSize(30, 30))
        self.browseLigandButton.setSizePolicy(sizePolicy)
        self.browseLigandButton.setMinimumWidth(100)
        self.browseLigandButton.clicked.connect(self.browseLigand)
        self.gridLayout3.addWidget(self.browseLigandButton, 0, 1, 1, 1, Qt.AlignRight)
        
        work_space_layout.addWidget(self.groupLigandWD)
        
        # Start button
        self.startButton = QPushButton(self)
        self.startButton.setObjectName("startButton")
        self.startButton.setText("V I S U A L I Z E!")
        self.startButton.setSizePolicy(sizePolicy3)
        self.startButton.setMinimumHeight(50)
        self.startButton.setMaximumHeight(50)
        self.startButton.clicked.connect(self.visualize)
        work_space_layout.addWidget(self.startButton, Qt.AlignCenter)
        
        # Finalize main layout
        central_widget_layout.addLayout(work_space_layout)
        central_widget.setLayout(central_widget_layout)
        self.setCentralWidget(central_widget)

    def center_above_screen(self):
        screen_geometry = QDesktopWidget().screenGeometry()
        window_geometry = self.frameGeometry()
        screen_center = screen_geometry.center()
        window_geometry.moveCenter(screen_center)
        self.move(window_geometry.topLeft())
        
    def on_downloadRequested(self, download):
        directory = str(QFileDialog.getExistingDirectory(None, "Select Directory"))
        filename = QtCore.QFileInfo(download.path()).fileName()
        download.setPath(QtCore.QDir(directory).filePath(filename))
        download.accept()

    def browseLigand(self):
        # Updated file filter for PDBQT files
        filename = QFileDialog.getOpenFileName(None, 'Open ligand file', "", "PDBQT File (*.pdbqt)")[0]
        if filename:
            self.ligandLabel.setText(filename)

    def browseTarget(self):
        filename = QFileDialog.getOpenFileName(None, 'Open target file', "", "PDB File (*.pdb)")[0]
        if filename:
            self.targetLabel.setText(filename)

    def visualize(self):
        if (self.targetLabel.text() == "Browse target file..." or 
            self.ligandLabel.text() == "Browse ligand file (PDBQT)..."):
            ID_utils.displayError("❌ Error: Please select both files first.")
            return
        
        # Start the HTTP server in the background
        self.httpserver = HttpDaemon()
        self.httpserver.start()

        self.worker = Worker(
            receptor_file=self.targetLabel.text(),
            ligand_file=self.ligandLabel.text(),
            results_folder=os.path.join(
                os.path.dirname(self.targetLabel.text()),
                "2D_interaction_results"
            )
        )
        self.worker.finished.connect(self.show_visualization)
        self.worker.errorOccurred.connect(self.showError)
        self.startButton.setText("Loading...")
        self.startButton.setDisabled(True)
        self.worker.start()

    def showError(self, message):
        ID_utils.displayError(message)

    def show_visualization(self):
        # Find the generated HTML file
        receptor_name = os.path.basename(self.targetLabel.text()).split('.')[0]
        ligand_name = os.path.basename(self.ligandLabel.text()).split('.')[0]
        html_path = os.path.join(
            os.path.dirname(self.targetLabel.text()),
            "2D_interaction_results",
            f"{ligand_name}_{receptor_name}.html"
        )
        
        if os.path.exists(html_path):
            with open(html_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
                
            self.webWindow = WebEngineWindow(html_content, self)
            self.webWindow.show()
            self.startButton.setDisabled(False)
            self.startButton.setText("V I S U A L I Z E!")
        else:
            ID_utils.displayError("❌ Error: Failed to generate visualization file")

if hasattr(Qt, 'AA_EnableHighDpiScaling'):
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
if hasattr(Qt, 'AA_UseHighDpiPixmaps'):
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)
