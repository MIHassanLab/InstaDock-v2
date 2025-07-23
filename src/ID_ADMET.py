import os
import re
import io
import sys
import random
import shutil
import base64
import zipfile
import hashlib
import ID_utils
import tempfile
import IDv2Resource
from glob import glob 
from rdkit import Chem
from pathlib import Path
from PyQt5 import QtWidgets, QtCore
from rdkit.Chem import AllChem, Draw
from pyqttoast import Toast, ToastPreset
from ID_Predict import predictADMET
from PyQt5.QtGui import QIcon, QPixmap, QFont
from jinja2 import Environment, FileSystemLoader
from PyQt5 import  QtWebEngineWidgets, QtNetwork
from http.server import HTTPServer, SimpleHTTPRequestHandler
from PyQt5.QtCore import (
    QSize, QSize, Qt, QThread, QProcess, QUrl, QTimer, QBuffer, 
    QByteArray, QIODevice
)
from PyQt5.QtWidgets import (
    QPushButton, QGridLayout, QHBoxLayout, QLabel, QMainWindow, 
    QVBoxLayout, QWidget, QSpacerItem, QSizePolicy, QGroupBox,
    QDesktopWidget, QToolButton, QFrame, QFileDialog, QProgressBar
)

HOST = '0.0.0.0'
PORT = random.randint(1024, 65535)

LIGAND_TAGS = [' UNL ', ' UNK ', ' LIG ', ' <0> ', ' DRG ', ' INH ', ' NAG ', ' SO4 ', ' ATP ', ' ADP ', ' AMP ', ' HEM ', ' FMN ', ' FAD ', ' NAD ', ' GDP ', ' GTP ', ' SAM ']

STD_RESIDUES = [' ALA ', ' ARG ', ' ASN ', ' ASP ', ' CYS ', ' GLN ', ' GLU ', ' GLY ', ' HIS ', ' ILE ', ' LEU ', ' LYS ', ' MET ', ' PHE ', ' PRO ', ' SER ', ' THR ', ' TRP ', ' TYR ', ' VAL ']

def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    if hasattr(sys, '_MEIPASS'):
        # When bundled by PyInstaller, _MEIPASS contains the temporary extraction path.
        return os.path.join(sys._MEIPASS, relative_path)
    return os.path.join(os.path.abspath("."), relative_path)

def get_resource_data_url(qt_resource_path):
    """Convert a Qt resource path (e.g., ':/error_placeholder.png') to a base64 data URL."""
    # Load the Qt resource directly
    pixmap = QPixmap(qt_resource_path)
    if pixmap.isNull():
        print(f"Error: Could not load resource {qt_resource_path}")
        return "data:image/png;base64,"  # Fallback empty image
    
    # Convert QPixmap to PNG bytes
    byte_array = QByteArray()
    buffer = QBuffer(byte_array)
    buffer.open(QIODevice.WriteOnly)
    pixmap.save(buffer, "PNG")
    buffer.close()
    
    # Encode to base64 data URL
    base64_data = byte_array.toBase64().data().decode()
    return f"data:image/png;base64,{base64_data}"

class CORSRequestHandler (SimpleHTTPRequestHandler):
    def end_headers (self):
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
        self.setWindowTitle("IDv2 ADMET Predictions")
        self.resize(800, 600)
        self.setWindowIcon(QIcon(':/IDLogo.png'))
        
        # Create a QWebEngineView widget
        self.view = QtWebEngineWidgets.QWebEngineView(self)
        self.setCentralWidget(self.view)
        
        # Set local content to access remote URLs if needed.
        settings = QtWebEngineWidgets.QWebEngineSettings.defaultSettings()
        settings.setAttribute(QtWebEngineWidgets.QWebEngineSettings.LocalContentCanAccessRemoteUrls, True)
        settings.setAttribute(QtWebEngineWidgets.QWebEngineSettings.LocalContentCanAccessFileUrls, True)

class admet(QMainWindow):
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
        self.setWindowTitle("IDv2 ADMET Predictor")
        self.resize(350, 370)
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
        self.setMinimumSize(QSize(350, 370))
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
<span style="font-size:14px; font-family:Arial;">
  This submodule delivers rapid, precise predictions for key physicochemical descriptors, comprehensive ADME profiles (absorption, distribution, metabolism, excretion), toxicity endpoints, and standard drug-likeness rules.  
  All in one, <b>30+</b> validated properties at your fingertips.
</span>
                                ''')
        self.introLabel.setWordWrap(True)
        sizePolicy2.setHeightForWidth(self.introLabel.sizePolicy().hasHeightForWidth())
        self.introLabel.setSizePolicy(sizePolicy2)
        self.introLabel.setMinimumHeight(175)
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
        
        self.startButton = QPushButton("Predict ADMET", self)
        self.startButton.setObjectName("continueButton")
        self.startButton.setMinimumSize(QSize(150, 50))
        self.startButton.clicked.connect(self.startPred)
        work_space_layout.addWidget(self.startButton)

        self.groupDownload = QGroupBox(self)
        self.groupDownload.setObjectName("groupDownload")
        sizePolicy3.setHeightForWidth(self.groupDownload.sizePolicy().hasHeightForWidth())
        self.groupDownload.setSizePolicy(sizePolicy3)
        self.groupDownload.setMinimumSize(QSize(280, 50))

        self.horizontalLayout = QHBoxLayout(self.groupDownload)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)

        self.downloadLabel = QLabel(self)
        self.downloadLabel.setObjectName("downloadLabel")
        self.downloadLabel.setText("To download AI models click this button  ➜➜➜")
        self.downloadLabel.setWordWrap(True)
        self.downloadLabel.setDisabled(True)
        self.horizontalLayout.addWidget(self.downloadLabel)
        
        self.downloadButton = QPushButton("  Download", self)
        self.downloadButton.setObjectName("downloadButton")
        self.downloadButton.setDisabled(False)
        downloadIcon = QPixmap(":/Download.png")
        self.downloadButton.setIcon(QIcon(downloadIcon))
        self.downloadButton.setIconSize(QSize(20, 20)) 
        self.downloadButton.clicked.connect(self.download)
        self.horizontalLayout.addWidget(self.downloadButton)

        work_space_layout.addWidget(self.groupDownload, Qt.AlignCenter)

        central_widget_layout.addLayout(work_space_layout)
        central_widget.setLayout(central_widget_layout)
        self.setCentralWidget(central_widget)

    def browse(self):
        self.browseWDLabel.setText(QFileDialog.getExistingDirectory(None, 'Select Ligand directory'))

    def get_models_path(self):
        if sys.platform.startswith("win"):
            base_dir = os.getenv("LOCALAPPDATA", str(Path.home()))
        else:
            base_dir = str(Path.home() / ".local" / "share")
        return Path(base_dir) / "IDv2Models"

    def _compute_file_hash(self, file_path):
        h = hashlib.sha256()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                h.update(chunk)
        return h.hexdigest()

    def _load_reference_hashes(self, ref_file_path):
        ref = {}
        if not os.path.exists(ref_file_path):
            return None
        with open(ref_file_path, "r") as f:
            for line in f:
                if "::" in line:
                    rel, hsh = line.strip().split("::", 1)
                    ref[rel] = hsh
        return ref

    def _verify_folder_integrity(self, models_path, reference_hashes):
        models_path = Path(models_path)
        local = {}
        for root, _, files in os.walk(models_path):
            for fn in files:
                p = Path(root) / fn
                rel = str(p.relative_to(models_path))
                local[rel] = self._compute_file_hash(str(p))

        # any missing/mismatch?
        for rel, expected in reference_hashes.items():
            if rel not in local or local[rel] != expected:
                return False
        # any extras?
        for rel in local:
            if rel not in reference_hashes:
                return False
        return True

    def download(self):
        models_path = self.get_models_path()

        # 1) If folder exists & non‑empty, verify its hashes
        if models_path.exists() and any(models_path.iterdir()):
            ref_file = resource_path("idv2_hashes.txt")
            ref_hashes = self._load_reference_hashes(ref_file)
            if ref_hashes and self._verify_folder_integrity(models_path, ref_hashes):
                self.show_toast([3000, "✅ Models already present and verified."])
                return
            # bad or tampered: blow it away
            shutil.rmtree(str(models_path))

        # 2) Prepare for download
        models_path.mkdir(parents=True, exist_ok=True)
        self.downloadButton.setDisabled(True)

        # 3) Build request with a “real” UA
		cfg = resource_path("models_config.json")
        with open(cfg, "r") as f:
            testURL = json.load(f)["admet_zip_url"]
			
        url = QtCore.QUrl(testURL)
        req = QtNetwork.QNetworkRequest(url)
        req.setRawHeader(b"User-Agent",
                         b"Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
                         b"AppleWebKit/537.36 (KHTML, like Gecko) "
                         b"Chrome/100.0.4896.127 Safari/537.36")

        # 4) Kick off download
        self.net_manager = QtNetwork.QNetworkAccessManager(self)
        self.reply = self.net_manager.get(req)

        # 5) Show progress dialog
        prog = QtWidgets.QProgressDialog("Downloading AI models...", "Cancel", 0, 100, self)
        prog.setWindowModality(QtCore.Qt.WindowModal)
        prog.setWindowTitle("IDv2 ADMET Download")
        prog.setStyleSheet("QLabel { font-weight: bold; }")
        prog.setAutoClose(False)
        prog.show()

        # 6) Wire up signals
        self.reply.downloadProgress.connect(
            lambda rec, tot: prog.setValue(int(rec * 100 / tot)) if tot else None
        )
        prog.canceled.connect(self.reply.abort)
        self.reply.finished.connect(lambda: self._on_download_finished(models_path, prog))

    def _on_download_finished(self, models_path, prog):
        # 1) Check for network errors
        err = self.reply.error()
        if err != QtNetwork.QNetworkReply.NoError:
            prog.close()
            ID_utils.displayError(f"❌ Download error: {self.reply.errorString()}")
            self.downloadButton.setDisabled(False)
            return

        # 2) Extract, stripping any top‑level folders
        prog.setLabelText("Extracting models…")
        raw = self.reply.readAll().data()
        try:
            with zipfile.ZipFile(io.BytesIO(raw)) as zf:
                files = [m for m in zf.infolist() if not m.is_dir()]
                # find common root prefix
                paths = [f.filename.split('/') for f in files]
                common = []
                for segs in zip(*paths):
                    if all(s == segs[0] for s in segs):
                        common.append(segs[0])
                    else:
                        break
                prefix_len = len(common)
                # extract each file
                for m in files:
                    parts = m.filename.split('/')
                    rel = parts[prefix_len:]
                    if not rel:
                        continue
                    target = models_path.joinpath(*rel)
                    target.parent.mkdir(parents=True, exist_ok=True)
                    with zf.open(m) as src, open(target, 'wb') as dst:
                        dst.write(src.read())
        except Exception as e:
            prog.close()
            ID_utils.displayError(f"❌ Extraction failed:\n{e}")
            self.downloadButton.setDisabled(False)
            return

        # 3) Now verify hashes with live progress updates
        prog.setLabelText("Verifying integrity…")
        # load your reference list
        ref_file = resource_path("idv2_hashes.txt")
        ref_hashes = self._load_reference_hashes(ref_file) or {}
        total = len(ref_hashes)
        prog.setRange(0, total)
        prog.setValue(0)

        ok = True
        checked = 0
        for rel_path, expected_hash in ref_hashes.items():
            # allow cancellation
            QtWidgets.QApplication.processEvents()
            if prog.wasCanceled():
                ok = False
                break

            full_path = models_path / rel_path
            if not full_path.exists():
                ok = False
                break

            actual_hash = self._compute_file_hash(str(full_path))
            if actual_hash != expected_hash:
                ok = False
                break

            checked += 1
            prog.setValue(checked)

        # 4) Wrap up
        prog.close()
        if not ok:
            ID_utils.displayError("❌ Post‑download integrity check failed.")
        else:
            self.show_toast([3000, f"✅ Models ready at:{models_path}"])

        self.downloadButton.setDisabled(False)

    def findligs(self, path):
        out_ligands = []
        unprocessed_ligands = []
        ligands_info = []  # List of tuples (smiles, filename)
        obabel_exe = "obabel"
        self.process = QProcess()

        # Gather all pdbqt, mol2, sdf, and smi files
        for ext in ("*.pdbqt", "*.mol2", "*.sdf", "*.smi"):
            for file in glob(os.path.join(path, ext)):
                if file.endswith("out.pdbqt"):
                    out_ligands.append(file)
                unprocessed_ligands.append(file)

        # If no ligand file is found, return None
        if not unprocessed_ligands:
            return None

        # Loop over each ligand file
        for ligand in unprocessed_ligands:
            ext = os.path.splitext(ligand)[1].lower()
            filename = os.path.basename(ligand)  # Get the filename

            try:
                if ext == ".pdbqt":
                    # Process pdbqt to get SMILES
                    with tempfile.NamedTemporaryFile(suffix=".mol2", delete=False) as tmp:
                        tmp_mol2 = tmp.name
                    args1 = ["-ipdbqt", ligand, "-omol2", "-O", tmp_mol2]
                    self.process.start(obabel_exe, args1)
                    self.process.waitForFinished(-1)

                    smi_input, input_flag = tmp_mol2, "-imol2"
                elif ext == ".mol2":
                    smi_input, input_flag = ligand, "-imol2"
                elif ext == ".sdf":
                    smi_input, input_flag = ligand, "-isdf"
                elif ext == ".smi":
                    # Read SMILES from .smi file
                    with open(ligand, 'r') as f:
                        content = f.read()
                    entries = [s.strip() for s in re.split(r"[\n,]", content) if s.strip()]
                    for entry in entries:
                        ligands_info.append((entry, filename))
                    continue  # Skip obabel conversion
                else:
                    continue

                # Generate SMILES using obabel
                args2 = [input_flag, smi_input, "-osmi"]
                self.process.start(obabel_exe, args2)
                self.process.waitForFinished(-1)

                output = self.process.readAllStandardOutput().data().decode().strip()
                smiles = output.split("\n")[0].split()[0]
                ligands_info.append((smiles, filename))

            finally:
                if ext == ".pdbqt" and 'tmp_mol2' in locals() and os.path.exists(tmp_mol2):
                    os.remove(tmp_mol2)

        return ligands_info
    
    def startPred(self):
        # 1) ensure a ligand folder is selected
        if not self.browseWDLabel.text():
            ID_utils.displayError("❌ Please select a ligand directory first.")
            return

        # 2) quick toast
        self.show_toast([3000, "Starting the prediction process..."])

        # 3) verify the models folder & its hashes
        models_dir = self.get_models_path()
        ref_hash_file = resource_path("idv2_hashes.txt")
        ref_hashes = self._load_reference_hashes(ref_hash_file)

        if (not models_dir.exists()
            or not any(models_dir.iterdir())
            or ref_hashes is None
            or not self._verify_folder_integrity(models_dir, ref_hashes)
        ):
            ID_utils.displayError("❌ Models missing or corrupted.\nPlease (re)download the AI models first.")
            return

        # 4) find ligands
        ligands_info  = self.findligs(self.browseWDLabel.text())
        if ligands_info  is None:
            ID_utils.displayError("❌ No ligands found in the selected directory.")
            return

        smiles_list = [smi for smi, fn in ligands_info]
        filenames_list = [fn for smi, fn in ligands_info]
        self.filenames_list = filenames_list

        # 5) kick off the ADMET prediction
        self.predWindow = predictADMET()
        self.predWindow.smiles(smiles_list)
        self.predWindow.progressBarSignal.connect(self.update_progress_bar)
        self.predWindow.finishedSignal.connect(self.continueADMET)
        self.predWindow.toastSignal.connect(self.show_toast)
        self.startButton.setText("Loading...")
        self.startButton.setDisabled(True)
        self.predWindow.start()
    
    def update_progress_bar(self, val):
        self.progress_bar.setValue(val)

    def generate_rdkit_image_base64(self, smiles, size=(500,500)):
        """Generate base64 image or return error placeholder"""
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            AllChem.Compute2DCoords(mol)
            img = Draw.MolToImage(mol, size=size, fitImage=True)
            buffer = io.BytesIO()
            img.save(buffer, format="PNG")
            encoded = base64.b64encode(buffer.getvalue()).decode("utf-8")
            return f"data:image/png;base64,{encoded}"
        else:
            # Return error placeholder directly
            return get_resource_data_url(':/ErrPlaceholder.png')
        
    def random_pastel_color(self):
        """Generate a random very light pastel color in hex."""
        import colorsys
        hue = random.random()  # 0..1
        saturation = 0.3
        lightness = 0.97
        (r, g, b) = colorsys.hls_to_rgb(hue, lightness, saturation)
        return f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}"

    def get_general_properties(self, prediction_dict, i):
        hba = prediction_dict["HBA"][i]
        hbd = prediction_dict["HBD"][i]
        mp = round(prediction_dict["MP"][i], 4)
        bp = round(prediction_dict["BP"][i], 4)
        return {
            "Molecular Weight": round(prediction_dict["MW"][i], 4),
            "TPSA": round(prediction_dict["TPSA"][i], 4),
            "HBD / HBA": f"{hbd} / {hba}",
            "Rotatable Bonds": prediction_dict["Rotatable Bonds"][i],
            "Nitro Groups": prediction_dict["Nitro groups"][i],
            "Halogen Groups": prediction_dict["Halogen groups"][i],
            "Aniline Groups": prediction_dict["Aniline groups"][i],
            "Benzene rings": prediction_dict["Benzene rings"][i],
            "FractionCSP3": round(prediction_dict["FractionCSP3"][i], 4),
            "Melting Point (℃)": f"{mp}",
            "Boiling Point (℃)": f"{bp}",
        }

    def get_rules(self, prediction_dict, i):
        return {
            "Lipinski Rule": prediction_dict["Lipinski Rule"][i],
            "Ghose rule": prediction_dict["Ghose rule"][i],
            "GSK Rule": prediction_dict["GSK Rule"][i],
            "Golden Triangle": prediction_dict["Golden Triangle"][i],
            "Egan's rule": prediction_dict["Egan's rule"][i],
            "Pfizer Rule": prediction_dict["Pfizer Rule"][i],
            "Veber": prediction_dict["Veber"][i]
        }

    def process_admet_data(self, prediction_dict, i):
        caco2 = round(prediction_dict['CACO2_REG'][i], 4)
        fu = round(prediction_dict['FU'][i], 4)
        ppb = round(prediction_dict['PPB'][i], 4)
        vd = round(prediction_dict['VD'][i], 4)
        cl = round(prediction_dict['CL'][i], 4)
        
        return [
            # Absorption
            {"category": "Absorption", "name": "Oral Bioavailability 50%", "value": "Positive" if prediction_dict['OB'][i] == 1 else "Negative"},
            {"category": "Absorption", "name": "Drug efflux pump (PGP) Inhibitor", "value": "Positive" if prediction_dict['PGP_INHIBITOR'][i] == 1 else "Negative"},
            {"category": "Absorption", "name": "Drug efflux pump (PGP) Substrate", "value": "Positive" if prediction_dict['PGP_SUBSTRATE'][i] == 1 else "Negative"},
            {"category": "Absorption", "name": "Intestinal absorption (HIA)", "value": "Positive" if prediction_dict['HIA'][i] == 1 else "Negative"},
            {"category": "Absorption", "name": "Human colorectal adenocarcinoma cells permeability (Caco2) logPapp", "value": f"{caco2}"},
            # Distribution
            {"category": "Distribution", "name": "Blood-Brain Barrier Permeability (BBB)", "value": "Positive" if prediction_dict['BBB'][i] == 1 else "Negative"},
            {"category": "Distribution", "name": "Fraction unbound in human plasma (%)", "value": f"{fu}"},
            {"category": "Distribution", "name": "Plasma protein binding (%)", "value": f"{ppb}"},
            {"category": "Distribution", "name": "Steady State Volume of Distribution (L/kg)", "value":  f"{vd}"},
            {"category": "Distribution", "name": "LogD (octanol-water dist. coeff.) at pH=7.4", "value": round(prediction_dict['LOGD'][i], 4)},
            {"category": "Distribution", "name": "LogP (octanol-water part. coeff.)", "value": round(prediction_dict['LOGP'][i], 4)},
            {"category": "Distribution", "name": "LogS (water solubility)", "value": round(prediction_dict['LOGS'][i], 4)},
            {"category": "Distribution", "name": "pKa Acid", "value": round(prediction_dict['PKA'][i], 4)},
            {"category": "Distribution", "name": "pKa Basic", "value": round(prediction_dict['PKB'][i], 4)},
            {"category": "Distribution", "name": "Breast cancer resistance protein Inhibitor", "value": "Positive" if prediction_dict['BCRP'][i] == 1 else "Negative"},
            # Metabolism
            {"category": "Metabolism", "name": "CYP 2D6 Inhibitor", "value": "Positive" if prediction_dict['CYP2D6_INHIBITOR'][i] == 1 else "Negative"},
            {"category": "Metabolism", "name": "CYP 2D6 Substrate", "value": "Positive" if prediction_dict['CYP2D6_SUBSTRATE'][i] == 1 else "Negative"},
            {"category": "Metabolism", "name": "CYP 3A4 Inhibitor", "value": "Positive" if prediction_dict['CYP3A4_INHIBITOR'][i] == 1 else "Negative"},
            {"category": "Metabolism", "name": "CYP 3A4 Substrate", "value": "Positive" if prediction_dict['CYP3A4_SUBSTRATE'][i] == 1 else "Negative"},
            {"category": "Metabolism", "name": "Hepatic uptake rate (OATP1B1 Inhibitor)", "value": "Positive" if prediction_dict['OATP1B1'][i] == 1 else "Negative"},
            {"category": "Metabolism", "name": "Hepatic uptake rate (OATP1B3 Inhibitor)", "value": "Positive" if prediction_dict['OATP1B3'][i] == 1 else "Negative"},
            # Excretion
            {"category": "Excretion", "name": "Organic Cation Transporter-2 (OCT2) Inhibitor", "value": "Positive" if prediction_dict['OCT2'][i] == 1 else "Negative"},
            {"category": "Excretion", "name": "Half-life of a drug ≥ 3 hours", "value": "Positive" if prediction_dict['T0.5'][i] == 1 else "Negative"},
            {"category": "Excretion", "name": "Clearance (ml/min/kg)", "value": f"{cl}"},
            # Toxicity
            {"category": "Toxicity", "name": "AMES Mutagenesis", "value": "Positive" if prediction_dict['AMES'][i] == 1 else "Negative"},
            {"category": "Toxicity", "name": "Carcinogenesis", "value": "Positive" if prediction_dict['CARCINOGENICITY'][i] == 1 else "Negative"},
            {"category": "Toxicity", "name": "Drug Induced Liver Injury", "value": "Positive" if prediction_dict['DILI'][i] == 1 else "Negative"},
            {"category": "Toxicity", "name": "Micronucleus Genotoxicity", "value": "Positive" if prediction_dict['MICRONUCLEUS_TOX'][i] == 1 else "Negative"}
        ]

    def generate_properties(self, smiles, prediction_dict):
        testList = []
        smi_counter = 1  # Counter for .smi entries

        for i, smile in enumerate(smiles):
            # Get filename from stored list
            if hasattr(self, 'filenames_list') and i < len(self.filenames_list):
                filename = self.filenames_list[i]
            else:
                filename = smile  # Fallback to smile

            # Determine display name
            if filename.lower().endswith('.smi'):
                display_name = f"Compound {smi_counter}"
                smi_counter += 1
            else:
                # Process filename: strip .pdbqt or _out.pdbqt, then base name + formatting
                fname_lower = filename.lower()
                if fname_lower.endswith('_out.pdbqt'):
                    base = filename[:-len('_out.pdbqt')]
                else:
                    base = os.path.splitext(filename)[0]

                # Remove special characters (non-alphanumeric)
                cleaned = re.sub(r'[^a-zA-Z0-9]', '', base)
                if cleaned:
                    # Capitalize first letter, rest as is
                    formatted = cleaned[0].upper() + cleaned[1:]
                else:
                    # Fallback if cleaned is empty
                    formatted = "Compound"
                display_name = formatted

            # Check if SMILES is valid first
            mol = Chem.MolFromSmiles(smile)
            if mol is None:
                # Create minimal entry with warning
                entry = {
                    "fileName": display_name,
                    "structure_url": self.generate_rdkit_image_base64(smile),
                    "smiles": smile,
                    "warning": "Warning: Invalid SMILES string, IDv2 failed to predict the ADMET properties",
                    "bgcolor": self.random_pastel_color(),
                    "radar_labels": ["LIPO", "FLEX", "SIZE", "POLAR", "INSOLU", "INSATU"],
                    "radar_values": [0, 0, 0, 0, 0, 0],
                    "radar_original_values": [0, 0, 0, 0, 0, 0],
                    "generalProps": {},
                    "rules": {},
                    "admetProps": []
                }
                testList.append(entry)
                continue

            # Proceed with valid SMILES
            original_values = [
                prediction_dict["LOGP"][i],
                prediction_dict["Rotatable Bonds"][i],
                prediction_dict["MW"][i],
                prediction_dict["TPSA"][i],
                prediction_dict["LOGS"][i],
                prediction_dict["FractionCSP3"][i]
            ]
            
            normalized_values = [
                (original_values[0] + 2) / 8,
                original_values[1] / 10,
                (original_values[2] - 100) / 400,
                original_values[3] / 150,
                (original_values[4] + 6) / 6.5,
                original_values[5]
            ]
            
            entry = {
                "fileName": display_name,
                "structure_url": self.generate_rdkit_image_base64(smile) or "N/A",
                "smiles": smile,
                "radar_labels": ["LIPO", "FLEX", "SIZE", "POLAR", "INSOLU", "INSATU"],
                "radar_values": normalized_values,
                "radar_original_values": original_values,
                "generalProps": self.get_general_properties(prediction_dict, i),
                "rules": self.get_rules(prediction_dict, i),
                "admetProps": self.process_admet_data(prediction_dict, i),
                "bgcolor": self.random_pastel_color()
            }
            testList.append(entry)

        return testList

    
    def continueADMET(self, admet_data):
        self.startButton.setText("Predict ADMET")
        self.startButton.setDisabled(False)
        error_placeholder_url = get_resource_data_url(':/ErrPlaceholder.png')
        smiles = admet_data['SMILES']
        compound_properties = self.generate_properties(smiles, admet_data)
        template_dir = os.path.dirname(resource_path('template.html'))
        env = Environment(loader=FileSystemLoader(template_dir))
        template = env.get_template('template.html')
        chart_js_url = "chart.js"
        output_html = template.render(data=compound_properties, error_placeholder_url=error_placeholder_url, chart_js_url=chart_js_url)
        with open("test.html", "w", encoding="utf-8") as f:
            f.write(output_html)
        
        # Create WebEngineWindow and load the saved HTML file
        self.webWindow = WebEngineWindow(self)
        self.webWindow.view.setUrl(QUrl.fromLocalFile(os.path.abspath("test.html")))
        
        # Connect download handler AFTER page is created
        if hasattr(self, 'on_downloadRequested'):
            page = self.webWindow.view.page()
            if page is not None:
                page.profile().downloadRequested.connect(self.on_downloadRequested)
        
        self.webWindow.show()

    def on_downloadRequested(self, download_item):
        """
        Called whenever the page initiates a download (e.g. user clicks a link
        to a .csv, .pdf, or JS code triggers a download).
        """
        # Suggest a filename based on the URL:
        default_path = os.path.join(
            os.getcwd(),
            os.path.basename(download_item.url().path())
        )
        # Ask the user where to save:
        save_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save CSV As",
            default_path,
            "CSV Files (*.csv)"
        )
        if save_path:
            download_item.setPath(save_path)
            download_item.accept()    # actually start the download
        else:
            download_item.cancel()    # user backed out

    def show_toast(self, text):
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

if hasattr(QtCore.Qt, 'AA_EnableHighDpiScaling'):
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
if hasattr(QtCore.Qt, 'AA_UseHighDpiPixmaps'):
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)