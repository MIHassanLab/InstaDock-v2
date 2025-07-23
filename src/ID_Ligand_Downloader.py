from PyQt5.QtCore import (
    QSize, Qt, QSize)
from PyQt5.QtGui import QIcon, QPixmap, QFont
from PyQt5.QtWidgets import (
    QPushButton, QTabWidget, QGridLayout, QHBoxLayout, QLabel, QMainWindow, QToolButton,
    QVBoxLayout, QWidget, QSpacerItem, QSizePolicy, QGroupBox, QProgressBar,
    QFrame, QDesktopWidget, QCheckBox, QPlainTextEdit, QFileDialog
)
from PyQt5 import QtCore, QtWidgets
import ID_utils
from pyqttoast import Toast, ToastPreset
from QSwitchControl import SwitchControl
import IDv2Resource

class CustomTitleBar(QWidget):
    def __init__(self, parent):
        super().__init__(parent)
        self.initial_pos = None
        title_bar_layout = QHBoxLayout(self)
        title_bar_layout.setContentsMargins(1, 1, 1, 1)
        title_bar_layout.setSpacing(2)

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


class DownloadUI(QMainWindow):
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
        self.setWindowIcon(QIcon(':/Logo.png'))
        self.setWindowTitle('''InstaDock-v2 (Downloader)''')
        self.resize(500, 500)
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
                            
        QComboBox {
            font-size: 14px; /* Matches the font size */
            padding: 4px 8px; /* Inner spacing for better usability */
            border: 1px solid #A9A9A9; /* Light border for light mode */
            border-radius: 5px; /* Rounded corners for modern look */
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, 
            stop:0 #F2F2F2, stop:0.5 #E6E6E6, stop:1 #D9D9D9); /* Light gradient background */
            color: black; /* Text color for readability */
        }
                            
        QLineEdit {
            font-size: 14px; /* Matches the font size */
            padding: 1px 1px; /* Inner spacing */
        }

        QTabWidget::pane {
            border: 1px solid #A9A9A9; /* Matching border color to QPushButton */
            background: #F2F2F2; /* Soft light background to keep it subtle */
        }

        QTabBar::tab {
            width: 157px;
            height: 35px;
            font-size: 14px;
            color: black; /* Matching text color */
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, 
                        stop:0 #F2F2F2, /* Start gradient with light grey */
                        stop:0.5 #EBEBEB, /* Mid-point for light grey */
                        stop:1 #DDDDDD); /* End gradient with a slightly darker grey */
            border: 1px solid #A9A9A9; /* Border to match QPushButton */
            border-bottom: 0px;
        }

        QTabBar::tab:hover {
            width: 157px;
            height: 35px;
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, 
                        stop:0 #E0E0E0, /* Slightly darker hover color */
                        stop:0.5 #D4D4D4, /* Smooth transition */
                        stop:1 #C8C8C8); /* End hover gradient */
            border: 1px solid #37cbe6; /* Highlighted border on hover */
            border-bottom: 0px;
        }

        QTabBar::tab:selected {
            width: 157px;
            height: 35px;
            background: white; /* Keeping the selected tab background white */
            border: 1px solid #37cbe6; /* Slight emphasis on the selected tab with matching border */
            border-bottom: 0px;
            color: black; /* Dark text for readability */
        }

        QCheckBox {
            font-size: 14px;
        }

        QCheckBox::indicator {
            width: 20px;
            height: 20px;
            border: 1px solid #A9A9A9;
            border-radius: 4px;
            }

        QCheckBox::indicator:checked {
            border: 5px solid  #12c6c2;
            width: 12px;
            height: 12px;
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
            width: 10px; 
            margin: 0.5px; 
            border-radius: 1.5px; 
        }
        
        """)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sizePolicy)
        self.setMinimumSize(QSize(500, 650))
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

        sizePolicy1 = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Maximum)
        sizePolicy1.setHorizontalStretch(0)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy2 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        sizePolicy2.setHorizontalStretch(0)
        sizePolicy2.setVerticalStretch(0)
        sizePolicy3 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        sizePolicy3.setHorizontalStretch(0)
        sizePolicy3.setVerticalStretch(0)

        self.introLabel = QLabel(self)
        self.introLabel.setObjectName("introLabel")
        self.introLabel.setTextFormat(Qt.RichText)
        self.introLabel.setWordWrap(True)
        self.introLabel.setText('''
        Welcome to the Ligand Download Module!<br/><br/>
        This tool enables you to download multiple ligand structures directly from major chemical databases: PubChem, ZINC, and IMPPAT. Simply enter the compound IDs or names into the respective fields, and the module will fetch the corresponding 3D structure files (and other formats where available) with real-time progress updates.<br/>
        Enjoy a streamlined and efficient workflow for building your ligand library!
        ''')

        self.introLabel.setMaximumSize(QSize(480,150))
        sizePolicy1.setHeightForWidth(self.introLabel.sizePolicy().hasHeightForWidth())
        self.introLabel.setSizePolicy(sizePolicy1)
        self.introLabel.setAlignment(Qt.AlignJustify)

        work_space_layout.addWidget(self.introLabel, Qt.AlignCenter)

        self.groupWD = QGroupBox(self)
        self.groupWD.setObjectName("groupWD")
        sizePolicy2.setHeightForWidth(self.groupWD.sizePolicy().hasHeightForWidth())
        self.groupWD.setSizePolicy(sizePolicy2)
        self.groupWD.setMinimumSize(QSize(0, 50))

        self.gridLayout1 = QGridLayout(self.groupWD)
        self.gridLayout1.setObjectName("gridLayout1")
        self.gridLayout1.setContentsMargins(0, 0, 0, 0)

        self.browseWDLabel = QLabel(self)
        self.browseWDLabel.setObjectName("browseWDLabel")
        self.browseWDLabel.setText("Browse folder to download the files...")
        self.browseWDLabel.setDisabled(True)
        self.browseWDLabel.setStyleSheet('''QLabel {
        padding: 1px 1px; /* Inner spacing */
        }''')
        sizePolicy2.setHeightForWidth(self.browseWDLabel.sizePolicy().hasHeightForWidth())
        self.browseWDLabel.setSizePolicy(sizePolicy2)
        self.browseWDLabel.setMinimumSize(QSize(400, 0))
        self.browseWDLabel.setEnabled(False)
        self.browseWDLabel.setFrameShape(QFrame.Box)
        self.gridLayout1.addWidget(self.browseWDLabel, 0, 0, 1, 2, Qt.AlignLeft)

        self.browseWDButton = QPushButton("  Browse", self)
        self.browseWDButton.setObjectName("browseWDButton")
        browseIcon = QPixmap(":/Browse.png")
        self.browseWDButton.setIcon(QIcon(browseIcon))
        self.browseWDButton.setIconSize(QSize(30, 30)) 
        sizePolicy2.setHeightForWidth(self.browseWDButton.sizePolicy().hasHeightForWidth())
        self.browseWDButton.setSizePolicy(sizePolicy2)
        self.browseWDButton.clicked.connect(self.browse)
        self.gridLayout1.addWidget(self.browseWDButton, 0, 1, 1, 1, Qt.AlignRight)

        work_space_layout.addWidget(self.groupWD)

        self.tabWidget = QTabWidget(self)
        self.tabWidget.setObjectName("tabWidget")
        sizePolicy2.setHeightForWidth(self.tabWidget.sizePolicy().hasHeightForWidth())
        self.tabWidget.setSizePolicy(sizePolicy2)
        self.tabWidget.setMinimumSize(QSize(0, 300))
        self.tabWidget.setUsesScrollButtons(False)

        self.pubchem = QWidget()
        self.pubchem.setObjectName("pubchem")
        self.verticalLayout = QVBoxLayout(self.pubchem)
        self.verticalLayout.setObjectName("verticalLayout")

        self.pubchemOptions = QGroupBox(self.pubchem)
        self.pubchemOptions.setObjectName("pubchemOptions")
        self.gridLayout2 = QGridLayout(self.pubchemOptions)
        self.gridLayout2.setObjectName("gridLayout2")
        
        self.tg1Group = QGroupBox(self.pubchemOptions)
        self.tg1Group.setObjectName("tg1Group")
        self.horizontalLayout1 = QHBoxLayout(self.tg1Group)
        self.horizontalLayout1.setObjectName("horizontalLayout1")

        self.cid = QLabel(self.tg1Group)
        self.cid.setObjectName("cid")
        self.cid.setText("Compound ID")
        sizePolicy2.setHeightForWidth(self.cid.sizePolicy().hasHeightForWidth())
        self.cid.setSizePolicy(sizePolicy2)
        self.horizontalLayout1.addWidget(self.cid)

        self.toggle_1 = SwitchControl(bg_color="#777777", circle_color="#DDD", active_color="#12c6c2", animation_curve=QtCore.QEasingCurve.InOutCubic, animation_duration=100, checked=True, change_cursor=True)
        self.toggle_1.stateChanged.connect(self.on_switch1_toggled)
        self.horizontalLayout1.addWidget(self.toggle_1)

        self.tg2Group = QGroupBox(self.pubchemOptions)
        self.tg2Group.setObjectName("tg1Group")
        self.horizontalLayout2 = QHBoxLayout(self.tg2Group)
        self.horizontalLayout2.setObjectName("horizontalLayout2")

        self.names = QLabel(self.tg2Group)
        self.names.setObjectName("names")
        self.names.setText("Names")
        sizePolicy2.setHeightForWidth(self.names.sizePolicy().hasHeightForWidth())
        self.names.setSizePolicy(sizePolicy2)
        self.horizontalLayout2.addWidget(self.names)

        self.toggle_2 = SwitchControl(bg_color="#777777", circle_color="#DDD", active_color="#12c6c2", animation_curve=QtCore.QEasingCurve.InOutCubic, animation_duration=100, checked=False, change_cursor=True)
        self.toggle_2.stateChanged.connect(self.on_switch2_toggled)
        self.horizontalLayout2.addWidget(self.toggle_2)

        self.gridLayout2.addWidget(self.tg1Group, 0, 0, 1, 1, Qt.AlignCenter)
        self.gridLayout2.addWidget(self.tg2Group, 0, 1, 1, 1, Qt.AlignCenter)
        
        self.verticalLayout.addWidget(self.pubchemOptions)

        self.textBrowserPC = QPlainTextEdit(self)
        self.textBrowserPC.setObjectName("textBrowserPC")
        self.textBrowserPC.setReadOnly(False)
        self.textBrowserPC.setPlaceholderText(
        'Enter one or more PubChem compound IDs or names (e.g., "Aspirin"), separated by newlines.'
        )
        self.textBrowserPC.setStyleSheet("QTextBrowser { font-size: 14px; }")
        self.textBrowserPC.setEnabled(True)

        self.verticalLayout.addWidget(self.textBrowserPC)

        self.tabWidget.addTab(self.pubchem, "PubChem")

        self.zinc = QWidget()
        self.zinc.setObjectName("zinc")
        self.verticalLayout = QVBoxLayout(self.zinc)
        self.verticalLayout.setObjectName("verticalLayout")

        self.textBrowserZN = QPlainTextEdit(self)
        self.textBrowserZN.setObjectName("textBrowserZN")
        self.textBrowserZN.setReadOnly(False)
        self.textBrowserZN.setPlaceholderText(
        'Enter one or more ZINC IDs, separated by commas, tabs, or newlines.'
        )
        self.textBrowserZN.setStyleSheet("QTextBrowser { font-size: 14px; }")
        self.textBrowserZN.setEnabled(True)

        self.verticalLayout.addWidget(self.textBrowserZN)

        self.tabWidget.addTab(self.zinc, "ZINC")

        self.imppat = QWidget()
        self.imppat.setObjectName("imppat")
        self.verticalLayout = QVBoxLayout(self.imppat)
        self.verticalLayout.setObjectName("verticalLayout")

        self.imppatOptions = QGroupBox(self.imppat)
        self.imppatOptions.setObjectName("imppatOptions")
        self.gridLayout2 = QGridLayout(self.imppatOptions)
        self.gridLayout2.setObjectName("gridLayout2")

        self.sdf = QCheckBox(self.imppatOptions)
        self.sdf.setObjectName("sdf")
        self.sdf.setText("SDF")
        self.sdf.setAutoExclusive(True)
        self.sdf.setChecked(True)
        sizePolicy1.setHeightForWidth(self.sdf.sizePolicy().hasHeightForWidth())
        self.sdf.setSizePolicy(sizePolicy1)
        self.gridLayout2.addWidget(self.sdf, 0, 0, 1, 1)

        self.mol2 = QCheckBox(self.imppatOptions)
        self.mol2.setObjectName("mol2")
        self.mol2.setText("MOL2")
        self.mol2.setAutoExclusive(True)
        sizePolicy1.setHeightForWidth(self.mol2.sizePolicy().hasHeightForWidth())
        self.mol2.setSizePolicy(sizePolicy1)
        self.gridLayout2.addWidget(self.mol2, 0, 1, 1, 1)

        self.pdb = QCheckBox(self.imppatOptions)
        self.pdb.setObjectName("pdb")
        self.pdb.setText("PDB")
        self.pdb.setAutoExclusive(True)
        sizePolicy1.setHeightForWidth(self.pdb.sizePolicy().hasHeightForWidth())
        self.pdb.setSizePolicy(sizePolicy1)
        self.gridLayout2.addWidget(self.pdb, 0, 2, 1, 1)

        self.verticalLayout.addWidget(self.imppatOptions)

        self.textBrowserIM = QPlainTextEdit(self)
        self.textBrowserIM.setObjectName("textBrowserIM")
        self.textBrowserIM.setReadOnly(False)
        self.textBrowserIM.setPlaceholderText(
        'Enter one or more IMPPAT IDs, separated by commas, tabs, or newlines.'
        )
        self.textBrowserIM.setStyleSheet("QTextBrowser { font-size: 14px; }")
        self.textBrowserIM.setEnabled(True)

        self.verticalLayout.addWidget(self.textBrowserIM)

        self.tabWidget.addTab(self.imppat, "IMPPAT")

        work_space_layout.addWidget(self.tabWidget, Qt.AlignCenter)

        self.downloadButton = QPushButton(self)
        self.downloadButton.setObjectName("downloadButton")
        self.downloadButton.setText("Download")
        self.downloadButton.clicked.connect(self.testDownload)
        self.downloadButton.setEnabled(True)

        work_space_layout.addWidget(self.downloadButton, Qt.AlignCenter)

        self.progressBar = QProgressBar(self)
        self.progressBar.setObjectName("progressBar")
        self.progressBar.setTextVisible(False)
        sizePolicy3.setHeightForWidth(self.progressBar.sizePolicy().hasHeightForWidth())
        self.progressBar.setSizePolicy(sizePolicy3)
        self.progressBar.setValue(0)
        work_space_layout.addWidget(self.progressBar, Qt.AlignCenter)

        central_widget_layout.addLayout(work_space_layout)
        central_widget.setLayout(central_widget_layout)
        self.setCentralWidget(central_widget)

    def browse(self):
        self.browseWDLabel.setText(QFileDialog.getExistingDirectory(None, 'Select working directory'))

    def testDownload(self):
        currTab = self.tabWidget.currentWidget().objectName()
        if self.browseWDLabel.text() == "Browse folder to download the files...":
            self.notifyUser("❌ Please select a working directory.")
        else:    
            if currTab == "pubchem":
                self.pubChemEngine = ID_utils.pubchemEngine()
                self.pubChemEngine.setWD(self.browseWDLabel.text())
                self.pubChemEngine.setFormat([self.toggle_1.isChecked(), self.toggle_2.isChecked()])
                if self.toggle_1.isChecked():
                    plain_text = self.textBrowserPC.toPlainText()
                    self.pubChemEngine.setIDs(plain_text)
                else:
                    plain_text = self.textBrowserPC.toPlainText()
                    self.pubChemEngine.setIDs(plain_text)
                self.pubChemEngine.finishedSignal.connect(self.notifyUser)
                self.pubChemEngine.progressBarSignalTotal.connect(self.update_total_progress_bar)
                self.pubChemEngine.progressBarSignal.connect(self.update_progress_bar)
                self.pubChemEngine.toastSignal.connect(self.show_toast)
                self.pubChemEngine.start()
            elif currTab == "zinc":
                self.zincEngine = ID_utils.zincEngine()
                self.zincEngine.setWD(self.browseWDLabel.text())
                self.zincEngine.setIDs(self.textBrowserZN.toPlainText())
                self.zincEngine.finishedSignal.connect(self.notifyUser)
                self.zincEngine.progressBarSignalTotal.connect(self.update_total_progress_bar)
                self.zincEngine.progressBarSignal.connect(self.update_progress_bar)
                self.zincEngine.toastSignal.connect(self.show_toast)
                self.zincEngine.start()
            elif currTab == "imppat":
                self.imppatEngine = ID_utils.imppatEngine()
                self.imppatEngine.setWD(self.browseWDLabel.text())
                self.imppatEngine.setIDs(self.textBrowserIM.toPlainText())
                self.imppatEngine.setFormat([
                    self.sdf.isChecked(), self.mol2.isChecked(), self.pdb.isChecked()
                ])
                self.imppatEngine.progressBarSignalTotal.connect(self.update_total_progress_bar)
                self.imppatEngine.progressBarSignal.connect(self.update_progress_bar)
                self.imppatEngine.finishedSignal.connect(self.notifyUser)
                self.imppatEngine.toastSignal.connect(self.show_toast)
                self.imppatEngine.start()

    def update_total_progress_bar(self, value):
        self.progressBar.setMaximum(value)

    def update_progress_bar(self, value):
        self.progressBar.setValue(value)

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

    def notifyUser(self, message):
        if message.startswith("❌"):
            ID_utils.displayError(message)
        elif message.startswith("✅"):
            ID_utils.displayMessage(message)

    def show_toast(self, text):
        toastfont = QFont('Arial', 12, QFont.Weight.Medium)
        toast = Toast()
        toast.setDuration(4000)
        toast.setText(text)
        toast.setTextFont(toastfont)
        toast.applyPreset(ToastPreset.INFORMATION)
        toast.setResetDurationOnHover(False)
        toast.setBorderRadius(3) 
        toast.setIcon(QPixmap(':/Logo.png'))
        toast.setShowCloseButton(False) 
        toast.setIconSize(QSize(40, 40)) 
        toast.setIconColor(None)
        toast.show()

if hasattr(QtCore.Qt, 'AA_EnableHighDpiScaling'):
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
if hasattr(QtCore.Qt, 'AA_UseHighDpiPixmaps'):
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)