import warnings
warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=DeprecationWarning)

import shutil
import ID_utils
import subprocess
import IDv2Resource
from ID_ADMET import admet
from ID_Lipinski import ro5
from ID_2D import idTwoDInteract
from ID_Visualizer import idVisualizer
from pyqttoast import Toast, ToastPreset
from ID_About_Multi import AboutMultiUI
from PyQt5.QtCore import QSize, Qt, QDir
from PyQt5.QtGui import QIcon, QPalette, QPixmap, QFontDatabase, QFont
from ID_Start_Calc_Multi import (multiCalculateS1, multiCalculateS2, multiCalculateS3, multiCalculateS4, startMultiTargetDocking, continueMultiID)
from ID_Prep_Multi import (prepareTargetM1, prepareTargetM2, prepareTargetM3, prepareLigand, prepareConfigM1, prepareConfigM2, prepareConfigM3, convertFiles, inspectProtein)
from PyQt5.QtWidgets import (
    QPushButton, QTabWidget, QGridLayout, QHBoxLayout, QLabel, QMainWindow, QToolButton,
    QVBoxLayout, QWidget, QSpacerItem, QSizePolicy, QGroupBox, QComboBox, QProgressBar,
    QFrame, QDesktopWidget, QLineEdit, QApplication, QStyle, QStylePainter, QStyleOptionComboBox, 
    QFileDialog, QFormLayout, QInputDialog
)
    
def gpu_available():
    if shutil.which("nvidia-smi") is None:
        return False
    try:
        out = subprocess.check_output(["nvidia-smi", "-L"], stderr=subprocess.DEVNULL)
        return bool(out.strip())
    except subprocess.CalledProcessError:
        return False

class ComboBox(QComboBox):
    def paintEvent(self, e):
        painter = QStylePainter(self)
        painter.setPen(self.palette().color(QPalette.Text))

        # draw the combobox frame, focusrect and selected etc.
        opt = QStyleOptionComboBox()
        self.initStyleOption(opt)
        painter.drawComplexControl(QStyle.CC_ComboBox, opt)

        if self.currentIndex() < 0:
            opt.palette.setBrush(
                QPalette.ButtonText,
                opt.palette.brush(QPalette.ButtonText).color().lighter(),
            )
            if self.placeholderText():
                opt.currentText = self.placeholderText()

        # draw the icon and text
        painter.drawControl(QStyle.CE_ComboBoxLabel, opt)

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

class MultiTargetUI(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.setupUi(self)
    def center_above_screen(self):
        screen_geometry = QDesktopWidget().screenGeometry()
        window_geometry = self.frameGeometry()
        screen_center = screen_geometry.center()
        window_geometry.moveCenter(screen_center)
        self.move(window_geometry.topLeft())

    def __init__(self):
        super().__init__()
        self.setWindowIcon(QIcon(':/IDLogo.png'))
        self.setWindowTitle("IDv2 (Multi Target)")
        self.resize(650, 800)
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
                           
        QLabel#browseTargetsLabel{
            background: transparent;
            font-size: 14px;
            border: 1px solid #A9A9A9;
            border-radius: 3px;
        }
        
        QLabel#browseLigandsLabel{
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
                           
        QComboBox::item:selected {
            background-color: #D4EDDA !important;
            color: black !important;
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
            width: 154px;
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
            width: 153px;
            height: 35px;
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, 
                        stop:0 #E0E0E0, /* Slightly darker hover color */
                        stop:0.5 #D4D4D4, /* Smooth transition */
                        stop:1 #C8C8C8); /* End hover gradient */
            border: 1px solid #37cbe6; /* Highlighted border on hover */
            border-bottom: 0px;
        }

        QTabBar::tab:selected {
            width: 153px;
            height: 35px;
            background: white; /* Keeping the selected tab background white */
            border: 1px solid #37cbe6; /* Slight emphasis on the selected tab with matching border */
            border-bottom: 0px;
            color: black; /* Dark text for readability */
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
                           
         QInputDialog {
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:1,
                stop:0 #F5F7FA,
                stop:0.5 #E4E7EB,
                stop:1 #D9DEE4);
            color: #000000;
            border: 1px solid #A9A9A9;
        }

        QInputDialog QLabel {
            background: transparent;
            font-size: 14px;
            font-family: "Arial";
        }

        QInputDialog QLineEdit {
            font-size: 14px;
            padding: 1px 1px;
        }
    
    """
)
        sizePolicy = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sizePolicy)
        self.setMinimumSize(QSize(650, 800))
        self.setWindowFlags(Qt.FramelessWindowHint)
        central_widget = QWidget()
        central_widget.setObjectName("Container")
        self.title_bar = CustomTitleBar(self)

        central_widget_layout = QVBoxLayout()
        central_widget_layout.setContentsMargins(0, 0, 0, 0)
        central_widget_layout.setAlignment(Qt.AlignTop)

        central_widget_layout.addWidget(self.title_bar)

        work_space_layout = QGridLayout()
        work_space_layout.setContentsMargins(11, 11, 11, 11)

        sizePolicy1 = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Minimum)
        sizePolicy1.setHorizontalStretch(0)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy2 = QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)
        sizePolicy2.setHorizontalStretch(0)
        sizePolicy2.setVerticalStretch(0)
        sizePolicy3 = QSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.MinimumExpanding)
        sizePolicy3.setHorizontalStretch(0)
        sizePolicy3.setVerticalStretch(0)
        sizePolicy4 = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        sizePolicy4.setHorizontalStretch(0)
        sizePolicy4.setVerticalStretch(0)
        sizePolicy5 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        sizePolicy5.setHorizontalStretch(0)
        sizePolicy5.setVerticalStretch(0)
        sizePolicy6 = QSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.Fixed)
        sizePolicy6.setHorizontalStretch(0)
        sizePolicy6.setVerticalStretch(0)

        self.groupTitle = QGroupBox(self)
        self.groupTitle.setObjectName("groupTitle")
        sizePolicy2.setHeightForWidth(self.groupTitle.sizePolicy().hasHeightForWidth())
        self.groupTitle.setSizePolicy(sizePolicy2)
        self.groupTitle.setMinimumSize(QSize(630, 0))

        self.verticalLayout1 = QVBoxLayout(self.groupTitle)
        self.verticalLayout1.setObjectName("verticalLayout1")
        self.verticalLayout1.setContentsMargins(0, 0, 0, 0)

        self.titleLogo = QLabel(self.groupTitle)
        self.titleLogo.setObjectName("titleLogo")
        sizePolicy4.setHeightForWidth(self.titleLogo.sizePolicy().hasHeightForWidth())
        self.titleLogo.setSizePolicy(sizePolicy4)
        pixmap = QPixmap(":/IDLogoWithText.png")
        scaled_pixmap = pixmap.scaled(
            630, 200,  # Adjust width and height as needed
            Qt.KeepAspectRatio,
            Qt.SmoothTransformation
        )
        self.titleLogo.setPixmap(scaled_pixmap)
        self.titleLogo.setScaledContents(True)
        self.titleLogo.setAlignment(Qt.AlignCenter)
        self.verticalLayout1.addWidget(self.titleLogo)

        self.titleLabel = QLabel(self)
        self.titleLabel.setObjectName("titleLabel")
        self.titleLabel.setTextFormat(Qt.RichText)
        self.titleLabel.setWordWrap(True)
        self.titleLabel.setText('''
<div style="text-align:justify; margin:1.5em 1em 0 1em; line-height:1.2;">
  InstaDock-v2 (IDv2) delivers dependable docking results, integrated ADMET insights, and streamlined batch workflows, enabling you to tackle complex screening campaigns with confidence.
</div>

<div style="text-align:justify; margin:1em 1em 0 1em; line-height:1.2;">
  The <strong>Multi Target</strong> mode enables simultaneous docking across multiple receptors, accelerating large-scale molecular modeling tasks.
</div>
                                ''')
        #self.titleLabel.setStyleSheet("QLabel { color : #FF9D00; font-size: 50pt; }")
        sizePolicy1.setHeightForWidth(self.titleLabel.sizePolicy().hasHeightForWidth())
        self.titleLabel.setSizePolicy(sizePolicy1)
        self.titleLabel.setMinimumSize(QSize(630, 0))
        self.titleLabel.setAlignment(Qt.AlignJustify)
        self.verticalLayout1.addWidget(self.titleLabel)

        work_space_layout.addWidget(self.groupTitle, 0, 0, 1, 2, Qt.AlignCenter)

        self.groupWD = QGroupBox(self)
        self.groupWD.setObjectName("groupWD")
        
        sizePolicy2.setHeightForWidth(self.groupWD.sizePolicy().hasHeightForWidth())
        self.groupWD.setSizePolicy(sizePolicy2)
        self.groupWD.setMinimumSize(QSize(627, 50))

        self.gridLayout = QGridLayout(self.groupWD)
        self.gridLayout.setObjectName("gridLayout")
        self.gridLayout.setContentsMargins(0, 0, 0, 0)

        self.groupDiff = QGroupBox(self)
        self.groupDiff.setObjectName("groupDiff")
        self.gridLayout1 = QGridLayout(self.groupDiff)
        self.gridLayout1.setObjectName("gridLayout1")
        self.gridLayout1.setContentsMargins(0, 0, 0, 0)

        self.browseTargetsLabel = QLabel(self)
        self.browseTargetsLabel.setObjectName("browseTargetsLabel")
        self.browseTargetsLabel.setText("Browse target directory...")
        self.browseTargetsLabel.setEnabled(False)
        self.browseTargetsLabel.setStyleSheet('''QLabel {
        padding: 1px 1px; /* Inner spacing */
        }''')
        self.browseTargetsLabel.setMinimumSize(QSize(400, 0))
        self.browseTargetsLabel.setFrameShape(QFrame.Box)
        self.gridLayout1.addWidget(self.browseTargetsLabel, 0, 0, 1, 2, Qt.AlignLeft)

        self.browseTargetsButton = QPushButton("  Browse", self)
        self.browseTargetsButton.setObjectName("browseTargetsButton")
        browseIcon = QPixmap(":/Browse.png")
        self.browseTargetsButton.setIcon(QIcon(browseIcon))
        self.browseTargetsButton.setIconSize(QSize(20, 20)) 
        sizePolicy3.setHeightForWidth(self.browseTargetsButton.sizePolicy().hasHeightForWidth())
        self.browseTargetsButton.setSizePolicy(sizePolicy3)
        self.browseTargetsButton.clicked.connect(self.browseTargets)
        self.gridLayout1.addWidget(self.browseTargetsButton, 0, 1, 1, 1, Qt.AlignRight)

        self.browseLigandsLabel = QLabel(self)
        self.browseLigandsLabel.setObjectName("browseLigandsLabel")
        self.browseLigandsLabel.setText("Browse ligand directory...")
        self.browseLigandsLabel.setEnabled(False)
        self.browseLigandsLabel.setStyleSheet('''QLabel {
        padding: 1px 1px; /* Inner spacing */
        }''')
        self.browseLigandsLabel.setMinimumSize(QSize(400, 0))
        self.browseLigandsLabel.setFrameShape(QFrame.Box)
        self.gridLayout1.addWidget(self.browseLigandsLabel, 1, 0, 1, 2, Qt.AlignLeft)

        self.browseLigandsButton = QPushButton("  Browse", self)
        self.browseLigandsButton.setObjectName("browseLigandsButton")
        browseIcon = QPixmap(":/Browse.png")
        self.browseLigandsButton.setIcon(QIcon(browseIcon))
        self.browseLigandsButton.setIconSize(QSize(20, 20)) 
        sizePolicy3.setHeightForWidth(self.browseLigandsButton.sizePolicy().hasHeightForWidth())
        self.browseLigandsButton.setSizePolicy(sizePolicy3)
        self.browseLigandsButton.clicked.connect(self.browseLigands)
        self.gridLayout1.addWidget(self.browseLigandsButton, 1, 1, 1, 1, Qt.AlignRight)

        self.gridLayout.addWidget(self.groupDiff, 0, 0, 1, 2)

        self.testGroup = QGroupBox(self)
        self.testGroup.setObjectName("testGroup")
        self.verticalLayout1 = QVBoxLayout(self.testGroup)
        self.verticalLayout1.setObjectName("verticalLayout1")
        self.verticalLayout1.setContentsMargins(0, 0, 0, 0)

        self.comboBox = ComboBox(self.groupWD)
        self.comboBox.setPlaceholderText("Select docking method…")
        self.comboBox.setObjectName("comboBox")

        # build the list dynamically
        dock_list = ["AutoDock Vina", "QuickVina-W"]
        if gpu_available():
            dock_list.append("AutoDock-GPU")

        self.comboBox.addItems(dock_list)
        self.comboBox.setCurrentIndex(1)
        self.comboBox.currentIndexChanged.connect(self.changedIndex)
        sizePolicy3.setHeightForWidth(self.comboBox.sizePolicy().hasHeightForWidth())
        self.comboBox.setSizePolicy(sizePolicy3)
        self.verticalLayout1.addWidget(self.comboBox)

        self.aboutButton = QPushButton("  About the tool", self)
        self.aboutButton.setObjectName("aboutButton")
        aboutIcon = QPixmap(":/About.png")
        self.aboutButton.setIcon(QIcon(aboutIcon))
        self.aboutButton.setIconSize(QSize(30, 30))
        self.aboutButton.setStyleSheet("text-align:left;")
        sizePolicy3.setHeightForWidth(self.aboutButton.sizePolicy().hasHeightForWidth())
        self.aboutButton.setSizePolicy(sizePolicy3)
        self.aboutButton.clicked.connect(self.connectAbout)
        self.verticalLayout1.addWidget(self.aboutButton)

        self.gridLayout.addWidget(self.testGroup, 0, 3, 1, 2)

        work_space_layout.addWidget(self.groupWD, 2, 0, 1, 2)

        self.groupTabs = QGroupBox(self)
        self.groupTabs.setObjectName("groupTabs")
        
        self.horizontalLayout1 = QHBoxLayout(self.groupTabs)
        self.horizontalLayout1.setObjectName("horizontalLayout1")
        self.horizontalLayout1.setContentsMargins(0, 0, 0, 0)

        self.tabWidget1 = QTabWidget(self.groupTabs)
        self.tabWidget1.setObjectName("tabWidget1")
        sizePolicy5.setHeightForWidth(self.tabWidget1.sizePolicy().hasHeightForWidth())
        self.tabWidget1.setSizePolicy(sizePolicy5)
        self.tabWidget1.setUsesScrollButtons(False)
                
        self.dock = QWidget()
        self.dock.setObjectName("dock")
        self.verticalLayout = QVBoxLayout(self.dock)
        self.verticalLayout.setObjectName("verticalLayout")
        
        self.blind = QPushButton("  Blind docking", self.dock)
        self.blind.setObjectName("blind")
        self.blind.setCheckable(True)
        self.blind.setChecked(True)
        self.blind.clicked.connect(self.connectBlind)
        blindIcon = QPixmap(":/Blind.png")
        self.blind.setIcon(QIcon(blindIcon))
        self.blind.setIconSize(QSize(30, 30))
        self.blind.setMinimumSize(QSize(120, 50))
        self.blind.setStyleSheet("text-align:left;")
        self.verticalLayout.addWidget(self.blind)

        self.siteSpecific = QPushButton("  Site specific docking",self.dock)
        self.siteSpecific.setObjectName("siteSpecific")
        self.siteSpecific.setAutoExclusive(True)    
        self.siteSpecific.setCheckable(True)
        self.siteSpecific.clicked.connect(self.connectSiteSpecific)
        siteIcon = QPixmap(":/Site.png")
        self.siteSpecific.setIcon(QIcon(siteIcon))
        self.siteSpecific.setIconSize(QSize(30, 30))
        self.siteSpecific.setMinimumSize(QSize(120, 50))
        self.siteSpecific.setStyleSheet("text-align:left;")
        self.verticalLayout.addWidget(self.siteSpecific)

        self.groupProteinThings = QGroupBox(self)
        self.groupProteinThings.setObjectName("groupProteinThings")
        sizePolicy2.setHeightForWidth(self.groupProteinThings.sizePolicy().hasHeightForWidth())
        self.groupProteinThings.setSizePolicy(sizePolicy2)
        self.groupProteinThings.setMinimumSize(QSize(120, 50))

        self.formlayout1 = QFormLayout(self.groupProteinThings)
        self.formlayout1.setObjectName("formlayout1")
        self.formlayout1.setContentsMargins(0, 0, 0, 0)

        self.bindingPred = QPushButton("  Binding site prediction", self.groupProteinThings)
        self.bindingPred.setObjectName("bindingPred")
        self.bindingPred.setCheckable(True)
        bindingIcon = QPixmap(":/Binding.png")
        self.bindingPred.setIcon(QIcon(bindingIcon))
        self.bindingPred.setIconSize(QSize(30, 30))
        self.bindingPred.setMinimumSize(QSize(120, 50))
        self.bindingPred.clicked.connect(self.connectbindingPred)
        self.bindingPred.setStyleSheet("text-align:left;")
        
        self.fixed = QPushButton("Fix PDB", self.groupProteinThings)
        self.fixed.setObjectName("fix")
        self.fixed.setMinimumSize(QSize(0, 50))
        self.fixed.setCheckable(True)
        sizePolicy6.setHeightForWidth(self.fixed.sizePolicy().hasHeightForWidth())
        self.fixed.setSizePolicy(sizePolicy6)

        self.formlayout1.addRow(self.bindingPred,self.fixed)

        self.verticalLayout.addWidget(self.groupProteinThings)

        self.tabWidget1.addTab(self.dock, "Dock")

        self.prep = QWidget()
        self.prep.setObjectName("prep")
        self.gridLayout2 = QGridLayout(self.prep)
        self.gridLayout2.setObjectName("gridLayout2")
        
        self.inspectTarget = QPushButton("  Inspect target\n  proteins", self.prep)
        self.inspectTarget.setObjectName("inspectTarget")
        self.inspectTarget.setMinimumSize(QSize(0, 50))
        inspectIcon = QPixmap(":/Inspect.png")
        self.inspectTarget.setIcon(QIcon(inspectIcon))
        self.inspectTarget.setIconSize(QSize(30, 30))
        self.inspectTarget.setStyleSheet("text-align:left;")
        self.inspectTarget.clicked.connect(self.inspectPDB)
        sizePolicy6.setHeightForWidth(self.inspectTarget.sizePolicy().hasHeightForWidth())
        self.inspectTarget.setSizePolicy(sizePolicy6)
        self.gridLayout2.addWidget(self.inspectTarget, 0, 0, 1, 1)

        self.prepTarget = QPushButton("  Prepare\n  Targets", self.prep)
        self.prepTarget.setObjectName("prepTarget")
        targetIcon = QPixmap(":/Target.png")
        self.prepTarget.setIcon(QIcon(targetIcon))
        self.prepTarget.setIconSize(QSize(30, 30))
        self.prepTarget.setStyleSheet("text-align:left;")
        self.prepTarget.setMinimumSize(QSize(0, 50))
        self.prepTarget.setSizePolicy(sizePolicy6)
        self.prepTarget.clicked.connect(self.connectTargetPrep)
        self.gridLayout2.addWidget(self.prepTarget, 0, 2, 1, 1)

        self.prepLigand = QPushButton("  Prepare\n  ligand(s)", self.prep)
        self.prepLigand.setObjectName("prepLigand")
        ligandIcon = QPixmap(":/Ligand.png")
        self.prepLigand.setIcon(QIcon(ligandIcon))
        self.prepLigand.setIconSize(QSize(30, 30))
        self.prepLigand.setMinimumSize(QSize(0, 50))
        self.prepLigand.setStyleSheet("text-align:left;")
        self.prepLigand.clicked.connect(self.connectLigandPrep)
        sizePolicy6.setHeightForWidth(self.prepLigand.sizePolicy().hasHeightForWidth())
        self.prepLigand.setSizePolicy(sizePolicy6)
        self.gridLayout2.addWidget(self.prepLigand, 1, 0, 1, 1)
        
        self.conf = QPushButton("  Generate\n  config files", self.prep)
        self.conf.setObjectName("conf")
        confIcon = QPixmap(":/Config.png")
        self.conf.setIcon(QIcon(confIcon))
        self.conf.setIconSize(QSize(30, 30))
        self.conf.setMinimumSize(QSize(0, 50))
        self.conf.setStyleSheet("text-align:left;")
        sizePolicy6.setHeightForWidth(self.conf.sizePolicy().hasHeightForWidth())
        self.conf.setSizePolicy(sizePolicy6)
        self.conf.clicked.connect(self.connectConf)
        self.gridLayout2.addWidget(self.conf, 1, 2, 1, 1)

        self.splitLib = QPushButton("  Split a\n  Library file", self.prep)
        self.splitLib.setObjectName("splitLib")
        libIcon = QPixmap(":/Split.png")
        self.splitLib.setIcon(QIcon(libIcon))
        self.splitLib.setIconSize(QSize(30, 30))
        self.splitLib.setMinimumSize(QSize(0, 50))
        self.splitLib.setStyleSheet("text-align:left;")
        self.splitLib.clicked.connect(self.connectSplitLib)
        self.gridLayout2.addWidget(self.splitLib, 2, 0, 1, 1)

        self.convert = QPushButton(" Convert files", self.prep)
        self.convert.setObjectName("convert")
        convertIcon = QPixmap(":/Convert.png")
        self.convert.setIcon(QIcon(convertIcon))
        self.convert.setIconSize(QSize(30, 30))
        self.convert.setMinimumSize(QSize(0, 50))
        self.convert.setStyleSheet("text-align:left;")
        self.convert.clicked.connect(self.connectConvert)
        self.gridLayout2.addWidget(self.convert, 2, 2, 1, 1)

        self.tabWidget1.addTab(self.prep, "Prepare")
        
        self.horizontalLayout1.addWidget(self.tabWidget1)

        self.tabWidget2 = QTabWidget(self.groupTabs)
        self.tabWidget2.setObjectName("tabWidget2")
        sizePolicy5.setHeightForWidth(self.tabWidget2.sizePolicy().hasHeightForWidth())
        self.tabWidget2.setSizePolicy(sizePolicy5)
        self.tabWidget2.setUsesScrollButtons(False)
        self.analyze = QWidget()
        self.analyze.setObjectName("analyze")
        self.gridLayout3 = QGridLayout(self.analyze)
        self.gridLayout3.setObjectName("gridLayout3")
        
        self.hit = QPushButton("  Hit identifier", self.analyze)
        self.hit.setObjectName("hit")
        hitIcon = QPixmap(":/Hit.png")
        self.hit.setIcon(QIcon(hitIcon))
        self.hit.setIconSize(QSize(30, 30))
        sizePolicy5.setHeightForWidth(self.hit.sizePolicy().hasHeightForWidth())
        self.hit.setSizePolicy(sizePolicy5)
        self.hit.clicked.connect(self.connectHit)
        self.hit.setStyleSheet("text-align:left;")
        self.gridLayout3.addWidget(self.hit, 0, 2, 1, 1)

        self.interact = QPushButton("  Target-Ligand\n  interactions",self.analyze)
        self.interact.setObjectName("interact")
        interactIcon = QPixmap(":/Interact.png")
        self.interact.setIcon(QIcon(interactIcon))
        self.interact.setIconSize(QSize(30, 30))
        sizePolicy5.setHeightForWidth(self.interact.sizePolicy().hasHeightForWidth())
        self.interact.setSizePolicy(sizePolicy5)
        self.interact.clicked.connect(self.connectInteract)
        self.gridLayout3.addWidget(self.interact, 0, 0, 1, 1)

        self.vis = QPushButton("Visualizer and\nComplex Maker", self.analyze)
        self.vis.setObjectName("vis")
        sizePolicy5.setHeightForWidth(self.vis.sizePolicy().hasHeightForWidth())
        self.vis.setSizePolicy(sizePolicy5)
        self.vis.clicked.connect(self.connectVis)
        self.gridLayout3.addWidget(self.vis, 1, 0, 1, 1)
        
        self.vinaSplit = QPushButton("  Vina splitter", self.analyze)
        self.vinaSplit.setObjectName("vinaSplit")
        SplitterIcon = QPixmap(":/Splitter.png")
        self.vinaSplit.setIcon(QIcon(SplitterIcon))
        self.vinaSplit.setIconSize(QSize(30, 30))
        sizePolicy5.setHeightForWidth(self.vinaSplit.sizePolicy().hasHeightForWidth())
        self.vinaSplit.setSizePolicy(sizePolicy5)
        self.vinaSplit.clicked.connect(self.vinaSplitter)
        self.vinaSplit.setStyleSheet("text-align:left;")
        self.gridLayout3.addWidget(self.vinaSplit, 1, 2, 1, 1)

        self.admet = QPushButton("  ADMET filter", self.analyze)
        self.admet.setObjectName("admet")
        admetIcon = QPixmap(":/Admet.png")
        self.admet.setIcon(QIcon(admetIcon))
        self.admet.setIconSize(QSize(30, 30))
        sizePolicy5.setHeightForWidth(self.admet.sizePolicy().hasHeightForWidth())
        self.admet.setSizePolicy(sizePolicy5)
        self.admet.clicked.connect(self.connectAdmet)
        self.admet.setStyleSheet("text-align:left;")
        self.gridLayout3.addWidget(self.admet, 2, 0, 1, 1)

        self.ro5 = QPushButton("  Lipinski Ro5", self.analyze)
        self.ro5.setObjectName("ro5")
        ro5Icon = QPixmap(":/Ro5.png")
        self.ro5.setIcon(QIcon(ro5Icon))
        self.ro5.setIconSize(QSize(30, 30))
        sizePolicy5.setHeightForWidth(self.ro5.sizePolicy().hasHeightForWidth())
        self.ro5.setSizePolicy(sizePolicy5)
        self.ro5.clicked.connect(self.connectRo5)
        self.ro5.setStyleSheet("text-align:left;")
        self.gridLayout3.addWidget(self.ro5, 2, 2, 1, 1)

        self.tabWidget2.addTab(self.analyze, "Analyze")
        
        self.settings = QWidget()
        self.settings.setObjectName("settings")
        self.gridLayout4 = QGridLayout(self.settings)
        self.gridLayout4.setObjectName("gridLayout4")

        self.conformerLabel = QLabel(self.settings)
        self.conformerLabel.setObjectName("conformerLabel")
        self.conformerLabel.setText("No. of conformers:\n(Default: 9)")
        self.conformerLabel.setDisabled(True)
        sizePolicy6.setHeightForWidth(self.conformerLabel.sizePolicy().hasHeightForWidth())
        self.conformerLabel.setSizePolicy(sizePolicy6)
        self.gridLayout4.addWidget(self.conformerLabel, 0, 0, 1, 2)

        self.conformers = QLineEdit("9", self.settings)
        self.conformers.setObjectName("conformers")
        self.conformers.setMinimumSize(QSize(0, 50))
        sizePolicy6.setHeightForWidth(self.conformers.sizePolicy().hasHeightForWidth())
        self.conformers.setSizePolicy(sizePolicy6)
        self.gridLayout4.addWidget(self.conformers, 0, 2, 1, 1)

        self.dpiLabel = QLabel(self.settings)
        self.dpiLabel.setObjectName("dpiLabel")
        self.dpiLabel.setText("Please select the DPI:\n(Default: 300)")
        self.dpiLabel.setDisabled(True)
        sizePolicy6.setHeightForWidth(self.dpiLabel.sizePolicy().hasHeightForWidth())
        self.dpiLabel.setSizePolicy(sizePolicy6)
        self.gridLayout4.addWidget(self.dpiLabel, 1, 0, 1, 2)

        self.dpi = QLineEdit("300", self.settings)
        self.dpi.setObjectName("dpi")
        self.dpi.setMinimumSize(QSize(0, 50))
        sizePolicy6.setHeightForWidth(self.dpi.sizePolicy().hasHeightForWidth())
        self.dpi.setSizePolicy(sizePolicy6)
        self.gridLayout4.addWidget(self.dpi, 1, 2, 1, 1)

        self.theme = QPushButton("  Light\n  Mode", self.settings)
        self.theme.setObjectName("vinaSplit")
        self.theme.setMinimumSize(QSize(0, 50))
        sizePolicy6.setHeightForWidth(self.theme.sizePolicy().hasHeightForWidth())
        self.theme.setSizePolicy(sizePolicy6)
        self.theme.setIcon(QIcon(":/Sun.png"))
        self.theme.setToolTip("Switch to Dark Mode")
        self.theme.clicked.connect(self.toggle_dark_mode)
        self.gridLayout4.addWidget(self.theme, 2, 0, 1, 1)

        self.resume = QPushButton(self.analyze)
        self.resume.setObjectName("resume")
        resumeIcon = QPixmap(":/Resume.png")
        self.resume.setIcon(QIcon(resumeIcon))
        self.resume.setIconSize(QSize(30, 30))
        self.resume.setStyleSheet("text-align:left; padding-right: 1px;")
        self.resume.setCheckable(True)
        self.resume.setMinimumSize(QSize(50, 50))
        self.resume.setMaximumSize(QSize(50, 50))
        self.resume.setToolTip("Click to resume docking after the last job")
        sizePolicy4.setHeightForWidth(self.resume.sizePolicy().hasHeightForWidth())
        self.resume.setSizePolicy(sizePolicy4)
        self.gridLayout4.addWidget(self.resume, 2, 1, 1, 1)

        self.log = QPushButton("  Enable Log", self.analyze)
        self.log.setObjectName("log")
        logIcon = QPixmap(":/Log.png")
        self.log.setIcon(QIcon(logIcon))
        self.log.setIconSize(QSize(30, 30))
        self.log.setStyleSheet("text-align:left;")
        self.log.setMinimumSize(QSize(0, 50))
        self.log.setCheckable(True)
        sizePolicy6.setHeightForWidth(self.log.sizePolicy().hasHeightForWidth())
        self.log.setSizePolicy(sizePolicy6)
        self.gridLayout4.addWidget(self.log, 2, 2, 1, 1)

        self.horizontalLayout1.addWidget(self.tabWidget2)
        self.tabWidget2.addTab(self.settings, "Settings")

        work_space_layout.addWidget(self.groupTabs, 3, 0, 1, 3)

        self.progress_bar = QProgressBar()
        self.progress_bar.setValue(0)
        self.progress_bar.setTextVisible(False)
        self.progress_bar.setMaximumSize(self.width(), 31)
        work_space_layout.addWidget(self.progress_bar, 4, 0, 1, 3)
        self.progressLabel = QLabel()
        self.progressLabel.setStyleSheet("text-align:center;color:darkred;font-size:16px;"
        "font-weight:extra-bold;")
        work_space_layout.addWidget(self.progressLabel, 4, 0, 1, 3)

        self.groupStart = QGroupBox(self)
        self.groupStart.setObjectName("groupLine")
        sizePolicy3.setHeightForWidth(self.groupStart.sizePolicy().hasHeightForWidth())
        self.groupStart.setSizePolicy(sizePolicy3)
        self.verticalLayout2 = QVBoxLayout(self.groupStart)
        self.verticalLayout2.setObjectName("verticalLayout2")
        self.verticalLayout2.setAlignment(Qt.AlignCenter) 
        self.verticalLayout2.setContentsMargins(0, 0, 0, 0)

        self.startButton = QPushButton("  S T A R T", self)
        self.startButton.setObjectName("startButton")
        self.startButton.setMinimumSize(QSize(300, 70))
        startIcon = QPixmap(":/Start.png")
        self.startButton.setIcon(QIcon(startIcon))
        self.startButton.setIconSize(QSize(50, 50))
        sizePolicy4.setHeightForWidth(self.startButton.sizePolicy().hasHeightForWidth())
        self.startButton.setSizePolicy(sizePolicy4)
        self.startButton.clicked.connect(self.startEngine)
        self.verticalLayout2.addWidget(self.startButton, alignment=Qt.AlignVCenter)
        
        work_space_layout.addWidget(self.groupStart, 5, 1, 1, 1)

        central_widget_layout.addLayout(work_space_layout)
        central_widget.setLayout(central_widget_layout)
        self.setCentralWidget(central_widget)

    def browseTargets(self):
        self.browseTargetsLabel.setText(QFileDialog.getExistingDirectory(None, 'Select Target directory'))
    def browseLigands(self):
        self.browseLigandsLabel.setText(QFileDialog.getExistingDirectory(None, 'Select Ligand directory'))

    def showResidueDialog(self):
        text, ok = QInputDialog.getMultiLineText(
            self,
            "Enter Binding‑Site Residues",
            "Comma‑separated (e.g. A256 [Hemoglobin], B102 [Myoglobin], 23C [Kinase]):",
            ""
        )
        if not ok or not text.strip():
            ID_utils.displayWarning(self, "Input required", "You must enter at least one residue for each target.")
            # re‑ask or send empty string
            text = ""
        # finally, send it back into the worker’s thread
        self.startIt.residuesReady.emit(text)

    def showResidueDialogTarget(self):
        text, ok = QInputDialog.getMultiLineText(
            self,
            "Enter Binding‑Site Residues",
            "Comma‑separated (e.g. A256 [Hemoglobin], B102 [Myoglobin], 23C [Kinase]):",
            ""
        )
        if not ok or not text.strip():
            ID_utils.displayWarning(self, "Input required", "You must enter at least one residue for each target.")
            # re‑ask or send empty string
            text = ""
        # finally, send it back into the worker’s thread
        self.targetPrep.residuesReady.emit(text)

    def showResidueDialogConf(self):
        text, ok = QInputDialog.getMultiLineText(
            self,
            "Enter Binding‑Site Residues",
            "Comma‑separated (e.g. A256 [Hemoglobin], B102 [Myoglobin], 23C [Kinase]):",
            ""
        )
        if not ok or not text.strip():
            ID_utils.displayWarning(self, "Input required", "You must enter at least one residue for each target.")
            # re‑ask or send empty string
            text = ""
        # finally, send it back into the worker’s thread
        self.config.residuesReady.emit(text)

    def connectAbout(self):
        self.about = AboutMultiUI()
        self.about.show()

    def changedIndex(self):
        if self.comboBox.currentText() == "AutoDock-GPU" or self.comboBox.currentText() == "AutoDock Vina":
            self.log.setDisabled(True)
            self.conf.setDisabled(True)
        else:
            self.log.setDisabled(False)
            self.log.setCheckable(True)
            self.conf.setDisabled(False)

    def startEngine(self):
        targetLabel = self.browseTargetsLabel.text()
        ligandLabel = self.browseLigandsLabel.text()
        if targetLabel == "Browse target directory..." or ligandLabel == "Browse ligand directory...":
            ID_utils.displayError("❌ Please select proper working directories.")
        else:
            self.show_toast([3000, "Starting the docking process..."])
            self.startIt = multiCalculateS1()
            self.startIt.finished.connect(self.startIt.deleteLater)
            self.startIt.setTargetWD(targetLabel)
            self.startIt.checkVals([
                self.siteSpecific.isChecked(), 
                self.bindingPred.isChecked()])
            if self.comboBox.currentText() == "AutoDock-GPU":
                self.startIt.dockMethod("GPU")
            else:
                self.startIt.dockMethod("CPU")
            self.startIt.requestResidueDialog.connect(self.showResidueDialog)
            self.startIt.toastSignal.connect(self.show_toast)
            self.startIt.finishedSignal.connect(self.startEngineE1)
            self.startIt.start()
            self.startButton.setEnabled(False)

    def startEngineE1(self, text):
        if "❌" in text[0]:
            ID_utils.displayError(text[0])
            self.startButton.setEnabled(True)
        elif "✅" in text[0]:
            self.continue1 = startMultiTargetDocking()
            self.continue1.finished.connect(self.continue1.deleteLater)
            self.continue1.checkVals([
                self.resume.isChecked(),
                self.log.isChecked()
            ])
            self.continue1.dockMethod(self.comboBox.currentText())
            self.continue1.setTargetWD(self.browseTargetsLabel.text())
            self.continue1.setLigandWD(self.browseLigandsLabel.text())
            self.continue1.setNConformers(self.conformers.text())
            self.continue1.toastSignal.connect(self.show_toast)
            self.continue1.progressBarSignal.connect(self.update_progress_bar)
            self.continue1.progressLabelSignal.connect(self.update_progress_label)        
            self.continue1.finishedSignal.connect(self.continueEngine)
            self.continue1.start()
        else:
            self.continue1 = multiCalculateS2()
            self.continue1.finished.connect(self.continue1.deleteLater)
            self.continue1.completeFile(text)
            self.continue1.finishedSignal.connect(self.startEngineE2)
            self.continue1.start()

    def startEngineE2(self, text):
        if "❌" in text[0]:
            ID_utils.displayError(text[0])
            self.startButton.setEnabled(True)
        else:
            self.continue2 = multiCalculateS3()
            self.continue2.finished.connect(self.continue2.deleteLater)
            self.continue2.setTargetWD(self.browseTargetsLabel.text())
            self.continue2.allThings(text)
            if self.comboBox.currentText() == "AutoDock-GPU":
                self.continue2.dockMethod("GPU")
            else:
                self.continue2.dockMethod("CPU")
            self.continue2.fixPDB(self.fixed.isChecked())
            self.continue2.toastSignal.connect(self.show_toast)
            self.continue2.finishedSignal.connect(self.startEngineE3)
            self.continue2.start()

    def startEngineE3(self, text):
        if "❌" in text:
            ID_utils.displayError(text)
            self.startButton.setEnabled(True)
        else:
            self.continue3 = multiCalculateS4()
            self.continue3.finished.connect(self.continue3.deleteLater)
            self.continue3.toastSignal.connect(self.show_toast)
            self.continue3.setLigandWD(self.browseLigandsLabel.text())
            self.continue3.finishedSignal.connect(self.startEngineE4)
            self.continue3.start()

    def startEngineE4(self, text):
        if "❌" in text:
            ID_utils.displayError(text)
            self.startButton.setEnabled(True)
        else:
            self.continue4 = startMultiTargetDocking()
            self.continue4.finished.connect(self.continue4.deleteLater)
            self.continue4.checkVals([
                self.resume.isChecked(),
                self.log.isChecked()
            ])
            self.continue4.toastSignal.connect(self.show_toast)
            self.continue4.dockMethod(self.comboBox.currentText())
            self.continue4.setTargetWD(self.browseTargetsLabel.text())
            self.continue4.setLigandWD(self.browseLigandsLabel.text())
            self.continue4.setNConformers(self.conformers.text())
            self.continue4.progressBarSignal.connect(self.update_progress_bar)
            self.continue4.progressLabelSignal.connect(self.update_progress_label)        
            self.continue4.finishedSignal.connect(self.continueEngine)
            self.continue4.start()

    def continueEngine(self, text):
        if "DONE" in text:
            self.continueIt = continueMultiID()
            self.continueIt.finished.connect(self.continueIt.deleteLater)
            self.continueIt.dockMethod(text)
            self.continueIt.toastSignal.connect(self.show_toast)
            self.continueIt.setTargetWD(self.browseTargetsLabel.text())
            self.continueIt.setLigandWD(self.browseLigandsLabel.text())
            self.continueIt.setDPI(self.dpi.text())
            self.continueIt.setLogging(self.log.isChecked())
            self.continueIt.progressBarSignal.connect(self.update_progress_bar)
            self.continueIt.setBlind(self.blind.isChecked())
            self.continueIt.setBinding(self.bindingPred.isChecked())
            self.continueIt.setSiteSpecific(self.siteSpecific.isChecked())
            self.continueIt.progressLabelSignal.connect(self.update_progress_label)        
            self.continueIt.finishedSignal.connect(self.showFinished)
            self.continueIt.start()
        elif "⚠️" in text:
            ID_utils.displayWarning(text)
            self.startButton.setEnabled(True)
        else:
            ID_utils.displayError(text)
            self.startButton.setEnabled(True)

    def update_progress_bar(self, val):
        self.progress_bar.setValue(val)

    def update_progress_label(self, text):
        self.progressLabel.setText(text)

    def showFinished(self, text):
        if "❌" in text[0]:
            ID_utils.displayError(text[0])
            self.startButton.setEnabled(True)
        else:
            if len(text) > 1:
                self.progress_bar.setValue(100)
                self.progressLabel.setText("<b><center>Progress: 100%")
                msg = f"<p align=center style=\"color:darkblue;font-size:20px;\"><b>Thank you for using InstaDock-v2!!</b></p><p align=center>The Virtual Screening process has finished!<br />The results of the entire process can be found under the 'Results' folder.<br /><span style=\"color:darkred;\">NOTE: For viewing the receptor-ligand interactions of the generated docked conformers, you can use the \"Visualizer and Complex Maker\" option under the \"Analyze\" section.</span></p>"
                ID_utils.displayComplete(msg)
                self.startButton.setEnabled(True)

    def connectBlind(self):
        if self.blind.isChecked():
            self.siteSpecific.setChecked(False)
        else:
            self.siteSpecific.setChecked(True)
            self.bindingPred.setChecked(False)

    def connectSiteSpecific(self):
        if self.siteSpecific.isChecked():
            self.blind.setChecked(False)
            self.bindingPred.setChecked(False)
        else:
            self.blind.setChecked(True)
            self.bindingPred.setChecked(True)
            
    def connectbindingPred(self):
        if self.bindingPred.isChecked():
            self.siteSpecific.setChecked(False)
        else:
            if self.blind.isChecked():
                pass
            elif self.blind.isChecked() == False and self.siteSpecific.isChecked() == False:
                self.blind.setChecked(True)
            else:
                self.siteSpecific.setChecked(True)

    def inspectPDB(self):
        self.inspectTarget.setEnabled(False)
        if self.comboBox.currentText() == "AutoDock-GPU":
            method = "GPU"
        else:
            method = "CPU"
        self.inspectMessage, self.inspectTitle = inspectProtein(self.browseTargetsLabel.text(), method)
        if self.inspectTitle == "OK":
            ID_utils.displayMess(self.inspectMessage)
        else:
            ID_utils.displayWarning(self.inspectMessage)
        self.inspectTarget.setEnabled(True)

    def connectTargetPrep(self):
        self.prepTarget.setEnabled(False)
        self.targetPrep = prepareTargetM1()
        self.targetPrep.finished.connect(self.targetPrep.deleteLater)
        self.targetPrep.setWD(self.browseTargetsLabel.text())
        if self.comboBox.currentText() == "AutoDock-GPU":
            self.targetPrep.dockMethod("GPU")
        else:
            self.targetPrep.dockMethod("CPU")
        self.targetPrep.checkVals([
            self.siteSpecific.isChecked(), 
            self.bindingPred.isChecked(),
        ])
        self.targetPrep.fixPDB(self.fixed.isChecked())
        self.targetPrep.requestResidueDialog.connect(self.showResidueDialogTarget)
        self.targetPrep.toastSignal.connect(self.show_toast)
        self.targetPrep.finishedSignal.connect(self.connectTargetPrepS2)
        self.targetPrep.start()

    def connectTargetPrepS2(self, text):
        if "❌" in text[0]:
            ID_utils.displayError(text[0])
            self.prepTarget.setEnabled(True)
        else:
            self.targetPrep = prepareTargetM2()
            self.targetPrep.finished.connect(self.targetPrep.deleteLater)
            self.targetPrep.totalFile(text)
            self.targetPrep.finishedSignal.connect(self.connectTargetPrepS3)
            self.targetPrep.toastSignal.connect(self.show_toast)
            self.targetPrep.start()
        
    def connectTargetPrepS3(self, text):
        if "❌" in text[0]:
            ID_utils.displayError(text[0])
            self.prepTarget.setEnabled(True)
        else:
            self.targetPrep = prepareTargetM3()
            self.targetPrep.finished.connect(self.targetPrep.deleteLater)
            self.targetPrep.allThings(text)
            self.targetPrep.setTargetWD(self.browseTargetsLabel.text())
            if self.comboBox.currentText() == "AutoDock-GPU":
                self.targetPrep.dockMethod("GPU")
            else:
                self.targetPrep.dockMethod("CPU")
            self.targetPrep.fixPDB(self.fixed.isChecked())
            self.targetPrep.toastSignal.connect(self.show_toast)
            self.targetPrep.finishedSignal.connect(self.finishPreppingTarget)
            self.targetPrep.start()
            
    def finishPreppingTarget(self, text):
        self.prepTarget.setEnabled(True)
        if "❌" in text:
            ID_utils.displayError(text)
        else:
            ID_utils.displayMessage(text)

    def connectLigandPrep(self):
        self.prepLigand.setEnabled(False)
        self.ligPrep = prepareLigand()
        self.ligPrep.finished.connect(self.ligPrep.deleteLater)
        self.ligPrep.setLigandWD(self.browseLigandsLabel.text())
        self.ligPrep.finishedSignal.connect(self.finishPreppingLigand)
        self.ligPrep.start()

    def finishPreppingLigand(self, text):
        self.prepLigand.setEnabled(True)
        if "❌" in text:
            ID_utils.displayError(text)
        else:
            ID_utils.displayMessage(text)

    def connectConf(self):
        self.conf.setEnabled(False)
        self.config = prepareConfigM1()
        self.config.finished.connect(self.config.deleteLater)
        self.config.setTargetWD(self.browseTargetsLabel.text())
        self.config.dockMethod(self.comboBox.currentText())
        self.config.requestResidueDialog.connect(self.showResidueDialogConf)
        self.config.checkVals([self.siteSpecific.isChecked(), self.bindingPred.isChecked()])
        self.config.toastSignal.connect(self.show_toast)
        self.config.finishedSignal.connect(self.connectConfS2)
        self.config.start()

    def connectConfS2(self, text):
        if "⚠️" in text[0]:
            ID_utils.displayWarning(text[0])
            self.conf.setEnabled(True)
        elif "❌" in text[0]:
            ID_utils.displayError(text[0])
            self.conf.setEnabled(True)
        else:
            self.config = prepareConfigM2()
            self.config.finished.connect(self.config.deleteLater)
            self.config.totalFile(text)
            self.config.finishedSignal.connect(self.connectConfS3)
            self.config.start()

    def connectConfS3(self, text):
        if "❌" in text[0]:
            ID_utils.displayError(text[0])
            self.conf.setEnabled(True)
        else:
            self.config = prepareConfigM3()
            self.config.finished.connect(self.config.deleteLater)
            self.config.setTargetWD(self.browseTargetsLabel.text())
            self.config.allThings(text)
            self.config.finishedSignal.connect(self.finishPreppingConfig)
            self.config.start()
        
    def finishPreppingConfig(self, text):
        self.conf.setEnabled(True)
        ID_utils.displayMessage(text)

    def connectSplitLib(self):
        self.splitLib = ID_utils.librarySplitter()
        self.splitLib.show()

    def vinaSplitter(self):
        self.splitLig = ID_utils.ligandSplitter()
        self.splitLig.show()

    def connectHit(self):
        self.tophits = ID_utils.topHits()
        self.tophits.show()

    def connectVis(self):
        self.vis = idVisualizer()
        self.vis.show()

    def connectInteract(self):
        self.interact = idTwoDInteract()
        self.interact.show()

    def connectAdmet(self):
        self.admet = admet()
        self.admet.show()

    def connectRo5(self):
        self.ro5 = ro5()
        self.ro5.show()

    def connectConvert(self):
        #Pass a disclaimer stating that these files might not be as compatible with the entire docking pipeline (since they are created using OpenBabel)
        ID_utils.displayWarning("Compatibility Warning: Files created with OpenBabel may cause unexpected\nissues in later docking stages.")
        self.converting = convertFiles()
        self.converting.setTargetWD(self.browseTargetsLabel.text())
        self.converting.setLigandWD(self.browseLigandsLabel.text())
        self.converting.toastSignal.connect(self.show_toast)
        self.converting.finishedSignal.connect(self.finishConvert)
        self.converting.start()
        self.convert.setEnabled(False)

    def finishConvert(self, text):
        self.convert.setEnabled(True)
        if "❌" in text:
            ID_utils.displayError(text)
        else:
            ID_utils.displayMessage(text)

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

    def load_fonts_from_dir(self, directory):
        families = set()
        for fi in QDir(directory).entryInfoList(["*.ttf"]):
            _id = QFontDatabase.addApplicationFont(fi.absoluteFilePath())
            families |= set(QFontDatabase.applicationFontFamilies(_id))
        return families

    def toggle_dark_mode(self):
        if self.is_dark_mode:
            # Switch to light mode
            self.setStyleSheet('''
            QMainWindow {
                background: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:1, 
                stop:0 #F5F7FA, 
                stop:0.5 #E4E7EB, 
                stop:1 #D9DEE4); /* Light and minimalistic gradient */
                color: #000000; /* Default dark text */
            }
            
            QLabel{
                background: transparent;
                font-size: 14px;
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
                width: 153px;
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
                width: 153px;
                height: 35px;
                background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, 
                            stop:0 #E0E0E0, /* Slightly darker hover color */
                            stop:0.5 #D4D4D4, /* Smooth transition */
                            stop:1 #C8C8C8); /* End hover gradient */
                border: 1px solid #37cbe6; /* Highlighted border on hover */
                border-bottom: 0px;
            }

            QTabBar::tab:selected {
                width: 153px;
                height: 35px;
                background: white; /* Keeping the selected tab background white */
                border: 1px solid #37cbe6; /* Slight emphasis on the selected tab with matching border */
                border-bottom: 0px;
                color: black; /* Dark text for readability */
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
            '''
            )
            self.theme.setText("  Light\n  Mode")
            self.theme.setIcon(QIcon(":/Sun.png"))
            self.theme.setToolTip("Switch to Dark Mode")
            self.is_dark_mode = False
        else:
            # Switch to dark mode
            self.setStyleSheet('''
            QMainWindow {
                background: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:1, 
                stop:0 #1E1E1E, 
                stop:0.5 #2A2A2A, 
                stop:1 #3A3A3A); /* Dark, minimalistic gradient for background */
                color: #E0E0E0; /* Light text for contrast */
            }

            QLabel {
                background: transparent;
                font-size: 14px;
                color: #E0E0E0; /* Light text for readability */
            }

            QLabel#browseWDLabel{
                background: transparent;
                font-size: 14px;
                border: 1px solid #555555;
                border-radius: 3px;
            }
                               
            QGroupBox {
                border: 0;
                background: transparent;
            }

            QPushButton {
                font-size: 14px; /* Matches the font size */
                padding: 2px 8px; /* Inner spacing */
                border: 1px solid #555555; /* Darker border for dark mode */
                border-radius: 3px; /* Slightly rounded corners for a modern look */
                color: #E0E0E0; /* Light text for contrast */
                background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, 
                stop:0 #6A6A6A, stop:0.5 #4B4B4B, stop:0.51 #3E3E3E, stop:1 #5A5A5A); /* Inverted gradient with darker shades at the bottom */
            }

            QPushButton:hover {
                border: 1px solid #37cbe6; /* Light blue border on hover */
                background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, 
                stop:0 #4F6F7F, stop:0.5 #6A7F8E, stop:1 #4F6F7F); /* Lighter blue gradient on hover with darker bottom */
            }

            QPushButton:pressed {
                padding: 3px 9px 4px 11px; /* Slight padding shift for pressed effect */
                border: 1px solid #37cbe6; /* Blue border when pressed */
                background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, 
                stop:0 #3A4F5A, stop:0.5 #4F6877, stop:1 #3E4F5C); /* Darker gradient when pressed with darker bottom */
            }

            QPushButton:checked {
                padding: 3px 9px 4px 11px; /* Slight padding shift for pressed effect */
                border: 1px solid #37cbe6; /* Blue border when checked */
                background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, 
                stop:0 #3A4F5A, stop:0.5 #4F6877, stop:1 #3E4F5C); /* Darker gradient when checked with darker bottom */
            }


            QPushButton#startButton {
                font-size: 35px;
            }

           QComboBox {
                font-size: 14px; /* Matches the font size */
                padding: 4px 8px; /* Inner spacing for better usability */
                border: 1px solid #555555; /* Darker border for dark mode */
                border-radius: 3px; /* Rounded corners for modern look */
                background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, 
                stop:0 #2C2C2C, stop:0.5 #333333, stop:1 #444444); /* Dark gradient background */
                color: #E0E0E0; /* Light text for contrast */
            }

            QLineEdit {
                font-size: 14px; /* Matches the font size */
                padding: 1px 1px; /* Inner spacing */
                color: #E0E0E0; /* Light text for contrast */
                background: #2A2A2A; /* Dark background to match dark mode */
                border: 1px solid #A9A9A9;
                border-radius: 3px;
            }

            QTabWidget::pane {
                border: 1px solid #A9A9A9; /* Matching border color to QPushButton */
                background: #3A3A3A; /* Subtle dark background */
            }

            QTabBar::tab {
                width: 153px;
                height: 35px;
                font-size: 14px;
                color: #E0E0E0; /* Light text for contrast */
                background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, 
                stop:0 #3A3A3A, stop:0.5 #444444, stop:1 #505050); /* Dark gradient for tabs */
                border: 1px solid #A9A9A9; /* Border to match QPushButton */
                border-bottom: 0px;
            }

            QTabBar::tab:hover {
                width: 153px;
                height: 35px;
                background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, 
                stop:0 #4A4A4A, stop:0.5 #555555, stop:1 #666666); /* Slightly lighter on hover */
                border: 1px solid #37cbe6; /* Highlighted border on hover */
                border-bottom: 0px;
            }

            QTabBar::tab:selected {
                width: 153px;
                height: 35px;
                background: #2A2A2A; /* Dark selected tab */
                border: 1px solid #37cbe6; /* Slight emphasis on the selected tab with matching border */
                border-bottom: 0px;
                color: #E0E0E0; /* Light text for readability */
            }

            QProgressBar { 
                border: 2px solid #37cbe6; 
                border-radius: 3px; 
                text-align: center;
                background: #444444; /* Dark background for the progress bar */
            }

            QProgressBar::chunk { 
                background-color: #5A9CCF; /* Lighter blue for the progress chunk */
                width: 2px; 
                margin: 0.5px; 
                border-radius: 1.5px; 
            }
    
            ''')
            self.theme.setText("  Dark\n  Mode")
            self.theme.setIcon(QIcon(":/Moon.png"))
            self.theme.setToolTip("Switch to Light Mode")
            self.is_dark_mode = True


if hasattr(Qt, 'AA_EnableHighDpiScaling'):
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
if hasattr(Qt, 'AA_UseHighDpiPixmaps'):
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)
