import warnings
warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=DeprecationWarning)

import os
import sys
import IDv2Resource
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtGui import QIcon, QPixmap
from PyQt5.QtCore import QSize, Qt, QSize 
from ID_Single_Main import SingleTargetUI
from ID_Multi_Main import MultiTargetUI
from QSwitchControl import SwitchControl
from ID_Ligand_Downloader import DownloadUI
from PyQt5.QtWidgets import (
    QPushButton, QGridLayout, QHBoxLayout, QLabel, QMainWindow, QToolButton,
    QVBoxLayout, QWidget, QSpacerItem, QSizePolicy, QGroupBox, QDesktopWidget, QApplication)

if getattr(sys, 'frozen', False):
    # When PyInstaller bundles everything, sys._MEIPASS contains the temporary folder path.
    basedir = sys._MEIPASS
    rdkit_data_path = os.path.join(basedir, "rdkit", "Chem")
    os.environ["RDKIT_DATA"] = rdkit_data_path

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


class StartingUI(QMainWindow):
    def center_above_screen(self):
        screen_geometry = QDesktopWidget().screenGeometry()
        window_geometry = self.frameGeometry()
        screen_center = screen_geometry.center()
        window_geometry.moveCenter(screen_center)
        self.move(window_geometry.topLeft())

    def __init__(self):
        super().__init__()
        self.setWindowIcon(QIcon(':/IDLogo.png'))
        self.setWindowTitle('''InstaDock-v2''')
        self.resize(520, 750)
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

        QGroupBox{
            border:0;
            background: transparent;
        }
        
        """)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sizePolicy)
        self.setMinimumSize(QSize(520, 750))
        self.setWindowFlags(Qt.FramelessWindowHint)
        central_widget = QWidget()
        central_widget.setObjectName("Container")
        self.title_bar = CustomTitleBar(self)

        central_widget_layout = QVBoxLayout()
        central_widget_layout.setContentsMargins(0, 0, 0, 0)
        central_widget_layout.setAlignment(Qt.AlignTop)

        central_widget_layout.addWidget(self.title_bar)

        work_space_layout = QGridLayout()
        work_space_layout.setContentsMargins(20, 20, 20, 20)

        sizePolicy1 = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Minimum)
        sizePolicy1.setHorizontalStretch(0)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy2 = QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)
        sizePolicy2.setHorizontalStretch(0)
        sizePolicy2.setVerticalStretch(0)

        self.titleLabel = QLabel(self)
        self.titleLabel.setObjectName("titleLabel")
        self.titleLabel.setTextFormat(Qt.RichText)  # Keeping RichText for future use if needed

        # Load and scale the image properly while maintaining aspect ratio
        pixmap = QPixmap(":/IDLogoWithText.png")
        scaled_pixmap = pixmap.scaled(
            420, 120,  # Adjust width and height as needed
            Qt.KeepAspectRatio,
            Qt.SmoothTransformation
        )
        self.titleLabel.setPixmap(scaled_pixmap)

        # Optional: Set background color and font size if needed
        # self.titleLabel.setStyleSheet("QLabel { color : #FF9D00; font-size: 50pt; }")

        sizePolicy1.setHeightForWidth(self.titleLabel.sizePolicy().hasHeightForWidth())
        self.titleLabel.setSizePolicy(sizePolicy1)
        self.titleLabel.setMinimumSize(QSize(0, 100))
        self.titleLabel.setAlignment(Qt.AlignCenter)
        
        work_space_layout.addWidget(self.titleLabel, 0, 0, 1, 2, Qt.AlignCenter)

        self.textLabel = QLabel(self)
        self.textLabel.setObjectName("textLabel")
        self.textLabel.setTextFormat(QtCore.Qt.RichText)
        self.textLabel.setText('''
<span style="font-size:14px; font-family:Arial, sans-serif;">
  <div style="text-align:center; font-weight:bold; font-size:16px; margin-bottom:1em;">
    Welcome to InstaDock-v2!
  </div>
  <p>
    We're excited to introduce InstaDock-v2, a ground-up redesign built around your feedback. This isn't just another update; it's a complete reimagining of our docking platform. With full GPU optimization, you'll enjoy dramatically faster performance on both Linux and Windows.
  </p>
  <p>
    InstaDock-v2 offers:
    <ul>
      <li>Integrated ADMET predictions for seamless early-stage assessment</li>
      <li>Advanced 3D visualization to explore binding poses in detail</li>
      <li>Comprehensive 2D interaction analysis to map key contacts</li>
    </ul>
    …and much more to support your discovery journey.
  </p>
  <p>
    Choose the mode that fits your workflow:
    <ul>
      <li><strong>Single Target Docking</strong><br>
          Perform in-depth docking studies on a single binding site—perfect for focused, high-precision analyses.
      </li>
      <li><strong>Multi Target Docking</strong><br>
          Screen your compounds against multiple proteins at once to accelerate discovery across diverse targets.
      </li>
    </ul>
  </p>
  <div style="text-align:center; margin-top:1.5em;">
    <span style="display:inline-block; padding:0.75em 1.5em; font-size:1em; font-weight:bold; border-radius:999px; background:linear-gradient(90deg, #4facfe 0%, #00f2fe 100%); color:#000; box-shadow:0 4px 6px rgba(0,0,0,0.1);">
      Thank you for choosing InstaDock-v2 to power your research!
    </span>
  </div>
</span>

''')

        self.textLabel.setWordWrap(True)
        sizePolicy2.setHeightForWidth(self.textLabel.sizePolicy().hasHeightForWidth())
        self.textLabel.setMinimumSize(QSize(0, 300))
        self.textLabel.setSizePolicy(sizePolicy2)
        self.textLabel.setAlignment(Qt.AlignJustify)

        work_space_layout.addWidget(self.textLabel, 1, 0, 15, 2, Qt.AlignCenter)

        self.tg1Group = QGroupBox(self)
        self.tg1Group.setObjectName("tg1Group")
        self.horizontalLayout1 = QtWidgets.QHBoxLayout(self.tg1Group)

        self.singleLabel = QLabel(self.tg1Group)
        self.singleLabel.setObjectName("singleLabel")
        self.singleLabel.setText("Single Target")
        sizePolicy1.setHeightForWidth(self.singleLabel.sizePolicy().hasHeightForWidth())
        self.singleLabel.setSizePolicy(sizePolicy1)
        self.horizontalLayout1.addWidget(self.singleLabel)

        self.toggle_1 = SwitchControl(bg_color="#777777", circle_color="#DDD", active_color="#500af2", animation_curve=QtCore.QEasingCurve.InOutCubic, animation_duration=100, checked=True, change_cursor=True)
        self.toggle_1.stateChanged.connect(self.on_switch1_toggled)
        self.horizontalLayout1.addWidget(self.toggle_1)

        self.tg2Group = QGroupBox(self)
        self.tg2Group.setObjectName("tg1Group")
        self.horizontalLayout2 = QtWidgets.QHBoxLayout(self.tg2Group)

        self.multiLabel = QLabel(self.tg2Group)
        self.multiLabel.setObjectName("multiLabel")
        self.multiLabel.setText("Multi Target")
        sizePolicy1.setHeightForWidth(self.multiLabel.sizePolicy().hasHeightForWidth())
        self.multiLabel.setSizePolicy(sizePolicy1)
        self.horizontalLayout2.addWidget(self.multiLabel)

        self.toggle_2 = SwitchControl(bg_color="#777777", circle_color="#DDD", active_color="#500af2", animation_curve=QtCore.QEasingCurve.InOutCubic, animation_duration=100, checked=False, change_cursor=True)
        self.toggle_2.stateChanged.connect(self.on_switch2_toggled)
        self.horizontalLayout2.addWidget(self.toggle_2)

        work_space_layout.addWidget(self.tg1Group, 17, 0, 1, 1, Qt.AlignCenter)
        work_space_layout.addWidget(self.tg2Group, 17, 1, 1, 1, Qt.AlignCenter)

        self.horizontalSpacer = QSpacerItem(40, 5, QSizePolicy.Expanding, QSizePolicy.Minimum)

        work_space_layout.addItem(self.horizontalSpacer, 18, 0, 1, 2)

        self.continueButton = QPushButton("Continue", self)
        self.continueButton.setObjectName("continueButton")
        self.continueButton.setMinimumSize(QSize(150, 50))
        self.continueButton.clicked.connect(self.connectWindow)

        work_space_layout.addWidget(self.continueButton, 19, 0, 1, 2, Qt.AlignCenter)

        self.downloadGroup = QGroupBox(self)
        self.downloadGroup.setObjectName("downloadGroup")
        self.horizontalLayout3 = QtWidgets.QHBoxLayout(self.downloadGroup)

        self.downloadLabel = QLabel(self.downloadGroup)
        self.downloadLabel.setObjectName("downloadLabel")
        self.downloadLabel.setText('''<span style="font-size:14px">
        Click here to download ligands for docking ➜
        </span>''')
        sizePolicy1.setHeightForWidth(self.downloadLabel.sizePolicy().hasHeightForWidth())
        self.downloadLabel.setSizePolicy(sizePolicy1)
        self.horizontalLayout3.addWidget(self.downloadLabel)

        self.downloadButton = QPushButton("  Download", self)
        self.downloadButton.setObjectName("downloadButton")
        downloadIcon = QPixmap(":/Download.png")
        self.downloadButton.setIcon(QIcon(downloadIcon))
        self.downloadButton.setIconSize(QSize(20, 20)) 
        self.downloadButton.clicked.connect(self.connectDownloader)
        self.horizontalLayout3.addWidget(self.downloadButton)

        work_space_layout.addWidget(self.downloadGroup, 20, 0, 1, 2, Qt.AlignJustify)

        work_space_layout.addWidget(self.textLabel, 1, 0, 1, 2, Qt.AlignCenter)
        central_widget_layout.addLayout(work_space_layout)
        central_widget.setLayout(central_widget_layout)
        self.setCentralWidget(central_widget)

    def connectWindow(self):
        if self.toggle_1.isChecked():
            self.close()
            self.singleTarget = SingleTargetUI()
            self.singleTarget.show()
        elif self.toggle_2.isChecked():
            self.close()
            self.multiTarget = MultiTargetUI()
            self.multiTarget.show()   

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

    def connectDownloader(self):
        self.downloadUI = DownloadUI()
        self.downloadUI.show()

if hasattr(QtCore.Qt, 'AA_EnableHighDpiScaling'):
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
if hasattr(QtCore.Qt, 'AA_UseHighDpiPixmaps'):
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)

if __name__ == "__main__":
    app = QApplication([])
    window = StartingUI()
    window.show()
    app.exec_()