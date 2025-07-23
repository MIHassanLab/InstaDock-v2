from PyQt5.QtGui import QIcon, QPixmap
from PyQt5.QtCore import QSize, Qt
from PyQt5.QtWidgets import ( QTabWidget, QGridLayout, QHBoxLayout, QLabel, QMainWindow, QToolButton,
    QVBoxLayout, QWidget, QSpacerItem, QSizePolicy, QDesktopWidget, 
)
    
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

class AboutMultiUI(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.setupUi(self)
    def center_above_screen(self):
        screen_geometry = QDesktopWidget().screenGeometry()
        window_geometry = self.frameGeometry()
        screen_center = screen_geometry.center()
        screen_center.setX(screen_center.x() + 200)
        screen_center.setY(screen_center.y() - 50)
        window_geometry.moveCenter(screen_center)
        self.move(window_geometry.topLeft())

    def __init__(self):
        super().__init__()
        self.setWindowIcon(QIcon(':/IDLogo.png'))
        self.setWindowTitle("IDv2 About section")
        self.resize(600, 600)
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
            width: 142px;
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
            width: 142px;
            height: 35px;
            background: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, 
                        stop:0 #E0E0E0, /* Slightly darker hover color */
                        stop:0.5 #D4D4D4, /* Smooth transition */
                        stop:1 #C8C8C8); /* End hover gradient */
            border: 1px solid #37cbe6; /* Highlighted border on hover */
            border-bottom: 0px;
        }

        QTabBar::tab:selected {
            width: 142px;
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
    
    """
)
        sizePolicy = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sizePolicy)
        self.setMinimumSize(QSize(600, 600))
        self.setMinimumSize(QSize(600, 600))
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

        aboutSection = '''
<p align="justify" style="margin:1em 0; line-height:1.4;">
  <strong>InstaDock‑v2 (IDv2)</strong> is a molecular docking platform designed for efficient, large-scale analysis of protein–ligand interactions. The <strong>Multi Target Docking</strong> module supports simultaneous docking across multiple protein targets, making it ideal for high-throughput screening, target profiling, and comparative binding studies. It leverages optimized docking algorithms for accuracy, supports GPU acceleration to boost performance, and includes advanced visualization tools to simplify the workflow.
</p>

<ul style="margin:0 2em 1em 2em; line-height:1.4;">
  <li>High‑throughput docking across multiple protein targets in a single run</li>
  <li>Automated comparison of ligand binding affinities across diverse targets</li>
  <li>Integrated 2D interaction analysis for visualizing conserved and target‑specific contacts</li>
  <li>Interactive 3D visualization with Mol* for each target–ligand complex</li>
  <li>Optional ADMET prediction add‑on for evaluating pharmacokinetic properties across targets</li>
</ul>

<p align="justify" style="margin:0 0 1em 0; line-height:1.4;">
  After a multi‑target docking run, IDv2 organizes results into a structured output folder that includes:
</p>

<ul style="margin:0 2em 1em 2em; line-height:1.4;">
  <li>Docking scores and ranked poses for each target–ligand combination</li>
  <li>Docked structures and conformers for detailed follow‑up analysis</li>
  <li>Binding affinity distribution charts across all target–ligand sets</li>
  <li>A comparative summary report for cross‑target performance evaluation</li>
</ul>

<p align="center" style="margin-top:1.5em; font-weight:bold; font-size:16px;">
  Discover more and view documentation at<br/>
  <a href="https://hassanlab.in/instadock-v2/" target="_blank">hassanlab.in/instadock-v2</a>
</p>
        '''

        self.tabWidget = QTabWidget(self)
        self.tabWidget.setObjectName("tabWidget")
        sizePolicy1 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        sizePolicy1.setHorizontalStretch(0)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy1.setHeightForWidth(self.tabWidget.sizePolicy().hasHeightForWidth())
        self.tabWidget.setSizePolicy(sizePolicy1)
        self.tabWidget.setUsesScrollButtons(False)

        self.about = QWidget()
        self.about.setObjectName("about")
        self.verticalLayout = QVBoxLayout(self.about)
        self.verticalLayout.setObjectName("verticalLayout")

        self.aboutlabel = QLabel(aboutSection)
        self.aboutlabel.setObjectName("aboutlabel")
        self.aboutlabel.setWordWrap(True)
        self.verticalLayout.addWidget(self.aboutlabel)

        self.tabWidget.addTab(self.about, "About Us")

        usageSection = '''<style> ol {line-height: 1.3;}</style>
        <p>Follow these steps to start docking in InstaDock-v2:</p>
        <p><b>Standard Workflow</b></p>
        <ol>
            <li>Browse for receptor and ligand directories.</li>
            <i>[Formats: PDB for Receptors, SDF & PDBQT for Ligands]</i>
            <li>Select the scoring function.</li>
            <li>Click <b>START</b> to begin docking.</li>
        </ol>

        <p><b>Advanced Workflow</b></p>
        <ol>
            <li>Browse for receptor and ligand directories.</li>
            <i>[Formats: PDB for Receptors, SDF & PDBQT for Ligands]</i>
            <li>Select the scoring function.</li>
            <li>Choose blind or site-specific docking in <b>Dock</b> tab.</li>
            <li>Enable binding site prediction if needed.</li>
            <li>Inspect and prepare receptor file (GPU processing may take time).</li>
            <li>Prepare ligand file(s).</li>
            <li>Generate configuration file.</li>
            <li>Set number of conformers in <b>Settings</b> tab.</li>
            <li>Toggle resume/log options if needed.</li>
            <li>Click <b>START</b> to begin docking.</li>
        </ol>

        <p align="center" style="color:#0c00f0; font-size:20px; font-family:Calibri;">
            <b>Happy Docking!!</b>
        </p>
        '''
        self.usage = QWidget()
        self.usage.setObjectName("usage")
        self.verticalLayout = QVBoxLayout(self.usage)
        self.verticalLayout.setObjectName("verticalLayout")

        self.usagelabel = QLabel(usageSection)
        self.usagelabel.setObjectName("aboutlabel")
        self.usagelabel.setWordWrap(True)
        self.verticalLayout.addWidget(self.usagelabel)

        self.tabWidget.addTab(self.usage, "How to use")

        citeSection = '''
    <div style="position: relative; padding-bottom: 80px;">
        <p>If you use InstaDock-v2 in your work, please cite:</p>
        <p align="justify" style="margin:1em 0; line-height:2;">
        <strong>Placeholder Citation for InstaDock-v2</strong><br/>
        <em>[TBD]. InstaDock-v2: An easy-to-use GPU-Accelerated Platform for Molecular Docking. Journal TBD. Year;Volume(Issue):pages. doi:10.xxxx/instadock-v2.</em>
        </p>

        <p align="justify" style="margin:1em 0; line-height:2;">
        <strong>For the original InstaDock GUI, cite:</strong><br/>
        <em>Mohammad T, Mathur Y, Hassan M. InstaDock: A single-click graphical user interface for molecular docking-based virtual high-throughput screening. Briefings in Bioinformatics. 2021;22(4):bbaa279. doi:10.1093/bib/bbaa279.</em>
        </p>

        <p align="center" style="margin-top:1.5em; font-style:italic;">
        Thank you for using InstaDock-v2! We hope you find it useful.
        </p>

    </div>

        '''
        self.cite = QWidget()
        self.cite.setObjectName("cite")
        self.verticalLayout = QVBoxLayout(self.cite)
        self.verticalLayout.setObjectName("verticalLayout")

        self.citelabel = QLabel(citeSection)
        self.citelabel.setObjectName("aboutlabel")
        self.citelabel.setWordWrap(True)
        self.verticalLayout.addWidget(self.citelabel)

        self.tabWidget.addTab(self.cite, "How to cite")

        contactSection = '''
<div style="position: relative; padding-bottom: 80px;">
  <p align="justify" style="margin:1em 0; line-height:1.5;">
    <strong>Stay Updated!</strong><br/>
    For the latest InstaDock-v2 releases, tutorials, and detailed documentation, visit 
    <a href="https://hassanlab.in/instadock-v2/" target="_blank">hassanlab.in/instadock-v2</a>. 
    To explore legacy materials, archived user guides, and updates for InstaDock‑v1, head over to 
    <a href="https://hassanlab.in/instadock" target="_blank">hassanlab.in/instadock</a>.
  </p>

  <p align="justify" style="margin:1em 0; line-height:1.5;">
    <strong>General Inquiries & Support</strong><br/>
    For questions, comments, collaboration proposals, or anything else, email our support team at 
    <a href="mailto:support@hassanlab.org">support@hassanlab.org</a>. We're here to help!
  </p>

  <p align="center" style="margin-top:1.5em; font-style:italic;">
    Thank you for using InstaDock-v2! <br/>Your feedback and engagement drive our ongoing innovation!
  </p>
</div>
'''

        self.contact = QWidget()
        self.contact.setObjectName("contact")
        self.verticalLayout = QVBoxLayout(self.contact)
        self.verticalLayout.setObjectName("verticalLayout")

        self.contactlabel = QLabel(contactSection)
        self.contactlabel.setObjectName("contactlabel")
        self.contactlabel.setWordWrap(True)
        self.verticalLayout.addWidget(self.contactlabel)

        self.verticalLayout.addStretch()

        h_logo_layout = QHBoxLayout()
        h_logo_layout.addStretch()

        self.logoLabel = QLabel(parent=self.contact)
        self.logoLabel.setObjectName("logoLabel")
        
        pixmap = QPixmap(":/VHL.png")
        scaled = pixmap.scaled(
            200, 200,                      
            Qt.KeepAspectRatio,
            Qt.SmoothTransformation
        )
        self.logoLabel.setPixmap(scaled)
        self.logoLabel.setAlignment(Qt.AlignRight | Qt.AlignBottom)

        h_logo_layout.addWidget(self.logoLabel)
        self.verticalLayout.addLayout(h_logo_layout)

        self.tabWidget.addTab(self.contact, "Contact Us")
        
        work_space_layout.addWidget(self.tabWidget, 0, 0, 1, 1)
        central_widget_layout.addLayout(work_space_layout)
        central_widget.setLayout(central_widget_layout)
        self.setCentralWidget(central_widget)