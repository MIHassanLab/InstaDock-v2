# InstaDock-v2: Automated Molecular Docking with Integrated ADMET Prediction

[![License](https://img.shields.io/badge/License-Academic%20Research-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)
[![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20Linux-lightgrey.svg)]()

An automated molecular docking tool with integrated ADMET prediction and analysis for Windows and Linux platforms.

---

## 🚀 Quick Start - Download Pre-built Executables

For easy installation without Python setup, download ready-to-use executables:

**[📥 Download InstaDock-v2 from hassanlab.in](https://hassanlab.in/instadock-v2)**

*Supported platforms: Windows 10/11, Linux (Ubuntu 20.04+)*

---

## 📑 Table of Contents
- [Features](#-features)
- [System Requirements](#️-system-requirements)
- [Installation Instructions](#-installation-instructions)
- [Quick Demo](#-quick-demo)
- [Academic Use and Citation](#-academic-use-and-citation)
- [License](#-license)
- [Authors & Support](#-authors--support)
- [Acknowledgments](#-acknowledgments)

## 📋 Features

- **Automated molecular docking** with multiple backend engines (AutoDock Vina, Quickvina-W, AutoDock-GPU)
- **High-throughput screening** capabilities for large compound libraries
- **Integrated ADMET prediction** for drug-like properties assessment
- **Cross-platform support** (Windows/Linux) with intuitive PyQt5 interface
- **2D molecular visualization** and interactive results analysis
- **Batch processing** capabilities for multiple ligands/targets

---

## ⚙️ System Requirements

### Operating Systems
- **Windows**: 10, 11 (64-bit)
- **Linux**: Ubuntu 20.04+, or equivalent distributions

### Hardware Requirements
- **Memory**: 8GB RAM minimum, 16GB recommended
- **Storage**: 2GB free disk space
- **CPU**: Multi-core processor recommended for optimal performance

### Software Dependencies

*****NONE*** if using pre-build executables**

**Essential Requirements:**
- **Python 3.10+** (for source installation only)
- **PyQt5** (GUI framework)

**External Binary Dependencies:**
The following third-party software must be installed separately by users if compiling from source:

| Software | Purpose | Download Link | Notes |
|----------|---------|---------------|-------|
| **QuickVina-W** | Primary docking engine | [github.com/Qvina/qvina-w](https://github.com/QVina/qvina/blob/master/bin/qvina-w) | Essential for all docking |
| **AutoDock Vina** | Additional docking engine | [vina.scripps.edu](http://vina.scripps.edu/download.html) | Essential for all docking |
| **AutoDock-GPU** | GPU-accelerated docking | [github.com/ccsb-scripps/AutoDock-GPU](https://github.com/ccsb-scripps/AutoDock-GPU) | Optional but recommended |
| **MGLTools** | Molecule preparation | [mgltools.scripps.edu](http://mgltools.scripps.edu/downloads) | Required for advanced features |
| **AutoGrid** | Molecule preparation | [github.com/ccsb-scripps/AutoGrid](https://github.com/ccsb-scripps/AutoGrid) | Required for advanced features |
| **Meeko** | Molecule preparation | [github.com/forlilab/Meeko](https://github.com/forlilab/Meeko) | Required for advanced features |
| **OpenBabel** | File format conversion | [openbabel.org](https://openbabel.org/wiki/Category:Installation) | Required for ADMET prediction |

---

## 🔧 Installation Instructions

### Option 1: Pre-built Executables (Recommended)

1. **Download InstaDock-v2 executable** from [hassanlab.in/instadock-v2](https://hassanlab.in/instadock-v2)
2. **Run InstaDock-v2** executable

**Installation time**: ~10-15 minutes (including dependencies)

### Option 2: From Source Code

#### Clone and Setup
```bash
git clone https://github.com/yourusername/InstaDock-v2.git
cd InstaDock-v2
pip install -r requirements.txt
```

#### Install External Binary Dependencies

**Important**: InstaDock-v2 requires several third-party molecular docking tools to function properly. These must be installed separately:

**For Windows Users:**
1. Download and install **AutoDock Vina** from the official website
2. Download and install **QuickVina-W** from the official website
3. Download and install **AutoDock GPU** from the official website
4. Download and install **Meeko** from the official website
5. Download and install **AutoGrid** from the official website
6. Download **MGLTools** and follow installation instructions
7. Install **OpenBabel** for Windows
8. Ensure all executables are accessible in your system PATH

**For Linux Users:**
```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install autodock-vina openbabel

# CentOS/RHEL
sudo yum install autodock-vina openbabel 

# Download and install QuickVina-w manually from the official website
# Download and install Autodock-GPU manually from the official website
# Download and install Meeko manually from the official website
# Download and install AutoGrid manually from the official website

```

#### Linux Compatibility Quick-Fix

1. Replace the resource_path helper with an OS-agnostic version:

```bash
from pathlib import Path
import os, sys

def resource_path(relative_path: str) -> str:
    """
    Return absolute path to <relative_path> regardless of
    ① PyInstaller one-file bundles (sys._MEIPASS)
    ② Plain source runs (script folder)
    Works unchanged on Windows, Linux, macOS.
    """
    base = getattr(sys, "_MEIPASS", Path(__file__).resolve().parent)
    return os.path.join(base, relative_path)
```

2. Add a binary locator that accommodates AppImage and system PATH:
```bash

python
import os, sys

def find_binary(name: str) -> str:
    """
    Locate an executable called <name>.

    Priority:
     1. PyInstaller bundle (sys._MEIPASS)
     2. AppImage mount point ($APPDIR/usr/bin)
     3. System PATH
    """
    # ① PyInstaller one-file
    if getattr(sys, "frozen", False):
        candidate = os.path.join(sys._MEIPASS, name)
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate

    # ② Inside AppImage
    appdir = os.environ.get("APPDIR")
    if appdir:
        candidate = os.path.join(appdir, "usr", "bin", name)
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate

    # ③ Fallback
    return name
```
3. Adjust path separators and executable names
    - Replace hard-coded \ with os.path.join() or forward slashes.
    - Call find_binary("vina") instead of "vina.exe".

## ADMET Model Configuration

InstaDock-v2 requires ADMET model weights packaged as a ZIP.  
1. Request your secure download URL by emailing support@hassanlab.org.  
2. Copy the URL into `models_config.json`.  
3. Run InstaDock-v2; the tool will fetch and unpack into automatically.

#### Run from Source
```bash
python src/ID_Main.py
```

---

## 📖 Binary Dependencies Setup Guide

### Why External Dependencies?

InstaDock-v2 integrates multiple specialized molecular docking engines and cheminformatics tools. Rather than redistributing these third-party binaries (which may have licensing restrictions), we require users to install them directly from official sources to ensure if compiling from source:

- **License compliance** with original software terms
- **Latest versions** with security updates and bug fixes
- **Proper attribution** to original developers
- **Reduced repository size** for faster downloads

### Installation Verification

After installing external dependencies, verify they work correctly:

```bash
# Test AutoDock Vina
vina --help

# Test OpenBabel
obabel -H

# Test MGLTools (if installed)
pythonsh -c "import AutoDockTools"
```

### Troubleshooting Common Issues

**Issue**: "Command not found" errors
**Solution**: Add binary locations to your system PATH environment variable

**Issue**: Permission denied on Linux
**Solution**: Ensure executables have proper permissions:
```bash
chmod +x /path/to/vina
chmod +x /path/to/autodock_gpu
chmod +x /path/to/openbabel
```

**Issue**: Library conflicts
**Solution**: Use virtual environments to isolate dependencies

**Issue**: InstaDock-v2 executable does not start on Linux
**Solution**: Ensure proper permission for the executable:
```bash
chmod +x /path/to/InstaDock-v2-x86_64.AppImage
```

---

## 🚦 Quick Demo

1. **Launch InstaDock-v2** (executable or `python src/ID_Main.py`)
2. **Load sample data** from the `data/sample_data/` directory
3. **Configure docking parameters** or use defaults
4. **Start docking calculation** 
5. **Analyze results** in the visualization panel

**Demo runtime**: ~5-10 minutes on standard desktop hardware

---

## 🤝 Academic Use and Citation

If you use InstaDock-v2 in your research, please cite:

```bibtex
@article{mathur2025instadock,
  title={InstaDock-v2: an automated molecular docking tool with integrated ADMET prediction and analysis},
  author={Mathur, Yash and Hassan, Md. Imtaiyaz},
  journal={TBD},
  year={TBD},
  note={TBD}
}
```

---

## 📄 License

This software is available for **academic and non-commercial research use only** under the Academic Research License. See [LICENSE](LICENSE) for full terms.

**Key Points:**
- ✅ Free for academic research and education
- ✅ Source code available for transparency and reproducibility
- ✅ Attribution required in publications

---

## 👥 Authors & Support

**Lead Developer**: Yash Mathur  
**Principal Investigator**: Md. Imtaiyaz Hassan  
**Affiliation**: Center for Interdisciplinary Research in Basic Sciences, Jamia Millia Islamia, New Delhi, India

### Getting Help

- **🐛 Bug Reports**: [GitHub Issues](https://github.com/MIHassanLab/InstaDock-v2/issues)
- **💬 Academic Support**: Contact [support@hassanlab.org](mailto:support@hassanlab.org)
- **🤝 Collaborations**: We welcome research partnerships and feedback

---

## 🙏 Acknowledgments

We thank the developers of AutoDock Vina, OpenBabel, PyQt5, and the broader computational biology community for providing the foundational tools that make InstaDock-v2 possible.

---

**Ready to start molecular docking?** [Download InstaDock-v2](https://hassanlab.in/instadock-v2) and begin your computational drug discovery journey!
