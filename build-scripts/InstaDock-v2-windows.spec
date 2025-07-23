# -*- mode: python ; coding: utf-8 -*-
import os
from PyInstaller.utils.hooks import collect_data_files

# Get repository root directory
repo_root = os.path.dirname(os.path.dirname(os.path.abspath(SPEC)))

# Platform-specific data files (Windows uses .pyw extensions)
datas = [
    (os.path.join(repo_root, 'src', 'ID_Predict.pyw'), '.'), 
    (os.path.join(repo_root, 'src', 'ID_2D.pyw'), '.'), 
    (os.path.join(repo_root, 'src', 'ID_About_Multi.pyw'), '.'), 
    (os.path.join(repo_root, 'src', 'ID_About_Single.pyw'), '.'), 
    (os.path.join(repo_root, 'src', 'ID_ADMET.pyw'), '.'), 
    (os.path.join(repo_root, 'src', 'ID_Start_Calc_Multi.pyw'), '.'), 
    (os.path.join(repo_root, 'src', 'ID_Start_Calc_Single.pyw'), '.'), 
    (os.path.join(repo_root, 'src', 'ID_Ligand_Downloader.pyw'), '.'), 
    (os.path.join(repo_root, 'src', 'ID_Lipinski.pyw'), '.'), 
    (os.path.join(repo_root, 'src', 'ID_Prep_Multi.pyw'), '.'), 
    (os.path.join(repo_root, 'src', 'ID_Prep_Single.pyw'), '.'), 
    (os.path.join(repo_root, 'src', 'ID_utils.pyw'), '.'), 
    (os.path.join(repo_root, 'src', 'IDv2Resource.py'), '.'), 
    (os.path.join(repo_root, 'src', 'ID_Visualizer.pyw'), '.'),  
    (os.path.join(repo_root, 'src', 'ID_Multi_Main.pyw'), '.'),
    (os.path.join(repo_root, 'src', 'ID_Single_Main.pyw'), '.'), 
    (os.path.join(repo_root, 'src', 'template.html'), '.'),  
    (os.path.join(repo_root, 'src', 'idv2_hashes.txt'), '.'),
    (os.path.join(repo_root, 'src', 'chart.js'), '.')
]

# Collect package data files
datas += collect_data_files('openbabel')
datas += collect_data_files('pyKVFinder')
datas += collect_data_files('pdbfixer')
datas += collect_data_files('prolif')

# Windows-specific binaries
windows_binaries = [
    (os.path.join(repo_root, 'binaries', 'windows', 'AutoDock-GPU.exe'), '.'), 
    (os.path.join(repo_root, 'binaries', 'windows', 'autogrid4.exe'), '.'), 
    (os.path.join(repo_root, 'binaries', 'windows', 'cyggcc_s-seh-1.dll'), '.'), 
    (os.path.join(repo_root, 'binaries', 'windows', 'cyggomp-1.dll'), '.'), 
    (os.path.join(repo_root, 'binaries', 'windows', 'cygstdc++-6.dll'), '.'), 
    (os.path.join(repo_root, 'binaries', 'windows', 'cygwin1.dll'), '.'), 
    (os.path.join(repo_root, 'binaries', 'windows', 'libomp140.x86_64.dll'), '.'),
    (os.path.join(repo_root, 'binaries', 'windows', 'mk_prepare_receptor.exe'), '.'), 
    (os.path.join(repo_root, 'binaries', 'windows', 'qvinaw.exe'), '.'), 
    (os.path.join(repo_root, 'binaries', 'windows', 'vina.exe'), '.'),
]

a = Analysis(
    [os.path.join(repo_root, 'src', 'ID_Main.py')],
    pathex=[repo_root],
    binaries=windows_binaries,
    datas=datas,
    hiddenimports=['MDAnalysis.lib.formats.cython_util', 'numpy', 'scipy', 'pandas', 
                   'torch', 'lightning', 'rdkit', 'multiprocessing.Pool', 
                   'multiprocessing.RLock', 'catboost', 'numpy.core.isnan', 
                   'numpy.core.count_nonzero'],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='InstaDock-v2',
    debug=False,
    bootloader_ignore_signals=False,
    strip=True,
    upx=False,
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=[os.path.join(repo_root, 'assets', 'IDLogo.ico')],
)
