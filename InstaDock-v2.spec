# -*- mode: python ; coding: utf-8 -*-
import os
from PyInstaller.utils.hooks import collect_data_files

# Get the directory containing this spec file
spec_root = os.path.dirname(os.path.abspath(SPEC))

datas=[
    (os.path.join(spec_root, 'src', 'ID_Predict.pyw'), '.'), 
    (os.path.join(spec_root, 'src', 'ID_2D.pyw'), '.'), 
    (os.path.join(spec_root, 'src', 'ID_About_Multi.pyw'), '.'), 
    (os.path.join(spec_root, 'src', 'ID_About_Single.pyw'), '.'), 
    (os.path.join(spec_root, 'src', 'ID_ADMET.pyw'), '.'), 
    (os.path.join(spec_root, 'src', 'ID_Start_Calc_Multi.pyw'), '.'), 
    (os.path.join(spec_root, 'src', 'ID_Start_Calc_Single.pyw'), '.'), 
    (os.path.join(spec_root, 'src', 'ID_Ligand_Downloader.pyw'), '.'), 
    (os.path.join(spec_root, 'src', 'ID_Lipinski.pyw'), '.'), 
    (os.path.join(spec_root, 'src', 'ID_Prep_Multi.pyw'), '.'), 
    (os.path.join(spec_root, 'src', 'ID_Prep_Single.pyw'), '.'), 
    (os.path.join(spec_root, 'src', 'ID_utils.pyw'), '.'), 
    (os.path.join(spec_root, 'src', 'IDv2Resource.py'), '.'), 
    (os.path.join(spec_root, 'src', 'ID_Visualizer.pyw'), '.'),  
    (os.path.join(spec_root, 'src', 'ID_Multi_Main.pyw'), '.'),
    (os.path.join(spec_root, 'src', 'ID_Single_Main.pyw'), '.'), 
    (os.path.join(spec_root, 'src', 'template.html'), '.'),  
    (os.path.join(spec_root, 'src', 'idv2_hashes.txt'), '.'),
    (os.path.join(spec_root, 'src', 'chart.js'), '.')
]

# Collect data files from packages
datas += collect_data_files('openbabel')
datas += collect_data_files('pyKVFinder')
datas += collect_data_files('pdbfixer')
datas += collect_data_files('prolif')

a = Analysis(
    [os.path.join(spec_root, 'src', 'ID_Main.py')],
    pathex=[spec_root],
    binaries=[
        (os.path.join(spec_root, 'binaries', 'AutoDock-GPU.exe'), '.'), 
        (os.path.join(spec_root, 'binaries', 'autogrid4.exe'), '.'), 
        (os.path.join(spec_root, 'binaries', 'cyggcc_s-seh-1.dll'), '.'), 
        (os.path.join(spec_root, 'binaries', 'cyggomp-1.dll'), '.'), 
        (os.path.join(spec_root, 'binaries', 'cygstdc++-6.dll'), '.'), 
        (os.path.join(spec_root, 'binaries', 'cygwin1.dll'), '.'), 
        (os.path.join(spec_root, 'binaries', 'libomp140.x86_64.dll'), '.'),
        (os.path.join(spec_root, 'binaries', 'mk_prepare_receptor.exe'), '.'), 
        (os.path.join(spec_root, 'binaries', 'qvinaw.exe'), '.'), 
        (os.path.join(spec_root, 'binaries', 'vina.exe'), '.'),
    ],
    datas=datas,
    hiddenimports=['MDAnalysis.lib.formats.cython_util', 'numpy', 'scipy', 'pandas', 'torch', 'lightning', 'rdkit',
        'multiprocessing.Pool', 'multiprocessing.RLock', 'catboost', 'numpy.core.isnan', 'numpy.core.count_nonzero',],
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
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=[os.path.join(spec_root, 'assets', 'IDLogo.ico')],
)