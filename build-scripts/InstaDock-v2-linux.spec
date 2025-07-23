# -*- mode: python ; coding: utf-8 -*-
import os
from PyInstaller.utils.hooks import collect_data_files

# Get repository root directory
repo_root = os.path.dirname(os.path.dirname(os.path.abspath(SPEC)))

# Collect pyKVFinder data
testData = collect_data_files('pyKVFinder', subdir='data')

# Platform-specific data files (Linux uses .py extensions)
linux_datas = [
    (os.path.join(repo_root, 'src', 'ID_Predict.py'), '.'),
    (os.path.join(repo_root, 'src', 'ID_2D.py'), '.'),
    (os.path.join(repo_root, 'src', 'ID_About_Multi.py'), '.'),
    (os.path.join(repo_root, 'src', 'ID_About_Single.py'), '.'),
    (os.path.join(repo_root, 'src', 'ID_ADMET.py'), '.'),
    (os.path.join(repo_root, 'src', 'ID_Start_Calc_Multi.py'), '.'),
    (os.path.join(repo_root, 'src', 'ID_Start_Calc_Single.py'), '.'),
    (os.path.join(repo_root, 'src', 'ID_Ligand_Downloader.py'), '.'),
    (os.path.join(repo_root, 'src', 'ID_Lipinski.py'), '.'),
    (os.path.join(repo_root, 'src', 'ID_Prep_Multi.py'), '.'),
    (os.path.join(repo_root, 'src', 'ID_Prep_Single.py'), '.'),
    (os.path.join(repo_root, 'src', 'ID_utils.py'), '.'),
    (os.path.join(repo_root, 'src', 'IDv2Resource.py'), '.'),
    (os.path.join(repo_root, 'src', 'ID_Visualizer.py'), '.'),
    (os.path.join(repo_root, 'src', 'ID_Multi_Main.py'), '.'),
    (os.path.join(repo_root, 'src', 'ID_Single_Main.py'), '.'),
    (os.path.join(repo_root, 'src', 'template.html'), '.'),
    (os.path.join(repo_root, 'src', 'idv2_hashes.txt'), '.'),
]

# Linux-specific binaries and system libraries
linux_binaries = [
    (os.path.join(repo_root, 'binaries', 'linux', 'obabel'), '.'),
    (os.path.join(repo_root, 'binaries', 'linux', 'adgpu-v1.6_linux_x64_ocl_128wi'), '.'),
    (os.path.join(repo_root, 'binaries', 'linux', 'autogrid4'), '.'),
    (os.path.join(repo_root, 'binaries', 'linux', 'mk_prepare_receptor'), '.'),
    (os.path.join(repo_root, 'binaries', 'linux', 'qvinaw'), '.'),
    (os.path.join(repo_root, 'binaries', 'linux', 'vina'), '.'),
]

# Add system libraries if they exist (with fallback handling)
system_libs = [
    "/usr/lib/x86_64-linux-gnu/librsvg-2.so.2",
    "/usr/lib/x86_64-linux-gnu/libpango-1.0.so.0"
]

for lib in system_libs:
    if os.path.exists(lib):
        linux_binaries.append((lib, '.'))

a = Analysis(
    [os.path.join(repo_root, 'src', 'ID_Main.py')],
    pathex=[repo_root],
    binaries=linux_binaries,
    datas=testData + linux_datas,
    hiddenimports=['MDAnalysis.lib.formats.cython_util', 'numpy', 'scipy', 'pandas',
                   'torch', 'lightning', 'multiprocessing.Pool', 'multiprocessing.RLock',
                   'catboost', 'numpy.core.isnan', 'numpy.core.count_nonzero'],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[os.path.join(repo_root, 'build-scripts', 'hook-pixbuf.py'),
                   os.path.join(repo_root, 'build-scripts', 'hook-libraries.py')],
    excludes=['tkinter'],
    noarchive=True,
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
    icon=[os.path.join(repo_root, 'assets', 'IDLogo.png')],
)
