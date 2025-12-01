# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import collect_all

datas = [('/Users/guido/Library/CloudStorage/GoogleDrive-guido.cupani@inaf.it/My Drive/GitHub/astrocook/astrocook/data', 'astrocook/data')]
binaries = []
hiddenimports = ['astrocook.legacy.vars']
tmp_ret = collect_all('asdf')
datas += tmp_ret[0]; binaries += tmp_ret[1]; hiddenimports += tmp_ret[2]
tmp_ret = collect_all('scienceplots')
datas += tmp_ret[0]; binaries += tmp_ret[1]; hiddenimports += tmp_ret[2]


a = Analysis(
    ['tests/launch_pyside_app.py'],
    pathex=['.'],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='Astrocook',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=['assets/icon_3d_HR.icns'],
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='Astrocook',
)
app = BUNDLE(
    coll,
    name='Astrocook.app',
    icon='assets/icon_3d_HR.icns',
    bundle_identifier='com.inaf.astrocook',
    info_plist={
        'CFBundleDocumentTypes': [
            {
                'CFBundleTypeName': 'Astrocook Session',
                'CFBundleTypeIconFile': 'assets/icon_file_3d_HR.icns', # Usa la stessa icona dell'app per i file
                'LSItemContentTypes': ['com.astrocook.acs2'],
                'CFBundleTypeRole': 'Editor',
                'LSHandlerRank': 'Owner',
                'CFBundleTypeExtensions': ['acs2']
            },
            {
                # Opzionale: Se vuoi associare anche i file .acs vecchi
                'CFBundleTypeName': 'Astrocook Session V1',
                'CFBundleTypeIconFile': 'assets/icon_file_3d_HR.icns',
                'LSItemContentTypes': ['com.astrocook.acs'],
                'CFBundleTypeRole': 'Editor',
                'LSHandlerRank': 'Owner',
                'CFBundleTypeExtensions': ['acs']
            }
        ],
        'UTExportedTypeDeclarations': [
            {
                'UTTypeIdentifier': 'com.astrocook.acs2',
                'UTTypeConformsTo': ['public.data'],
                'UTTypeDescription': 'Astrocook Session',
                'UTTypeIconFile': 'assets/icon_file_3d_HR.icns',
                'UTTypeTagSpecification': {
                    'public.filename-extension': ['acs2']
                }
            },
            {
                'UTTypeIdentifier': 'com.astrocook.acs',
                'UTTypeConformsTo': ['public.data'],
                'UTTypeDescription': 'Astrocook Session V1',
                'UTTypeIconFile': 'assets/icon_file_3d_HR.icns',
                'UTTypeTagSpecification': {
                    'public.filename-extension': ['acs']
                }
            }
        ],
        # Opzionale: Forza la modalità scura/chiara se serve
        'NSHighResolutionCapable': 'True'
    }
)
