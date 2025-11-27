#!/bin/bash
rm -rf build dist release
# Nota: NON cancellare più *.spec!

# Lancia PyInstaller usando il file spec
pyinstaller --clean --noconfirm Astrocook.spec

#!/bin/bash

# 1. Pulizia
echo "🧹  Clean previous builds..."
rm -rf build dist 
# Nota: Rimuoviamo anche eventuali DMG vecchi per evitare errori
rm -f Astrocook.dmg

# 2. Creazione dell'App (PyInstaller)
echo "🚀  Building App..."
pyinstaller --clean --noconfirm Astrocook.spec

# Controllo se la compilazione è andata a buon fine
if [ ! -d "dist/Astrocook.app" ]; then
    echo "❌  Error: App creation failed."
    exit 1
fi

# 3. Creazione del DMG
echo "📦  Creating DMG package..."

# Verifica se create-dmg è installato
if ! command -v create-dmg &> /dev/null; then
    echo "⚠️  create-dmg not found. Install it with 'brew install create-dmg'"
    echo "The App was created in dist/Astrocook.app, but the DMG was not generated."
    exit 0
fi

# Creiamo una cartella per l'output finale se non esiste
mkdir -p release

# Comando magico per creare il DMG
create-dmg \
  --volname "Astrocook Installer" \
  --volicon "assets/icon_disk_3d_HR.icns" \
  --window-pos 200 120 \
  --window-size 500 320 \
  --icon-size 100 \
  --icon "Astrocook.app" 100 150 \
  --background "assets/logo_background_flat_HR.png" \
  --hide-extension "Astrocook.app" \
  --app-drop-link 380 150 \
  "release/Astrocook.dmg" \
  "dist/Astrocook.app"

echo "✅  Done! Your DMG is ready in: release/Astrocook.dmg"