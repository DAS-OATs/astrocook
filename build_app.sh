#!/bin/bash
rm -rf build dist *.spec  # Pulisce le vecchie build
pyinstaller --clean \
            --noconfirm \
            --name="Astrocook" \
            --windowed \
            --icon="logo/icon_3d.icns" \
            --paths="." \
            --collect-all asdf \
            --collect-all scienceplots \
            --add-data "$(pwd)/astrocook/data:astrocook/data" \
            --hidden-import="astrocook.v1.vars" \
            tests/launch_pyside_app.py