#!/bin/bash
docker run -v$PWD/:/data -e DISPLAY=docker.for.mac.host.internal:0 -it astrocook python3 -u /metauser/astrocook/ac_gui.py
