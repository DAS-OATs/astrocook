#!/bin/bash

echo "You need XQuartz and to tick XQuartz -> Preferences -> Security -> Allow connections from network clients"

xhost + 127.0.0.1

docker run -v$PWD/:/data -v$PWD/../../:/metauser/astrocook -e DISPLAY=docker.for.mac.host.internal:0 -it astrocook python3 -u /metauser/astrocook/ac_gui.py
