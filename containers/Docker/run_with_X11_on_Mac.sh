#!/bin/bash

echo "Note: you need XQuartz and to tick XQuartz -> Preferences -> Security -> Allow connections from network clients"

xhost + 127.0.0.1

docker run -it -e DISPLAY=docker.for.mac.host.internal:0 sarusso/astrocook:9ddab14 python3 /metauser/astrocook/ac_gui.py


# p.s. https://github.com/docker/for-mac/issues/2965
