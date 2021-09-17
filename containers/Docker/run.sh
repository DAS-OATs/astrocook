#!/bin/bash
docker run -v$PWD/:/data -p5900:5900 -p8590:8590 -it astrocook
