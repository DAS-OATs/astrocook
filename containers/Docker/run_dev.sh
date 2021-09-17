#!/bin/bash

if [ ! -f "Dockerfile" ]; then
    # TODO: This check is weak: improve me!
    echo "Please run this script from the containers/Docker folder"
    exit 1
fi

docker run -v$PWD/:/data -v$PWD/../../:/metauser/astrocook  -p5900:5900 -p8590:8590 -it astrocook
