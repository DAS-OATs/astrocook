#!/bin/bash
if ! command -v sphinx-build &> /dev/null; then
    echo "Error: sphinx-build not found. Install it with: pip install -e \".[docs]\""
    exit 1
fi

cd docs
make clean
make html
cd ..