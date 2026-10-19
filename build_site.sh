#!/bin/bash
# Builds the documentation website into the directory given as argument:
# V2 docs at the root, the frozen V1 docs (docs_backup/) under /v1/.
# Used by deploy_docs.sh (by hand) and by .github/workflows/docs.yml.
set -e

OUT_DIR="$1"
if [ -z "$OUT_DIR" ]; then
    echo "Usage: ./build_site.sh <output directory>"
    exit 1
fi
if ! command -v sphinx-build &> /dev/null; then
    echo "Error: sphinx-build not found. Install it with: pip install -e \".[docs]\""
    exit 1
fi
mkdir -p "$OUT_DIR"

# 1. API stubs and HTML
./generate_sphinx.sh
make -C docs clean html

# 2. Assemble the site
cp -r docs/_build/html/* "$OUT_DIR"
cp switcher.json "$OUT_DIR/"
mkdir -p "$OUT_DIR/v1"
cp -r docs_backup/* "$OUT_DIR/v1/"
# The Jekyll config at the root makes GitHub Pages load the V1 theme
cp docs_backup/_config.yml "$OUT_DIR/"

# 3. GitHub Pages runs Jekyll, which ignores folders starting with "_":
# rename them and patch the references (perl -pi works on macOS and Linux)
cd "$OUT_DIR"
for d in _static _images _sources; do
    [ -d "$d" ] && mv "$d" "${d#_}"
done
find . \( -name '*.html' -o -name '*.js' \) -type f \
    -exec perl -pi -e 's/_static/static/g; s/_images/images/g; s/_sources/sources/g' {} +
