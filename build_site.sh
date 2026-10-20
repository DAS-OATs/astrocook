#!/bin/bash
# Builds the documentation website into the directory given as argument.
#   ./build_site.sh <output directory> [stable|dev]
# stable (default): docs of the checked-out release for the root of the
#   site, with the frozen V1 docs (docs_backup/) under /v1/;
# dev: docs of master only, published under /dev/.
# Used by deploy_docs.sh (by hand) and by .github/workflows/docs.yml.
set -e

OUT_DIR="$1"
TARGET="${2:-stable}"
if [ -z "$OUT_DIR" ] || { [ "$TARGET" != stable ] && [ "$TARGET" != dev ]; }; then
    echo "Usage: ./build_site.sh <output directory> [stable|dev]"
    exit 1
fi
if ! command -v sphinx-build &> /dev/null; then
    echo "Error: sphinx-build not found. Install it with: pip install -e \".[docs]\""
    exit 1
fi
mkdir -p "$OUT_DIR"

# 1. API stubs and HTML; the version switcher marks the dev docs as "dev"
# (docs/conf.py), the release docs with the package version
./generate_sphinx.sh
if [ "$TARGET" = dev ]; then
    export ASTROCOOK_DOCS_VERSION=dev
fi
make -C docs clean html

# 2. Assemble the site
cp -r docs/_build/html/* "$OUT_DIR"
if [ "$TARGET" = stable ]; then
    cp switcher.json "$OUT_DIR/"
    mkdir -p "$OUT_DIR/v1"
    cp -r docs_backup/* "$OUT_DIR/v1/"
    # The Jekyll config at the root makes GitHub Pages load the V1 theme
    cp docs_backup/_config.yml "$OUT_DIR/"
fi

# 3. GitHub Pages runs Jekyll, which ignores folders starting with "_":
# rename them and patch the references (perl -pi works on macOS and Linux).
# _modules holds the highlighted source code behind the "[source]" links.
cd "$OUT_DIR"
for d in _static _images _sources _modules; do
    if [ -d "$d" ]; then mv "$d" "${d#_}"; fi
done
find . \( -name '*.html' -o -name '*.js' \) -type f \
    -exec perl -pi -e 's/_static/static/g; s/_images/images/g; s/_sources/sources/g; s{(?<!\w)_modules/}{modules/}g' {} +
