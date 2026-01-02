#!/bin/bash

# --- Configuration ---
DOCS_DIR="docs"
BUILD_OUTPUT="docs/_build/html"
TARGET_BRANCH="gh-pages"
TARGET_FOLDER="dev"
SOURCE_PACKAGE="astrocook" 

# Stop the script if any command fails
set -e

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}=== Starting Documentation Deployment (Jekyll-Friendly Mode) ===${NC}"

# 0. Safety Check
if [[ -n $(git status -s) ]]; then
    echo "Error: Your working directory is not clean."
    echo "Please commit or stash your changes before running this script."
    exit 1
fi

CURRENT_BRANCH=$(git branch --show-current)

# 1. Regenerate API stubs
echo -e "${BLUE}Step 1: Regenerating API stubs...${NC}"
sphinx-apidoc -f -e --no-toc -o "$DOCS_DIR/api" "$SOURCE_PACKAGE" \
    astrocook/core/atomic_data.py \
    astrocook/core/photometry.py \
    astrocook/core/session_manager.py \
    astrocook/core/spectrum_operations.py \
    astrocook/core/system_list_migration.py \
    astrocook/core/utils.py \
    astrocook/gui/*.py \
    astrocook/io/*.py

rm docs/api/astrocook.rst docs/api/astrocook.settings.rst docs/api/modules.rst

# 2. Build the HTML
echo -e "${BLUE}Step 2: Building HTML...${NC}"
cd "$DOCS_DIR"
make clean
make html
cd ..

# 3. Post-Processing: Rename folders to bypass Jekyll ignore rules
# This is crucial for hosting on gh-pages without a .nojekyll file
echo -e "${BLUE}Step 3: Patching HTML for Jekyll compatibility...${NC}"

# Create a temp dir
TEMP_DIR=$(mktemp -d)
cp -r "$BUILD_OUTPUT"/* "$TEMP_DIR"

# Go into the temp dir to perform operations
pushd "$TEMP_DIR" > /dev/null

# A. Rename directories if they exist
[ -d "_static" ] && mv "_static" "static"
[ -d "_images" ] && mv "_images" "images"
[ -d "_sources" ] && mv "_sources" "sources"

# B. Patch the HTML files (MacOS compatible sed)
# We look for references to _static, _images, etc., and remove the underscore
find . -name '*.html' -type f -exec sed -i '' 's/_static/static/g' {} +
find . -name '*.html' -type f -exec sed -i '' 's/_images/images/g' {} +
find . -name '*.html' -type f -exec sed -i '' 's/_sources/sources/g' {} +

# C. Patch the Javascript search index (attempts to fix search box)
find . -name '*.js' -type f -exec sed -i '' 's/_static/static/g' {} +
find . -name '*.js' -type f -exec sed -i '' 's/_images/images/g' {} +
find . -name '*.js' -type f -exec sed -i '' 's/_sources/sources/g' {} +

popd > /dev/null

# 4. Switch to gh-pages and Deploy
echo -e "${BLUE}Step 4: Switching to $TARGET_BRANCH...${NC}"
git checkout "$TARGET_BRANCH"

echo -e "${BLUE}Step 5: Updating '$TARGET_FOLDER' directory...${NC}"
mkdir -p "$TARGET_FOLDER"
rm -rf "$TARGET_FOLDER"/*
cp -r "$TEMP_DIR"/* "$TARGET_FOLDER"/

echo -e "${BLUE}Step 6: Committing and Pushing...${NC}"
git add "$TARGET_FOLDER"
git commit -m "Deploy Jekyll-patched docs to dev folder" || echo "No changes to commit"
git push origin "$TARGET_BRANCH"

# 5. Cleanup
echo -e "${BLUE}Step 7: Returning to $CURRENT_BRANCH...${NC}"
rm -rf "$TEMP_DIR"
git checkout "$CURRENT_BRANCH"

echo -e "${GREEN}=== Deployment Complete! ===${NC}"
echo "Your patched dev docs are at: https://das-oats.github.io/astrocook/$TARGET_FOLDER/"