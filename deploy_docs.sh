#!/bin/bash

# Builds the website (build_site.sh) and pushes it to gh-pages by hand.
# The same happens automatically after each merge to master
# (.github/workflows/docs.yml); this script stays as a fallback.

# --- Configuration ---
TARGET_BRANCH="gh-pages"

# Stop the script if any command fails
set -e

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}=== Starting Documentation Deployment (Jekyll-Friendly Mode) ===${NC}"

# 0. Safety Check
if ! command -v sphinx-apidoc &> /dev/null; then
    echo "Error: sphinx-apidoc not found. Install it with: pip install -e \".[docs]\""
    exit 1
fi
if [[ -n $(git status -s) ]]; then
    echo "Error: Your working directory is not clean."
    echo "Please commit or stash your changes before running this script."
    exit 1
fi

CURRENT_BRANCH=$(git branch --show-current)

# 1-3. Build the site (V2 docs, V1 docs under /v1/, Jekyll patches)
echo -e "${BLUE}Steps 1-3: Building the site...${NC}"
TEMP_DIR=$(mktemp -d)
./build_site.sh "$TEMP_DIR"

# 4. Switch to gh-pages and Deploy
echo -e "${BLUE}Step 4: Switching to $TARGET_BRANCH...${NC}"
git checkout "$TARGET_BRANCH"

echo -e "${BLUE}Step 5: Updating root directory...${NC}"
# Copy the new V2 files to the root, overwriting the old ones
cp -r "$TEMP_DIR"/* .

echo -e "${BLUE}Step 6: Committing and Pushing...${NC}"
git add .
git commit -m "Deploy V2 docs to root" || echo "No changes to commit"
git push origin "$TARGET_BRANCH"

# 5. Cleanup
echo -e "${BLUE}Step 7: Returning to $CURRENT_BRANCH...${NC}"
rm -rf "$TEMP_DIR"
git checkout "$CURRENT_BRANCH"

echo -e "${GREEN}=== Deployment Complete! ===${NC}"
echo "Your patched V2 docs are at: https://das-oats.github.io/astrocook/"