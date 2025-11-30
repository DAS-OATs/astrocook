#!/bin/bash

# --- Configuration ---
DOCS_DIR="docs"
BUILD_OUTPUT="docs/_build/html"
TARGET_BRANCH="gh-pages"
TARGET_FOLDER="dev"
SOURCE_PACKAGE="astrocook" # Change this if your source is just "astrocook"

# Stop the script if any command fails
set -e

# Colors for pretty output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}=== Starting Documentation Deployment ===${NC}"

# 0. Safety Check: Ensure git is clean
if [[ -n $(git status -s) ]]; then
    echo "Error: Your working directory is not clean."
    echo "Please commit or stash your changes before running this script."
    exit 1
fi

# 1. Get the current branch name so we can switch back later
CURRENT_BRANCH=$(git branch --show-current)
echo -e "Current branch is: ${GREEN}$CURRENT_BRANCH${NC}"

# 2. Regenerate API stubs (Optional but recommended)
echo -e "${BLUE}Step 1: Regenerating API stubs...${NC}"
sphinx-apidoc -f -o "$DOCS_DIR/api" "$SOURCE_PACKAGE"

# 3. Build the HTML
echo -e "${BLUE}Step 2: Building HTML...${NC}"
cd "$DOCS_DIR"
make clean
make html
cd ..

# 4. Move build artifacts to a temporary location
# We do this because switching branches might hide/delete the build folder
echo -e "${BLUE}Step 3: Stashing build artifacts...${NC}"
TEMP_DIR=$(mktemp -d)
cp -r "$BUILD_OUTPUT"/* "$TEMP_DIR"

# 5. Switch to gh-pages
echo -e "${BLUE}Step 4: Switching to $TARGET_BRANCH...${NC}"
git checkout "$TARGET_BRANCH"

# 6. Update the dev folder
echo -e "${BLUE}Step 5: Updating '$TARGET_FOLDER' directory...${NC}"
mkdir -p "$TARGET_FOLDER"
# Remove old contents of dev to ensure deleted files are gone
rm -rf "$TARGET_FOLDER"/*
# Copy new contents
cp -r "$TEMP_DIR"/* "$TARGET_FOLDER"/

# 7. Commit and Push
echo -e "${BLUE}Step 6: Committing and Pushing...${NC}"
git add "$TARGET_FOLDER"
git commit -m "Deploy updated dev documentation from branch $CURRENT_BRANCH" || echo "No changes to commit"
git push origin "$TARGET_BRANCH"

# 8. Cleanup and Return
echo -e "${BLUE}Step 7: Returning to $CURRENT_BRANCH...${NC}"
rm -rf "$TEMP_DIR"
git checkout "$CURRENT_BRANCH"

echo -e "${GREEN}=== Deployment Complete! ===${NC}"
echo "Check your site at: https://das-oats.github.io/astrocook/$TARGET_FOLDER/"