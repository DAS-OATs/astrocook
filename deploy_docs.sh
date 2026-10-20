#!/bin/bash
# Builds the website (build_site.sh) and pushes it to gh-pages by hand:
#   ./deploy_docs.sh stable   docs of the checked-out release at the root
#                             (check out its tag first: git checkout v2.0.3)
#   ./deploy_docs.sh dev      docs of the checked-out master under /dev/
# The same happens automatically (.github/workflows/docs.yml) when a
# release is published and after each merge to master; this script stays
# as a fallback. gh-pages is updated in a separate worktree, so files of
# this working directory (sessions, caches) cannot end up on the site.

set -e

TARGET="$1"
if [ "$TARGET" != stable ] && [ "$TARGET" != dev ]; then
    echo "Usage: ./deploy_docs.sh stable|dev"
    exit 1
fi
if ! command -v sphinx-apidoc &> /dev/null; then
    echo "Error: sphinx-apidoc not found. Install it with: pip install -e \".[docs]\""
    exit 1
fi
if [[ -n $(git status -s) ]]; then
    echo "Error: Your working directory is not clean."
    echo "Please commit or stash your changes before running this script."
    exit 1
fi

TEMP_DIR=$(mktemp -d)
trap 'git worktree remove --force "$TEMP_DIR/gh-pages" 2>/dev/null; rm -rf "$TEMP_DIR"' EXIT

# 1. Build the site; switcher.json at the root always comes from master
./build_site.sh "$TEMP_DIR/site" "$TARGET"
git fetch origin master gh-pages
git show origin/master:switcher.json > "$TEMP_DIR/switcher.json"

# 2. Replace the root (but /dev/) or /dev/ on gh-pages
git worktree add "$TEMP_DIR/gh-pages" origin/gh-pages
cd "$TEMP_DIR/gh-pages"
if [ "$TARGET" = stable ]; then
    find . -mindepth 1 -maxdepth 1 ! -name .git ! -name dev -exec rm -rf {} +
    cp -r "$TEMP_DIR/site"/. .
else
    rm -rf dev
    mkdir dev
    cp -r "$TEMP_DIR/site"/. dev/
fi
cp "$TEMP_DIR/switcher.json" switcher.json

# 3. Commit and push
git add -A
git commit -m "Deploy $TARGET docs by hand" || echo "No changes to commit"
git push origin HEAD:gh-pages

echo "Done: https://das-oats.github.io/astrocook/"
