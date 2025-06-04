#!/bin/bash
# Git operations script

cd "/Users/fvrodriguez/Library/CloudStorage/GoogleDrive-frodrig4@ncsu.edu/My Drive/repos/BzeaSeq"

echo "Current status:"
git status --short

echo -e "\nCommitting updated codemcp.toml..."
git add codemcp.toml
git commit -m "Add git commands to codemcp.toml configuration"

echo -e "\nPulling from origin/main..."
git pull origin main

echo -e "\nPushing to origin/main..."
git push origin main

echo -e "\nFinal status:"
git status --short

# Clean up
rm git_operations.sh
