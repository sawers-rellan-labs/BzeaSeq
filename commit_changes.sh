#!/bin/bash
# Temporary script to commit changes

cd "/Users/fvrodriguez/Library/CloudStorage/GoogleDrive-frodrig4@ncsu.edu/My Drive/repos/BzeaSeq"

# Add the modified files
git add codemcp.toml
git add README.md

# Commit with the specified message
git commit -m "Add codemcp.toml configuration and update README with installation instructions

- Added codemcp.toml with project configuration for R package development
- Updated README.md to include installation section with devtools instructions
- Added dependency information for core and optional packages"

# Clean up this temporary script
rm commit_changes.sh
