#!/bin/bash

INPUT_FILE="${1:-checksums.tsv}"
DOWNLOAD_DIR="${2:-.}"

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Verify inputs
[[ ! -f "$INPUT_FILE" ]] && { echo -e "${RED}Error: Input file not found${NC}" >&2; exit 1; }
[[ ! -d "$DOWNLOAD_DIR" ]] && { echo -e "${RED}Error: Directory not found${NC}" >&2; exit 1; }

echo -e "Checking files in: $DOWNLOAD_DIR"
echo -e "Using checksums from: $INPUT_FILE\n"

all_good=true
declare -i total=0 verified=0 missing=0 failed=0

while IFS=$'\t' read -r filename _ _ _ expected_md5 _; do
    # Skip empty lines
    [[ -z "$filename" ]] && continue
    
    ((total++))
    filepath="$DOWNLOAD_DIR/$filename"
    
    echo -n "Checking $filename... "
    
    if [[ ! -f "$filepath" ]]; then
        echo -e "${YELLOW}Not Found${NC}"
        ((missing++))
        all_good=false
        continue
    fi
    
    actual_md5=$(md5sum "$filepath" | awk '{print $1}')
    
    if [[ "$actual_md5" == "$expected_md5" ]]; then
        echo -e "${GREEN}OK${NC}"
        ((verified++))
    else
        echo -e "${RED}FAIL${NC}"
        echo -e "  Expected: $expected_md5"
        echo -e "  Actual:   $actual_md5"
        ((failed++))
        all_good=false
    fi
done < "$INPUT_FILE"

# Summary
echo -e "\nResults:"
echo -e "${GREEN}Verified: $verified${NC}"
echo -e "${RED}Failed: $failed${NC}"
echo -e "${YELLOW}Missing: $missing${NC}"
echo -e "Total files checked: $total"

exit $((all_good ? 0 : 1))