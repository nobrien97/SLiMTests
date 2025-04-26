#!/bin/bash

# Usage: ./getMissingFiles.sh input.txt

grep -oE '[0-9]+_[0-9]+_[0-9]+' "$1" \
    | sort -u
