#!/bin/bash
set -e
cd "$(dirname "$0")"
python3 construct_and_analyze.py
find . -type f ! -name SHA256SUMS.txt -print0 | sort -z | xargs -0 shasum -a 256 > SHA256SUMS.txt
echo "done: $(wc -l < SHA256SUMS.txt) files hashed"
