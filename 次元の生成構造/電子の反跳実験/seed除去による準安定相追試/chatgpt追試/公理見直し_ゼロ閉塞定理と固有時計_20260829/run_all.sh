#!/bin/bash
set -e
cd "$(dirname "$0")"
python3 analyze_axioms.py 2>&1 | grep -v Warning | tee results/analyze_axioms.log
find . -type f ! -name SHA256SUMS.txt -print0 | sort -z | xargs -0 shasum -a 256 > SHA256SUMS.txt
echo "done: $(wc -l < SHA256SUMS.txt) files hashed"
