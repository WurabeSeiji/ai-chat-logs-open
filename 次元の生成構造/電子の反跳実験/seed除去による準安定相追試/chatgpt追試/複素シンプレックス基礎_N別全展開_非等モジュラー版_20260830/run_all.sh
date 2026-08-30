#!/bin/bash
set -e
set -o pipefail
cd "$(dirname "$0")"
for n in 3 4 5 6 7 8 9 10 11 12 13 14 15 16; do python3 program/make_N.py $n; done
python3 program/assemble.py
bash build_pdf.sh
find . -type f ! -name SHA256SUMS.txt -print0 | sort -z | xargs -0 shasum -a 256 > SHA256SUMS.txt
echo "done"
