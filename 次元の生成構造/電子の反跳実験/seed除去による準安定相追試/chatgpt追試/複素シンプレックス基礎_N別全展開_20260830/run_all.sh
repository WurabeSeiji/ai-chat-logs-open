#!/bin/bash
set -e
set -o pipefail
cd "$(dirname "$0")"
python3 program/make_N3.py
python3 program/make_N4.py
for n in 5 6 7 8 9 10 11 12 13 14 15 16; do python3 program/make_Ngeneric.py $n; done
python3 program/assemble.py
find . -type f ! -name SHA256SUMS.txt -print0 | sort -z | xargs -0 shasum -a 256 > SHA256SUMS.txt
echo "done"
