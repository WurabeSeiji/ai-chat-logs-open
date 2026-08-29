#!/bin/bash
set -e
set -o pipefail
cd "$(dirname "$0")"
export MPLBACKEND=Agg
echo "== PASS1 =="; python3 program/pass1_make_parents.py 2>&1 | tee results/pass1.log
echo "== PASS2 =="
for N in 4 5 6 7 8; do python3 program/pass2_run.py $N 2>&1 | tee results/pass2_N$N.log; done
echo "== PASS3 =="; python3 program/pass3_figures.py 2>&1 | tee results/pass3.log
find . -type f ! -name SHA256SUMS.txt -print0 | sort -z | xargs -0 shasum -a 256 > SHA256SUMS.txt
echo "ALLDONE $(wc -l < SHA256SUMS.txt) files"
