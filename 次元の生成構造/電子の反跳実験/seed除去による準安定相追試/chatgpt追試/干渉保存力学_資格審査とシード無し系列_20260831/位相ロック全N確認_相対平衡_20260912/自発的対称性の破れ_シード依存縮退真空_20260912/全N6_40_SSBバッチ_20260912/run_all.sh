#!/bin/bash
set -e
cd "$(dirname "$0")"
python3 run_ssb_batch_N6_N40_20260912.py 2>/dev/null | tee run.log
shasum -a 256 *.py *.md per_N/*.txt > SHA256SUMS.txt
echo "done. per_N/ と ssb_batch_summary.csv 参照"
