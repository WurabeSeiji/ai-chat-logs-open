#!/bin/bash
# SSB検証（同一床＋異なる微小シード→縮退真空）の再現。厳密one_stepの制御実験。
set -e
cd "$(dirname "$0")"
python3 verify_ssb_seed_vacua_20260912.py 6 2>/dev/null | tee run.log
shasum -a 256 *.py *.md > SHA256SUMS.txt
echo "done. see REPORT.md"
