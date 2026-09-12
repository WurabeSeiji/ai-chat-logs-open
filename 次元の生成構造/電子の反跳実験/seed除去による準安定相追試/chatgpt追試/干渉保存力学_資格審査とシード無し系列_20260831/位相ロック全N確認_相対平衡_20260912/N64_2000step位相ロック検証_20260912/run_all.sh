#!/bin/bash
# N=64 2000step 位相ロック検証の再現。読出しのみ・物理無変更。
set -e
cd "$(dirname "$0")"
python3 analyze_phaselock_N64_2000_20260912.py 2>/dev/null | tee run.log
shasum -a 256 *.py *.png REPORT.md > SHA256SUMS.txt
echo "done. see REPORT.md / fig_phaselock_N64_2000.png"
