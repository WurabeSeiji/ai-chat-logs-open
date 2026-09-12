#!/bin/bash
# 全N位相ロック確認の再現。読出しのみ・物理無変更。
set -e
cd "$(dirname "$0")"
python3 analyze_phaselock_allN_20260912.py | tee run.log
python3 figure_and_convergence_20260912.py 2>/dev/null | tee -a run.log
shasum -a 256 *.py *.csv *.png REPORT.md > SHA256SUMS.txt
echo "done. see REPORT.md / fig_phaselock_allN.png / phaselock_allN.csv"
