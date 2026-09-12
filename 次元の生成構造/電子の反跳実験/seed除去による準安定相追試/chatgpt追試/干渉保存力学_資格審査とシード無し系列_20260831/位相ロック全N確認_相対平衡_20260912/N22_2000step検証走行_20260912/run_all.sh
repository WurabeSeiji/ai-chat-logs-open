#!/bin/bash
# N=22 2000step 検証走行の再現（走行→解析）。物理正本を忠実コピーで exec。
set -e
cd "$(dirname "$0")"
python3 wrapper_run_N22_2000_v1.py 2>&1 | tee 実行ログ_20260912.log
python3 analyze_phaselock_N22_2000_20260912.py 2>/dev/null | tee 分析ログ_20260912.log
shasum -a 256 *.py *.png REPORT.md results/*.npz > SHA256SUMS.txt
echo "done. see REPORT.md / fig_phaselock_N22_2000.png"
