#!/bin/bash
# N=30..39 2000step バッチ検証の再現（走行→一括解析）。物理正本を忠実コピーで exec。
set -e
cd "$(dirname "$0")"
python3 wrapper_run_N30_N39_2000_v1.py 2>&1 | tee 実行ログ_20260912.log
python3 analyze_batch_N30_N39_20260912.py 2>/dev/null | tee 分析ログ_20260912.log
shasum -a 256 *.py *.png *.csv REPORT.md results/*.npz > SHA256SUMS.txt
echo "done. REPORT.md / phaselock_batch_N30_N39.csv / fig_batch_lockstep_vs_pred.png"
