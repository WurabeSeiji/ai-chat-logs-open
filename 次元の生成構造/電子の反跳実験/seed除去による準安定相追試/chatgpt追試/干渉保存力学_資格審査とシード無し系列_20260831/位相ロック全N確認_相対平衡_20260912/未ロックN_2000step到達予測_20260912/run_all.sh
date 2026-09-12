#!/bin/bash
# 未ロックN抽出＋2000step到達予測。読出しのみ・走行なし。
set -e
cd "$(dirname "$0")"
python3 predict_lock_by_2000_20260912.py 2>/dev/null | tee run.log
shasum -a 256 *.py *.png REPORT.md > SHA256SUMS.txt
echo "done. see REPORT.md / fig_predict_lock_by_2000.png"
