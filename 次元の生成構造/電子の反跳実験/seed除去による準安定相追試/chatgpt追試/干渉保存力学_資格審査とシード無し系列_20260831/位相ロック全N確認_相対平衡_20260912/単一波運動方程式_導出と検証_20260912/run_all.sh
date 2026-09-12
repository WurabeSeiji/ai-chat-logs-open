#!/bin/bash
# 単一波EOM導出＋90°格子必然性の検証を再現。読出しのみ・物理無変更。
set -e
cd "$(dirname "$0")"
python3 verify_single_wave_eom_20260912.py 6 22 40 2>/dev/null | tee 検証ログ_20260912.log
python3 verify_90deg_necessity_20260912.py 6 22 40 2>/dev/null | tee 検証ログ_90度必然性_20260912.log
shasum -a 256 *.py *.md > SHA256SUMS.txt
echo "done. see REPORT.md"
