#!/bin/bash
set -e
cd "$(dirname "$0")"
python3 analyze_inflation_linear_instability_20260912.py 6 8 10 2>/dev/null | tee 解析ログ_20260912.log
shasum -a 256 *.py *.md > SHA256SUMS.txt
echo "done. see REPORT.md"
