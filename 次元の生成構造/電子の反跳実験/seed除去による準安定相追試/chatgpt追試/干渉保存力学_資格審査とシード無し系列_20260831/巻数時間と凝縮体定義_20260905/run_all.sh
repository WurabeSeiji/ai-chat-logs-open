#!/bin/bash
# 巻数時間仮説の検定（読み出しのみ）
set -euo pipefail
cd "$(dirname "$0")"
python3 analyze_winding_rates_v1.py
python3 check_aliased_winding_v1.py
python3 check_clock_ratio_v1.py
echo "run_all done"
