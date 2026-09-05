#!/bin/bash
# N=5 10000歩走行＋対照ゲート＋走行後分析
set -euo pipefail
cd "$(dirname "$0")"
python3 run_N5_den5_steps10000_v1.py
python3 check_control_first501_v1.py
python3 analyze_N5_long10000_v1.py
echo "run_all done"
