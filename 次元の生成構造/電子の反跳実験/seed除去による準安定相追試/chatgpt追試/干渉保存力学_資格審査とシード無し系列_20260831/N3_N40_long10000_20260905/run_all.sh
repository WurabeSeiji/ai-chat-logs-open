#!/bin/bash
# 全系列10000歩スイープ: 走行→対照ゲート→軸結晶化→終状態最終検査
set -euo pipefail
cd "$(dirname "$0")"
python3 run_N3_N40_steps10000_v1.py
python3 check_control_sweep10000_v1.py
python3 analyze_axis_crystallization_v1.py
python3 check_final_states_10000_v1.py
python3 check_phase_spacing_statistics_v1.py
echo "run_all done"
