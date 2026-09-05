#!/bin/bash
# N=3..40 段1+2+3 スイープ一式の再現（親生成→スイープ→入力ゲート→図化）
set -e
cd "$(dirname "$0")"
python3 make_static_parents_N3_N40_v1.py
python3 run_N3_N40_stage123_v1.py
python3 check_sweep_inputs_v1.py
python3 plot_complex_plane_N3_N40_stage123_v1.py
echo 'SWEEP PACKAGE ALL DONE'
