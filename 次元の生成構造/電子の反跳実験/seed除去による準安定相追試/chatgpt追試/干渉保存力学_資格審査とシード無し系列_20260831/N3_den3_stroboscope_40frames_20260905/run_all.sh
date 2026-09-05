#!/bin/bash
# N=3・Δτ=2π/3 ストロボ40コマ図化＋所見の数値検証（読み出しのみ）
set -euo pipefail
cd "$(dirname "$0")"
python3 plot_complex_plane_N3_den3_frames_v1.py
python3 check_strobe_alternation_v1.py
python3 check_phase_differences_v1.py
python3 check_rotation_lock_v1.py
echo "run_all done"
