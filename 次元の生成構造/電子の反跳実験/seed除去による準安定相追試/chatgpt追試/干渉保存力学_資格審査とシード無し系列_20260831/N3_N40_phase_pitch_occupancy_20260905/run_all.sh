#!/bin/bash
# 位相ピッチ・占有・振幅保護の全系列調査（読み出しのみ）
set -euo pipefail
cd "$(dirname "$0")"
python3 analyze_phase_pitch_occupancy_v1.py
python3 check_N5_long_geometry_v1.py
python3 analyze_transition_migration_v1.py
python3 check_axis_structure_v1.py
python3 check_closure_factorization_v1.py
python3 check_certified_factorization_v1.py
python3 check_parity_class_v1.py
echo "run_all done"
