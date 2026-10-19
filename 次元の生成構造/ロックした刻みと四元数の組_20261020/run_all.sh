#!/bin/bash
set -e
cd "$(dirname "$0")"
python3 run_N40_locked_step_20261020.py control
python3 run_N40_locked_step_20261020.py q8
python3 run_N40_locked_step_20261020.py q4
python3 analyze_locked_pairing_20261020.py
python3 plot_locked_inflation_complex_20261020.py
shasum -a 256 run_N40_locked_step_20261020.py analyze_locked_pairing_20261020.py plot_locked_inflation_complex_20261020.py figures/* results/*.json results/*.md results/*.npz > SHA256SUMS.txt
