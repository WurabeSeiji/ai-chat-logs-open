#!/bin/bash
set -e
cd "$(dirname "$0")"
python3 run_N40_locked_step_20261020.py control
python3 run_N40_locked_step_20261020.py q8
python3 run_N40_locked_step_20261020.py q4
python3 analyze_locked_pairing_20261020.py
shasum -a 256 run_N40_locked_step_20261020.py analyze_locked_pairing_20261020.py results/*.json results/*.md results/*.npz > SHA256SUMS.txt
