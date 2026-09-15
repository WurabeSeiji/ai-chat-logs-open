#!/bin/sh
cd "$(dirname "$0")" || exit 1
python3 regen_states_and_gate.py && python3 run_svd_analysis.py
