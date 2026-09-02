#!/bin/sh
# Full reproduction: Stage A -> aggregate -> Stage B -> aggregate -> figures.
# Runs are checkpointed per (N, D, stage); completed runs are skipped.
set -e
cd "$(dirname "$0")"
python3 program/run_sweep.py --stage A
python3 program/analyze_sweep.py
python3 program/run_sweep.py --stage B
python3 program/analyze_sweep.py
python3 program/plot_sweep.py
