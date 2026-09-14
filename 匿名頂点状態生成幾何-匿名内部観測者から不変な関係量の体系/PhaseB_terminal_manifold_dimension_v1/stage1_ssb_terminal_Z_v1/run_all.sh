#!/bin/sh
set -eu
python3 run_all_save_terminal_Z.py
python3 compare_against_canonical.py
python3 evaluate_effective_samples.py
