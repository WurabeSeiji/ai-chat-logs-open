#!/bin/sh
set -eu
python3 sens_generate_window_states.py
python3 run_sens.py
python3 run_sens2_window.py
