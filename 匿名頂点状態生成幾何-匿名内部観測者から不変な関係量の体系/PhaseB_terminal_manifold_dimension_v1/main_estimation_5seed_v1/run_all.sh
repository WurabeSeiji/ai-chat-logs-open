#!/bin/sh
set -eu
python3 run_main_estimation.py
python3 make_figures.py
