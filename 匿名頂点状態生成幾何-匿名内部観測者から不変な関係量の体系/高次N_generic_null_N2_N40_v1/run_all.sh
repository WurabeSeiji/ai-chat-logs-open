#!/bin/sh
set -eu
python3 run_highN_generic_null_N2_N40.py
python3 make_plots.py
python3 make_sha256.py
