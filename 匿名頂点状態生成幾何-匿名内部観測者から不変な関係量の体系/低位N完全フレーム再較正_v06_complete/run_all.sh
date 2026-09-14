#!/bin/sh
set -eu
python3 run_lowN_v06_recalibration.py
python3 make_sha256.py
