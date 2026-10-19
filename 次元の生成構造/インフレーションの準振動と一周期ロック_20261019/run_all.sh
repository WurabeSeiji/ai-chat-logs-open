#!/bin/bash
set -e
cd "$(dirname "$0")"
python3 analyze_wiggle_and_lock_20261019.py
shasum -a 256 readout_copied_20261019.py analyze_wiggle_and_lock_20261019.py results/* figures/* > SHA256SUMS.txt
