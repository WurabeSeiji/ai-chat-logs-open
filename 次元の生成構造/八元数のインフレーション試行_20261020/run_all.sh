#!/bin/bash
set -e
cd "$(dirname "$0")"
python3 octonion_inflation_trials_20261020.py
shasum -a 256 octonion_inflation_trials_20261020.py results/* > SHA256SUMS.txt
