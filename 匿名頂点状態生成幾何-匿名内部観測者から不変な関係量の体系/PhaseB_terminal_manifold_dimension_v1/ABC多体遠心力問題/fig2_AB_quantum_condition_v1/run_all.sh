#!/bin/sh
set -e
cd "$(dirname "$0")"
python3 fig2_AB_quantum_condition.py
shasum -a 256 fig2_AB_quantum_condition.py \
              fig2_AB_quantum_condition.svg \
              fig2_AB_quantum_condition.png \
              results.txt > SHA256SUMS.txt
echo "done: SHA256SUMS.txt updated"
