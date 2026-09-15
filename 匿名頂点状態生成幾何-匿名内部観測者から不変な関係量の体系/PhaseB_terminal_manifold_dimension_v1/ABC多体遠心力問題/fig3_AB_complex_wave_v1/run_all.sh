#!/bin/sh
set -e
cd "$(dirname "$0")"
python3 fig3_AB_complex_wave.py
shasum -a 256 fig3_AB_complex_wave.py \
              fig3_AB_complex_wave.svg \
              fig3_AB_complex_wave.png \
              results.txt > SHA256SUMS.txt
echo "done: SHA256SUMS.txt updated"
