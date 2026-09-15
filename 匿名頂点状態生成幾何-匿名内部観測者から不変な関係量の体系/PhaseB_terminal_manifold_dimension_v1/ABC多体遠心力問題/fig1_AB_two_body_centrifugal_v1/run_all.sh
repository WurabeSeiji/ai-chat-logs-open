#!/bin/sh
set -e
cd "$(dirname "$0")"
python3 fig1_AB_two_body_centrifugal.py
shasum -a 256 fig1_AB_two_body_centrifugal.py \
              fig1_AB_two_body_centrifugal.svg \
              fig1_AB_two_body_centrifugal.png \
              results.txt > SHA256SUMS.txt
echo "done: SHA256SUMS.txt updated"
