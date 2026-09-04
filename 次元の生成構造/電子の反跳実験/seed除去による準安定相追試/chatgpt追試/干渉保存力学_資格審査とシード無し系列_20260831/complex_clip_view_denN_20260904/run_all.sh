#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
python3 program/plot_complex_clip_denN.py \
  --src "../ChatGPT_denominator_controls_N3_N33_legacyparent_20260903/results_2000steps" \
  --out results 2>&1 | tee results/run_all.log
