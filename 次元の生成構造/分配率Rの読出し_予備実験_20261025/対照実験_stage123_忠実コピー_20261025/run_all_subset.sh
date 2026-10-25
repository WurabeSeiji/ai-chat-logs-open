#!/bin/bash
# 対照実験 案B（2026-10-25）: 代表 12 N × den∈{N,124} × 500 step。親は既存（make_static_parents は既存を検証のみ）。
set -e
cd "$(dirname "$0")"
PY="$HOME/.venvs/stage123_control_py39/bin/python3"
START=$(date +%s)
"$PY" make_static_parents_N3_N40_v1.py 2> results/log_parents_verify_stderr.txt | tee results/log_parents_verify.txt | tail -3 || true
"$PY" run_subset_stage123_v1.py 2> results/log_subset_sweep_stderr.txt | tee results/log_subset_sweep.txt
"$PY" check_subset_inputs_v1.py | tee results/log_subset_gate.txt
"$PY" analyze_subset_summary_v1.py > results/log_subset_analyze.txt
"$PY" plot_complex_plane_subset_v1.py | tee results/log_subset_plot.txt
"$PY" compare_with_original_v1.py | tee results/log_compare.txt
"$PY" compare_curves_with_original_v1.py | tee results/log_compare_curves.txt
echo "elapsed: $(( $(date +%s) - START )) s"
echo 'CONTROL SUBSET B ALL DONE'
