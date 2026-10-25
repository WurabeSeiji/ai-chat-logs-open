#!/bin/bash
# 対照実験（2026-10-25）: 原本 N3_N40_stage123_sweep_20260905 の忠実コピーを、原本と同じ Python 3.9.6 / numpy 2.0.2 の venv で再走行し、原本の生成物と照合する。
set -e
cd "$(dirname "$0")"
PY="$HOME/.venvs/stage123_control_py39/bin/python3"   # python3.9.6 -m venv; pip install numpy==2.0.2 matplotlib==3.9.4
"$PY" make_static_parents_N3_N40_v1.py 2> results/log_parents_stderr.txt | tee results/log_parents.txt
"$PY" run_N3_N40_stage123_v1.py 2> results/log_sweep_stderr.txt | tee results/log_sweep.txt
"$PY" check_sweep_inputs_v1.py | tee results/log_gate.txt
"$PY" analyze_sweep_summary_v1.py > results/log_analyze.txt
"$PY" plot_complex_plane_N3_N40_stage123_v1.py | tee results/log_plot.txt
"$PY" compare_with_original_v1.py | tee results/log_compare.txt
echo 'CONTROL PACKAGE ALL DONE'
