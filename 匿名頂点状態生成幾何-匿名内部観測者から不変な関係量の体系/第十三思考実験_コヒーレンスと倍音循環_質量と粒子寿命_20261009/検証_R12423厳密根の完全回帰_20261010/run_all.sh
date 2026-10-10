#!/bin/sh
# 厳密根 R_{124,23} の完全回帰: 記録 → 図化 を順に実行する。
# 前提: Python 3 と numpy, matplotlib。原本プログラムは同じ親フォルダの
#       検証_R070交換の正規化と振幅_20261009/original_copy/ から取る（無ければ手で置く）。
set -e
cd "$(dirname "$0")"
ORIG=original_copy/run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1.py
if [ ! -f "$ORIG" ]; then
  mkdir -p original_copy
  cp "../検証_R070交換の正規化と振幅_20261009/original_copy/run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1.py" "$ORIG"
fi
echo "f815320f5632ae1b23ccade3a53b01e9110ad770da8407e2b10fa6065ef1695c  $ORIG" | sha256sum -c -
# [S7] の summary.json（対照に使う）。無ければ対照はスキップされる。
if [ ! -f results/reference_summary_R070_S7.json ]; then
  mkdir -p results
  cp "../検証_R070交換の正規化と振幅_20261009/results/summary.json" results/reference_summary_R070_S7.json || true
fi
python3 run_exact_root_recurrence_R12423.py > results/run_stdout.txt
python3 plot_exact_root_recurrence_R12423.py > figures/plot_stdout.txt
echo "done. compare with SHA256SUMS (sha256sum -c SHA256SUMS); CSV may differ in the last digit across numpy versions"
