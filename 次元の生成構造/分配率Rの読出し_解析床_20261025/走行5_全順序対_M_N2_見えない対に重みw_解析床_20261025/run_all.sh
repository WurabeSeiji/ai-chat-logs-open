#!/bin/bash
# 走行5（解析床、2026-10-25）: 見えない対（添字を共有しない対）に重み w。w ∈ {0, 0.3}。w=0 は走行1 との一致確認（対照）。
set -e
cd "$(dirname "$0")"
PY="$HOME/.venvs/stage123_control_py39/bin/python3"
NS="3 4 5 6 8 10 12 15 20 25"
R1=../走行1_全順序対_M_N2_隣接あり_解析床_20261025
START=$(date +%s)
for W in 0 0.3; do
  mkdir -p floors_w$W results_w$W
  "$PY" floor_analytic_v5_w.py all floors_w$W $W $NS 2>results_w$W/log_floors_stderr.txt | tee results_w$W/log_floors.txt
  "$PY" run_sweep_from_parents_v5_w.py floors_w$W 'floor_90deg_all_N{N:05d}_analytic.npz' results_w$W "RUN5 analytic (2026-10-25): analytic floor on all ordered pairs (M=N^2), invisible pairs weighted w=$W (touching=1), one_step unchanged, seedless; den=N and 124" $W $NS 2>results_w$W/log_sweep_stderr.txt | tee results_w$W/log_sweep.txt
  "$PY" plot_complex_plane_from_results_v1.py results_w$W results_w$W RUN5_w$W $NS | tee results_w$W/log_plot.txt
done
"$PY" compare_curves_two_runs_v1.py results_w0/timeseries_64bit_with124_N3_N40.csv $R1/results_all/timeseries_64bit_with124_N3_N40.csv w0 run1 results_w0/compare_w0_vs_run1 $NS | tee results_w0/log_compare.txt
"$PY" compare_curves_two_runs_v1.py results_w0.3/timeseries_64bit_with124_N3_N40.csv results_w0/timeseries_64bit_with124_N3_N40.csv w0.3 w0 results_w0.3/compare_w03_vs_w0 $NS | tee results_w0.3/log_compare.txt
echo "elapsed: $(( $(date +%s) - START )) s"
echo 'RUN5 ALL DONE'
