#!/bin/bash
# 走行3 やり直し（2026-10-25、走行系固定）：初期値＝走行1 の床を乱数で入れ替えた npz。走行系＝走行1 の run_sweep_from_parents_v1.py（バイト同一）。
set -e
cd "$(dirname "$0")"
PY="$HOME/.venvs/stage123_control_py39/bin/python3"
NS="3 4 5 6 8 10 12 15 20 25"
R1=../../分配率Rの読出し_解析床_20261025/走行1_全順序対_M_N2_隣接あり_解析床_20261025
START=$(date +%s)
cmp run_sweep_from_parents_v1.py $R1/run_sweep_from_parents_v1.py && echo "走行系 unchanged (byte-identical to RUN1)"
"$PY" make_initial_shuffled_v1.py $R1/floors_all floors_init_shuffled $NS | tee results/log_initial_values.txt
"$PY" run_sweep_from_parents_v1.py floors_init_shuffled 'floor_90deg_all_N{N:05d}_analytic.npz' results 'RUN3 redo (2026-10-25): INITIAL VALUE = RUN1 analytic floor with array randomly permuted (npz); DYNAMICS = RUN1 program unchanged, seedless; den=N and 124' $NS 2>results/log_sweep_stderr.txt | tee results/log_sweep.txt
"$PY" compare_curves_two_runs_v1.py results/timeseries_64bit_with124_N3_N40.csv $R1/results_all/timeseries_64bit_with124_N3_N40.csv init_shuffled run1 results/compare_vs_run1 $NS | tee results/log_compare.txt
"$PY" plot_complex_plane_from_results_v1.py results results RUN3redo_init_shuffled $NS | tee results/log_plot.txt
echo "elapsed: $(( $(date +%s) - START )) s"
echo 'RUN3 REDO ALL DONE'
