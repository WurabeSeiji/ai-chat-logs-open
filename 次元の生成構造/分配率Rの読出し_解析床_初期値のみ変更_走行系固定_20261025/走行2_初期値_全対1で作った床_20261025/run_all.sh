#!/bin/bash
# 走行2 やり直し（2026-10-25、走行系固定）：初期値＝全対 1 の A で解いた床。走行系＝走行1 の run_sweep_from_parents_v1.py（バイト同一、A は ea,eb から標準規則）。
set -e
cd "$(dirname "$0")"
PY="$HOME/.venvs/stage123_control_py39/bin/python3"
NS="3 4 5 6 8 10 12 15 20 25"
R1=../../分配率Rの読出し_解析床_20261025/走行1_全順序対_M_N2_隣接あり_解析床_20261025
START=$(date +%s)
cmp run_sweep_from_parents_v1.py $R1/run_sweep_from_parents_v1.py && echo "走行系 unchanged (byte-identical to RUN1)"
"$PY" make_initial_floor_noadj_v1.py all floors_init_noadj $NS 2>results/log_floors_stderr.txt | tee results/log_floors.txt
"$PY" run_sweep_from_parents_v1.py floors_init_noadj 'floor_90deg_all_N{N:05d}_analytic.npz' results 'RUN2 redo (2026-10-25): INITIAL VALUE = floor solved with all-pairs-1 A; DYNAMICS = RUN1 program unchanged (standard adjacency from positions, one_step unchanged), seedless; den=N and 124' $NS 2>results/log_sweep_stderr.txt | tee results/log_sweep.txt
"$PY" compare_curves_two_runs_v1.py results/timeseries_64bit_with124_N3_N40.csv $R1/results_all/timeseries_64bit_with124_N3_N40.csv init_noadj run1 results/compare_vs_run1 $NS | tee results/log_compare.txt
"$PY" plot_complex_plane_from_results_v1.py results results RUN2redo_init_noadj $NS | tee results/log_plot.txt
echo "elapsed: $(( $(date +%s) - START )) s"
echo 'RUN2 REDO ALL DONE'
