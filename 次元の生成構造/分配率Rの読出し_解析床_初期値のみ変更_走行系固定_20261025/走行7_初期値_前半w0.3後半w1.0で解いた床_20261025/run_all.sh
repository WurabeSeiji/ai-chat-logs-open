#!/bin/bash
# 走行7（2026-10-25、走行系固定）：初期値＝見えない対の重みを前半 0.3／後半 1.0 で切り替えた A で解いた床。走行系＝走行1 の run_sweep_from_parents_v1.py（バイト同一）。
set -e
cd "$(dirname "$0")"
PY="$HOME/.venvs/stage123_control_py39/bin/python3"
NS="3 4 5 6 8 10 12 15 20 25"
R1=../../分配率Rの読出し_解析床_20261025/走行1_全順序対_M_N2_隣接あり_解析床_20261025
R5=../走行5_初期値_反対辺に重みwで解いた床_w0_w0.3_20261025
START=$(date +%s)
cmp run_sweep_from_parents_v1.py $R1/run_sweep_from_parents_v1.py && echo "走行系 unchanged (byte-identical to RUN1)"
"$PY" make_initial_floor_mixed_v1.py all floors_init_mixed 0.3 1.0 $NS 2>results/log_floors_stderr.txt | tee results/log_floors.txt
"$PY" run_sweep_from_parents_v1.py floors_init_mixed 'floor_90deg_all_N{N:05d}_analytic.npz' results 'RUN7 (2026-10-25): INITIAL VALUE = floor solved with invisible-pair weight 0.3 (first half, row i<N/2) / 1.0 (second half); DYNAMICS = RUN1 program unchanged, seedless; den=N and 124' $NS 2>results/log_sweep_stderr.txt | tee results/log_sweep.txt
"$PY" compare_curves_two_runs_v1.py results/timeseries_64bit_with124_N3_N40.csv $R1/results_all/timeseries_64bit_with124_N3_N40.csv init_mixed run1 results/compare_vs_run1 $NS | tee results/log_compare.txt
"$PY" compare_curves_two_runs_v1.py results/timeseries_64bit_with124_N3_N40.csv $R5/results_w0.3/timeseries_64bit_with124_N3_N40.csv init_mixed init_w0.3 results/compare_vs_w03 $NS | tee results/log_compare_w03.txt
"$PY" plot_complex_plane_from_results_v1.py results results RUN7_init_mixed $NS | tee results/log_plot.txt
echo "elapsed: $(( $(date +%s) - START )) s"
echo 'RUN7 ALL DONE'
