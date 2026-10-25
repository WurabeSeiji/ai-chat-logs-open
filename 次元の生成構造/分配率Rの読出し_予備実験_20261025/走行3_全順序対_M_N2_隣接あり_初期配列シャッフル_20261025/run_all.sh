#!/bin/bash
# 走行3（2026-10-25）: 走行1 のプログラムに「親を読んだ直後に配列を乱数で入れ替える」1 行を足しただけ。親は走行1 の parents_all をそのまま使う。
set -e
cd "$(dirname "$0")"
PY="$HOME/.venvs/stage123_control_py39/bin/python3"
NS="3 4 5 6 8 10 12 15 20 25"
R1=../走行1_全順序対_M_N2_隣接重みあり_20261025
START=$(date +%s)
"$PY" run_sweep_from_parents_v1_shuffle.py $R1/parents_all 'parent_dense_all_N{N:05d}_20261025.npz' results_shuffle 'RUN3 (2026-10-25): RUN1 parents, initial array randomly permuted (one added line), adjacency 0/1 from positions, one_step unchanged; den=N and 124' $NS 2> results_shuffle/log_sweep_stderr.txt | tee results_shuffle/log_sweep.txt
"$PY" compare_curves_two_runs_v1.py results_shuffle/timeseries_64bit_with124_N3_N40.csv $R1/results_all/timeseries_64bit_with124_N3_N40.csv shuffled run1 results_shuffle/compare_shuffled_vs_run1 $NS | tee results_shuffle/log_compare.txt
"$PY" plot_complex_plane_from_results_v1.py results_shuffle results_shuffle RUN3_shuffled $NS | tee results_shuffle/log_plot.txt
echo "elapsed: $(( $(date +%s) - START )) s"
echo 'RUN3 ALL DONE'
