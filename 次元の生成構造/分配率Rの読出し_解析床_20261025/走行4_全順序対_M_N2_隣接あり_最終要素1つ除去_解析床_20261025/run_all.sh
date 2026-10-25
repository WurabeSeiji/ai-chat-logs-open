#!/bin/bash
# 走行4（解析床、2026-10-25）: 走行1 のプログラムに「初期化後の配列の最終要素を 1 つ除く」1 行を足しただけ。床は走行1 の floors_all をそのまま使う。
set -e
cd "$(dirname "$0")"
PY="$HOME/.venvs/stage123_control_py39/bin/python3"
NS="3 4 5 6 8 10 12 15 20 25"
R1=../走行1_全順序対_M_N2_隣接あり_解析床_20261025
START=$(date +%s)
"$PY" run_sweep_from_parents_v1_droplast.py $R1/floors_all 'floor_90deg_all_N{N:05d}_analytic.npz' results_droplast 'RUN4 analytic (2026-10-25): RUN1 analytic floors (M=N^2) with the LAST element removed (one added line), adjacency 0/1 from positions, one_step unchanged, seedless; den=N and 124' $NS 2>results_droplast/log_sweep_stderr.txt | tee results_droplast/log_sweep.txt
"$PY" compare_curves_two_runs_v1.py results_droplast/timeseries_64bit_with124_N3_N40.csv $R1/results_all/timeseries_64bit_with124_N3_N40.csv droplast run1 results_droplast/compare_droplast_vs_run1 $NS | tee results_droplast/log_compare.txt
"$PY" plot_complex_plane_from_results_v1.py results_droplast results_droplast RUN4_droplast_analytic $NS | tee results_droplast/log_plot.txt
echo "elapsed: $(( $(date +%s) - START )) s"
echo 'RUN4 ANALYTIC ALL DONE'
