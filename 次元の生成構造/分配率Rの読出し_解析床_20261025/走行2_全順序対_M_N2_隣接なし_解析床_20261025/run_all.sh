#!/bin/bash
# 走行2（解析床、2026-10-25）: 全順序対 M=N²、隣接重みづけなし（全対 1）。床の構成→走行→走行1（解析床・隣接あり）との曲線照合→複素平面図
set -e
cd "$(dirname "$0")"
PY="$HOME/.venvs/stage123_control_py39/bin/python3"
NS="3 4 5 6 8 10 12 15 20 25"
R1=../走行1_全順序対_M_N2_隣接あり_解析床_20261025
START=$(date +%s)
"$PY" floor_analytic_v2_noadj.py all floors_all_noadj $NS 2>results_all_noadj/log_floors_stderr.txt | tee results_all_noadj/log_floors.txt
"$PY" run_sweep_from_parents_v2_noadj.py floors_all_noadj 'floor_90deg_all_N{N:05d}_analytic.npz' results_all_noadj 'RUN2 analytic (2026-10-25): analytic 90deg floor on all ordered pairs (M=N^2), NO adjacency weighting (all pairs 1), one_step unchanged, seedless; den=N and 124' $NS 2>results_all_noadj/log_sweep_stderr.txt | tee results_all_noadj/log_sweep.txt
"$PY" compare_curves_two_runs_v1.py results_all_noadj/timeseries_64bit_with124_N3_N40.csv $R1/results_all/timeseries_64bit_with124_N3_N40.csv noadj adj results_all_noadj/compare_noadj_vs_run1_adj $NS | tee results_all_noadj/log_compare.txt
"$PY" plot_complex_plane_from_results_v1.py results_all_noadj results_all_noadj RUN2_allpairs_noadj_analytic $NS | tee results_all_noadj/log_plot.txt
echo "elapsed: $(( $(date +%s) - START )) s"
echo 'RUN2 ANALYTIC ALL DONE'
