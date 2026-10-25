#!/bin/bash
# 走行2（2026-10-25）: 全順序対 M=N²、隣接重みづけなし（全対 1）。親生成→スイープ→走行1（隣接あり）との曲線照合→複素平面図。
set -e
cd "$(dirname "$0")"
PY="$HOME/.venvs/stage123_control_py39/bin/python3"
NS="3 4 5 6 8 10 12 15 20 25"
START=$(date +%s)
"$PY" parents_dense_v2_noadj.py all parents_all_noadj $NS 2> results_all_noadj/log_parents_stderr.txt | tee results_all_noadj/log_parents.txt
"$PY" run_sweep_from_parents_v2_noadj.py parents_all_noadj 'parent_dense_all_N{N:05d}_20261025.npz' results_all_noadj 'RUN2 (2026-10-25): all ordered pairs incl. self (M=N^2), NO adjacency weighting (all pairs 1), one_step unchanged; den=N and 124' $NS 2> results_all_noadj/log_sweep_stderr.txt | tee results_all_noadj/log_sweep.txt
"$PY" compare_curves_two_runs_v1.py results_all_noadj/timeseries_64bit_with124_N3_N40.csv ../走行1_全順序対_M_N2_隣接重みあり_20261025/results_all/timeseries_64bit_with124_N3_N40.csv noadj adj results_all_noadj/compare_noadj_vs_run1_adj $NS | tee results_all_noadj/log_compare.txt
"$PY" plot_complex_plane_from_results_v1.py results_all_noadj results_all_noadj RUN2_allpairs_noadj $NS | tee results_all_noadj/log_plot.txt
echo "elapsed: $(( $(date +%s) - START )) s"
echo 'RUN2 ALL DONE'
