#!/bin/bash
# 走行1（解析床、2026-10-25）: 論文3 の床の再生成（上三角、対照）→ N² 床の構成 → 両者を同じ one_step で走行 → 比較 → 複素平面図
set -e
cd "$(dirname "$0")"
PY="$HOME/.venvs/stage123_control_py39/bin/python3"
NS="3 4 5 6 8 10 12 15 20 25"
START=$(date +%s)
"$PY" floor_analytic_v1.py triu floors_triu $NS 2>results_control/log_floors_triu_stderr.txt | tee results_control/log_floors_triu.txt
"$PY" control_vs_paper3_manifest_v1.py | tee results_control/log_control_vs_paper3.txt
"$PY" floor_analytic_v1.py all floors_all $NS 2>results_all/log_floors_all_stderr.txt | tee results_all/log_floors_all.txt
"$PY" run_sweep_from_parents_v1.py floors_triu 'floor_90deg_triu_N{N:05d}_analytic.npz' results_triu 'BASELINE (2026-10-25): paper-3 analytic 90deg floor, upper-triangular M=N(N-1)/2, one_step unchanged, seedless; den=N and 124' $NS 2>results_triu/log_sweep_stderr.txt | tee results_triu/log_sweep.txt
"$PY" run_sweep_from_parents_v1.py floors_all 'floor_90deg_all_N{N:05d}_analytic.npz' results_all 'RUN1 analytic (2026-10-25): analytic 90deg floor on all ordered pairs incl. self (M=N^2), adjacency 0/1 kept, one_step unchanged, seedless; den=N and 124' $NS 2>results_all/log_sweep_stderr.txt | tee results_all/log_sweep.txt
"$PY" compare_curves_two_runs_v1.py results_all/timeseries_64bit_with124_N3_N40.csv results_triu/timeseries_64bit_with124_N3_N40.csv all_N2 triu results_all/compare_all_vs_triu_analytic $NS | tee results_all/log_compare.txt
"$PY" plot_complex_plane_from_results_v1.py results_triu results_triu BASELINE_triu_analytic $NS | tee results_triu/log_plot.txt
"$PY" plot_complex_plane_from_results_v1.py results_all results_all RUN1_allpairs_analytic $NS | tee results_all/log_plot.txt
echo "elapsed: $(( $(date +%s) - START )) s"
echo 'RUN1 ANALYTIC ALL DONE'
