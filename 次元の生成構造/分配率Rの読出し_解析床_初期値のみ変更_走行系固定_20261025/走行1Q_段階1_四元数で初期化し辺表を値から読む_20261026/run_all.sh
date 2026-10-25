#!/bin/bash
# 走行1Q 段階1（2026-10-26、走行系固定）：初期値＝走行1 の床を四元数 q=e^(iπv/N) r e^(jπw/N) で符号化し、辺表 ea, eb と波 Z0 を q の値だけから復号したもの。
# 走行系＝走行1 の run_sweep_from_parents_v1.py（バイト同一、A は ea,eb から標準規則）。
set -e
cd "$(dirname "$0")"
PY="$HOME/.venvs/stage123_control_py39/bin/python3"
NS="3 4 5 6 8 10 12 15 20 25"
R1=../../分配率Rの読出し_解析床_20261025/走行1_全順序対_M_N2_隣接あり_解析床_20261025
START=$(date +%s)
cmp run_sweep_from_parents_v1.py $R1/run_sweep_from_parents_v1.py && echo "走行系 unchanged (byte-identical to RUN1)"
"$PY" make_initial_quaternion_v1.py $R1/floors_all floors_init_quaternion $NS 2>results/log_floors_stderr.txt | tee results/log_floors.txt
"$PY" run_sweep_from_parents_v1.py floors_init_quaternion 'floor_90deg_all_N{N:05d}_analytic.npz' results 'RUN1Q stage1 (2026-10-26): INITIAL VALUE = RUN1 floor encoded as quaternion q=e^(i pi v/N) r e^(j pi w/N), ea/eb/Z0 decoded from q values only; DYNAMICS = RUN1 program unchanged, seedless; den=N and 124' $NS 2>results/log_sweep_stderr.txt | tee results/log_sweep.txt
"$PY" compare_curves_two_runs_v1.py results/timeseries_64bit_with124_N3_N40.csv $R1/results_all/timeseries_64bit_with124_N3_N40.csv init_quaternion run1 results/compare_vs_run1 $NS | tee results/log_compare.txt
"$PY" plot_hperp_grid_from_timeseries_v1.py results 'RUN1Q stage1 (2026-10-26): RUN1 floor via quaternion, ea/eb/Z0 decoded from values; dynamics = RUN1 unchanged; den=N and 124' | tee results/log_grid.txt
"$PY" plot_complex_plane_from_results_v1.py results results RUN1Q_stage1_quaternion $NS | tee results/log_plot.txt
echo "elapsed: $(( $(date +%s) - START )) s"
echo 'RUN1Q STAGE1 ALL DONE'
