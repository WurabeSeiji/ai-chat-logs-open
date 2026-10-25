#!/bin/bash
# 走行1Q 段階2（2026-10-26）：四元数のまま走らせる。control＝走行1 の床を 1–i 面に置いたもの（走行1 と一致するはず）、main＝q＝e^(iφ0) i^v r j^w。
set -e
cd "$(dirname "$0")"
PY="$HOME/.venvs/stage123_control_py39/bin/python3"
NS="3 4 5 6 8 10 12 15 20 25"
S1=../走行1Q_段階1_四元数で初期化し辺表を値から読む_20261026
R1=../../分配率Rの読出し_解析床_20261025/走行1_全順序対_M_N2_隣接あり_解析床_20261025
START=$(date +%s)
cmp make_initial_quaternion_v1.py $S1/make_initial_quaternion_v1.py && echo "stage1 decoder unchanged"
"$PY" make_state_quaternion_v2.py $S1/floors_init_quaternion states_main states_control $NS 2>results_main/log_states_stderr.txt | tee results_main/log_states.txt
for RUN in control main; do
  "$PY" run_quaternion_sweep_v1.py states_$RUN 'floor_90deg_all_N{N:05d}_analytic.npz' results_$RUN "RUN1Q stage2 $RUN (2026-10-26): quaternion dynamics K^(c)=A(conj(u_e)u_f)_c, c=i,j,k; den=N and 124" $NS 2>results_$RUN/log_sweep_stderr.txt | tee results_$RUN/log_sweep.txt
  "$PY" compare_curves_two_runs_v1.py results_$RUN/timeseries_64bit_with124_N3_N40.csv $R1/results_all/timeseries_64bit_with124_N3_N40.csv q_$RUN run1 results_$RUN/compare_vs_run1 $NS | tee results_$RUN/log_compare.txt
  "$PY" plot_hperp_grid_from_timeseries_v1.py results_$RUN "RUN1Q stage2 $RUN (2026-10-26): quaternion dynamics, Hperp/H outside span of initial 4 columns; den=N and 124" | tee results_$RUN/log_grid.txt
  "$PY" plot_component_dormant_v1.py results_$RUN "RUN1Q stage2 $RUN: per-component part outside Pi" | tee results_$RUN/log_comp.txt
  for P in i j k; do "$PY" plot_quaternion_planes_v1.py results_$RUN results_$RUN $P RUN1Q_stage2_$RUN $NS | tee results_$RUN/log_plot_$P.txt; done
done
echo "elapsed: $(( $(date +%s) - START )) s"
echo 'RUN1Q STAGE2 ALL DONE'
