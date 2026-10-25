#!/bin/bash
# 走行6 やり直し（2026-10-25、走行系固定）：初期値＝反対辺に重み w で解いた床、w ∈ {0.5,0.7,0.9,0.95,0.99,1.0}。走行系＝走行1 の run_sweep_from_parents_v1.py（バイト同一）。集計で w=0・0.3（走行5′）も並べる。
set -e
cd "$(dirname "$0")"
PY="$HOME/.venvs/stage123_control_py39/bin/python3"
NS="3 4 5 6 8 10 12 15"
R1=../../分配率Rの読出し_解析床_20261025/走行1_全順序対_M_N2_隣接あり_解析床_20261025
START=$(date +%s)
cmp run_sweep_from_parents_v1.py $R1/run_sweep_from_parents_v1.py && echo "走行系 unchanged (byte-identical to RUN1)"
for W in 0.5 0.7 0.9 0.95 0.99 1.0; do
  mkdir -p floors_init_w$W results_w$W
  "$PY" make_initial_floor_w_v1.py all floors_init_w$W $W $NS 2>results_w$W/log_floors_stderr.txt | tee results_w$W/log_floors.txt
  "$PY" run_sweep_from_parents_v1.py floors_init_w$W 'floor_90deg_all_N{N:05d}_analytic.npz' results_w$W "RUN6 redo (2026-10-25): INITIAL VALUE = floor solved with invisible-pair weight w=$W; DYNAMICS = RUN1 program unchanged, seedless; den=N and 124" $NS 2>results_w$W/log_sweep_stderr.txt | tee results_w$W/log_sweep.txt
done
"$PY" summarize_w_sweep_initial_v1.py | tee results_summary/log_summary.txt
echo "elapsed: $(( $(date +%s) - START )) s"
echo 'RUN6 REDO ALL DONE'
