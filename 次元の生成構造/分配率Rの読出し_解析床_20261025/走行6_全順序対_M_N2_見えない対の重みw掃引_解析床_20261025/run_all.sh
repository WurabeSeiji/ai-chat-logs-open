#!/bin/bash
# 走行6（解析床、2026-10-25）: 見えない対の重み w を掃引。w ∈ {0.5, 0.7, 0.9, 0.95, 0.99}。プログラムは走行5 の同一コピー。集計で w=0（走行5）・0.3（走行5）・1.0（走行2）も並べる。
set -e
cd "$(dirname "$0")"
PY="$HOME/.venvs/stage123_control_py39/bin/python3"
NS="3 4 5 6 8 10 12 15"
START=$(date +%s)
for W in 0.5 0.7 0.9 0.95 0.99; do
  mkdir -p floors_w$W results_w$W
  "$PY" floor_analytic_v5_w.py all floors_w$W $W $NS 2>results_w$W/log_floors_stderr.txt | tee results_w$W/log_floors.txt
  "$PY" run_sweep_from_parents_v5_w.py floors_w$W 'floor_90deg_all_N{N:05d}_analytic.npz' results_w$W "RUN6 analytic (2026-10-25): analytic floor on all ordered pairs (M=N^2), invisible pairs weighted w=$W, one_step unchanged, seedless; den=N and 124" $W $NS 2>results_w$W/log_sweep_stderr.txt | tee results_w$W/log_sweep.txt
done
"$PY" summarize_w_sweep_v1.py | tee results_summary/log_summary.txt
echo "elapsed: $(( $(date +%s) - START )) s"
echo 'RUN6 ALL DONE'
