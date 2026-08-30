#!/bin/bash
# 全パス再実行（親生成・予測は走行前に固定される順序）。全体 ~12 分。
set -e; set -o pipefail
cd "$(dirname "$0")"
python3 program/pass1_family.py            > results/pass1.log
python3 program/pass2_jacobian.py          > results/pass2.log
python3 program/pass2b_nullspace_probe.py  > results/pass2b.log
python3 program/pass3_parents.py           > results/pass3.log          # 親＋走行前予測（parents_predictions.csv）
python3 program/pass7_balanced_random.py   > results/pass7.log          # 恒等式検証＋乱数均衡親＋走行前予測（balanced_random_parents.csv）
for t in N6_eps0.00_k2 N6_eps0.30_k2 N6_eps0.60_k2 N6_eps0.90_k2 N8_eps0.00_k2 N8_eps0.60_k2 N8_eps0.60_k3 N9_eps0.00_k2 N9_eps0.60_k2; do
  python3 program/run_dynamics.py $t; done                 > results/run_dynamics.log
for N in 5 6 7 8; do for s in 0 1 2 3 4; do python3 program/run_dynamics.py random_N${N}_s$s; done; done > results/run_dynamics_random.log
gzip -f data/*/treatment_linear124_amplitude_aware_timeseries.csv                 # 走行 CSV は gzip 保存（読出しは両対応）
python3 program/pass3b_monodromy_calibration.py > results/pass3b.log   # 走行後に追加した較正（本文 §7.3）
python3 program/pass5_analysis.py          > results/pass5.log
python3 program/pass6_figures.py
bash build_pdf.sh
find . -type f ! -name SHA256SUMS.txt -print0 | sort -z | xargs -0 shasum -a 256 > SHA256SUMS.txt
echo "done"
