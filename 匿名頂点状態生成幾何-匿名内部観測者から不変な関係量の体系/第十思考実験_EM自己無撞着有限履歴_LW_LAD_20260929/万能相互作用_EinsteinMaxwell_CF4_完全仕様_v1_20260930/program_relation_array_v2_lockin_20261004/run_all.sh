#!/bin/sh
# 関係配列写像 v2：対照走行（v1 と一致確認）→ 引き込み A/B × 有限履歴 N の走行 → 要約・図・SHA
set -e
cd "$(dirname "$0")"
mkdir -p results
V1=../program_relation_array_v1_20261004/results

# ---- 対照走行：--lock none が v1 と同一であること（先頭 20 列の文字列一致）
python3 relation_array_node_map_v2.py --delta 0.0 --f 1.0 --nodes 20000 --lock none --out results/control_circ_f1.csv | tee results/control_circ_f1.log
python3 compare_control.py "$V1/run_circ_f1.csv" results/control_circ_f1.csv | tee results/control_compare.txt
python3 relation_array_node_map_v2.py --delta 0.0 --f 0.0 --nodes 2000 --lock none --out results/control_circ_f0.csv | tee results/control_circ_f0.log
python3 compare_control.py "$V1/run_circ_f0.csv" results/control_circ_f0.csv | tee -a results/control_compare.txt

# ---- 実験：引き込み A（時計側）/ B（軌道側）× 有限履歴 N
for N in 62 93 124 248; do
  for L in A B; do
    python3 relation_array_node_map_v2.py --delta 0.0 --f 1.0 --nodes 50000 --lock $L --hist $N --out results/lock_${L}_N${N}.csv | tee results/lock_${L}_N${N}.log
  done
done

python3 summarize_v2.py > results/SUMMARY.md
python3 plot_lockin_v2.py
# csv は大きい（8 走行 × 50000 節点で約 150 MB）ので gzip して保存する。読むときは gunzip。
gzip -9 -f results/*.csv
shasum -a 256 relation_array_node_map_v2.py summarize_v2.py compare_control.py plot_lockin_v2.py run_all.sh README.md results/*.csv.gz results/*.log results/*.txt results/SUMMARY.md results/*.svg results/*.png > SHA256SUMS.txt
