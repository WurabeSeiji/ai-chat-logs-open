#!/bin/sh
# 関係配列写像 v3：対照走行（v1 と一致確認）→ 決定済みの引き込み（lock on）の走行 → 要約・図・SHA
set -e
cd "$(dirname "$0")"
mkdir -p results
V1=../program_relation_array_v1_20261004/results

# ---- 対照走行：--lock off が v1 と同一であること（先頭 20 列の文字列一致）
python3 relation_array_node_map_v3.py --delta 0.0 --f 1.0 --nodes 20000 --lock off --out results/control_circ_f1.csv | tee results/control_circ_f1.log
python3 compare_control.py "$V1/run_circ_f1.csv" results/control_circ_f1.csv | tee results/control_compare.txt
python3 relation_array_node_map_v3.py --delta 0.0 --f 0.0 --nodes 2000 --lock off --out results/control_circ_f0.csv | tee results/control_circ_f0.log
python3 compare_control.py "$V1/run_circ_f0.csv" results/control_circ_f0.csv | tee -a results/control_compare.txt

# ---- 決定済みの引き込み
python3 relation_array_node_map_v3.py --delta 0.0   --f 0.0 --nodes 2000   --lock on --out results/lock_circ_f0.csv  | tee results/lock_circ_f0.log
python3 relation_array_node_map_v3.py --delta 0.0   --f 1.0 --nodes 100000 --lock on --out results/lock_circ_f1.csv  | tee results/lock_circ_f1.log
python3 relation_array_node_map_v3.py --delta 0.065 --f 1.0 --nodes 100000 --lock on --out results/lock_ecc_f1.csv   | tee results/lock_ecc_f1.log

python3 summarize_v3.py > results/SUMMARY.md
python3 plot_v3.py
diff -u ../program_relation_array_v1_20261004/relation_array_node_map_v1.py relation_array_node_map_v3.py > results/diff_v1_v3.txt || true
gzip -9 -f results/*.csv
shasum -a 256 relation_array_node_map_v3.py summarize_v3.py compare_control.py plot_v3.py run_all.sh README.md results/*.csv.gz results/*.log results/*.txt results/SUMMARY.md results/*.svg results/*.png > SHA256SUMS.txt
