#!/bin/sh
# 関係配列写像 v1：三走行（円・離心・円＋放射）。本走行ではない（短い）。
set -e
cd "$(dirname "$0")"
mkdir -p results
python3 relation_array_node_map_v1.py --delta 0.0   --f 0.0 --nodes 2000  --out results/run_circ_f0.csv   | tee results/run_circ_f0.log
python3 relation_array_node_map_v1.py --delta 0.065 --f 0.0 --nodes 2000  --out results/run_ecc_f0.csv    | tee results/run_ecc_f0.log
python3 relation_array_node_map_v1.py --delta 0.0   --f 1.0 --nodes 20000 --out results/run_circ_f1.csv   | tee results/run_circ_f1.log
python3 summarize_v1.py > results/SUMMARY.md
shasum -a 256 relation_array_node_map_v1.py summarize_v1.py run_all.sh README.md results/*.csv results/*.log results/SUMMARY.md > SHA256SUMS.txt
