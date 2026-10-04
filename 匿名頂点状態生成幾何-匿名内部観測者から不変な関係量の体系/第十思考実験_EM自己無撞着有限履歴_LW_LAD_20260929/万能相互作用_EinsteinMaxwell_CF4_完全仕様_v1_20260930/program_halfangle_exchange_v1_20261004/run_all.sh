#!/bin/sh
# 半角位相子の交換（固定点 4）：対照（Kepler 閉じた式）→ 走行 → 要約・図・差分・SHA
set -e
cd "$(dirname "$0")"
mkdir -p results
P=halfangle_exchange_node_map.py

# ---- 対照：Kepler のみ、放射なし。閉じた式と比較
python3 $P --em kepler --f 0.0 --delta 0.0   --nodes 2000 --out results/run_kepler_circ_f0.csv | tee results/run_kepler_circ_f0.log
python3 $P --em kepler --f 0.0 --delta 0.065 --nodes 2000 --out results/run_kepler_ecc_f0.csv  | tee results/run_kepler_ecc_f0.log
{ echo "== circ =="; python3 kepler_reference.py results/run_kepler_circ_f0.csv; echo "== ecc =="; python3 kepler_reference.py results/run_kepler_ecc_f0.csv; } | tee results/kepler_check.txt

# ---- 対照：近点移動の一次摂動を v7 §11.1 の表と比較
python3 precession_check.py | tee results/precession_check.txt

# ---- E–M 全項（近点移動は一次摂動）、放射なし：近点移動の読み
python3 $P --em full --f 0.0 --delta 0.0   --nodes 2000 --out results/run_full_circ_f0.csv | tee results/run_full_circ_f0.log
python3 $P --em full --f 0.0 --delta 0.065 --nodes 2000 --out results/run_full_ecc_f0.csv  | tee results/run_full_ecc_f0.log

# ---- 放射あり（背景への片側の交換）。ロックの因子は掛けない（読み (a)(b) は出力のみ）
python3 $P --em full --f 1.0 --delta 0.0   --nodes 100000 --out results/run_full_circ_f1.csv | tee results/run_full_circ_f1.log
python3 $P --em full --f 1.0 --delta 0.065 --nodes 100000 --out results/run_full_ecc_f1.csv  | tee results/run_full_ecc_f1.log

python3 summarize.py > results/SUMMARY.md
python3 plot.py
gzip -9 -f results/*.csv
shasum -a 256 $P kepler_reference.py precession_check.py summarize.py plot.py run_all.sh README.md results/*.csv.gz results/*.log results/*.txt results/SUMMARY.md results/*.svg results/*.png > SHA256SUMS.txt
