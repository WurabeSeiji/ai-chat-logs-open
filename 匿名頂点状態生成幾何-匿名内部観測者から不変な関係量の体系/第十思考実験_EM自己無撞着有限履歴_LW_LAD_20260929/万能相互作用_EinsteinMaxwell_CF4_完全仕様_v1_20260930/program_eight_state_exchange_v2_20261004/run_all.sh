#!/bin/sh
# 八状態 v2（時計の文字盤を状態に、T は状態から読む、放射は乗法更新）：対照 → 走行 → 要約・図・SHA
set -e
cd "$(dirname "$0")"
mkdir -p results
P=eight_state_exchange_v2.py

python3 $P --em kepler --f 0.0 --delta 0.0   --nodes 1240 --out results/run_kepler_circ_f0.csv | tee results/run_kepler_circ_f0.log
python3 $P --em kepler --f 0.0 --delta 0.065 --nodes 1240 --out results/run_kepler_ecc_f0.csv  | tee results/run_kepler_ecc_f0.log
{ echo "== circ =="; python3 kepler_reference.py results/run_kepler_circ_f0.csv; echo "== ecc =="; python3 kepler_reference.py results/run_kepler_ecc_f0.csv; } | tee results/kepler_check.txt

python3 $P --em full --f 0.0 --delta 0.0   --nodes 12400 --out results/run_full_circ_f0.csv | tee results/run_full_circ_f0.log
python3 $P --em full --f 0.0 --delta 0.065 --nodes 12400 --out results/run_full_ecc_f0.csv  | tee results/run_full_ecc_f0.log

python3 $P --em full --f 1.0 --delta 0.0   --nodes 124000 --out results/run_full_circ_f1.csv | tee results/run_full_circ_f1.log
python3 $P --em full --f 1.0 --delta 0.065 --nodes 124000 --out results/run_full_ecc_f1.csv  | tee results/run_full_ecc_f1.log

python3 summarize.py > results/SUMMARY.md
python3 plot.py
gzip -9 -f results/*.csv
shasum -a 256 $P kepler_reference.py summarize.py plot.py run_all.sh README.md results/*.csv.gz results/*.log results/*.txt results/SUMMARY.md results/*.svg results/*.png > SHA256SUMS.txt
