#!/bin/sh
# 八状態の指数表現・固定交換係数 θ₀ = 23π/124（固定点 5）：対照 → 走行 → 要約・図・SHA
set -e
cd "$(dirname "$0")"
mkdir -p results
P=eight_state_exchange.py

# ---- 対照：Kepler のみ、放射なし。閉じた式と比較（1240 節点 = 10 再帰 = 230 周）
python3 $P --em kepler --f 0.0 --delta 0.0   --nodes 1240 --out results/run_kepler_circ_f0.csv | tee results/run_kepler_circ_f0.log
python3 $P --em kepler --f 0.0 --delta 0.065 --nodes 1240 --out results/run_kepler_ecc_f0.csv  | tee results/run_kepler_ecc_f0.log
{ echo "== circ =="; python3 kepler_reference.py results/run_kepler_circ_f0.csv; echo "== ecc =="; python3 kepler_reference.py results/run_kepler_ecc_f0.csv; } | tee results/kepler_check.txt

# ---- E–M 全項（1PN〜は節点ごとの積和のキック）、放射なし：近点移動を E–M の解析値と比較（12400 節点 = 2300 周）
python3 $P --em full --f 0.0 --delta 0.0   --nodes 12400 --out results/run_full_circ_f0.csv | tee results/run_full_circ_f0.log
python3 $P --em full --f 0.0 --delta 0.01  --nodes 12400 --out results/run_full_e002_f0.csv | tee results/run_full_e002_f0.log
python3 $P --em full --f 0.0 --delta 0.065 --nodes 12400 --out results/run_full_ecc_f0.csv  | tee results/run_full_ecc_f0.log
python3 precession_check.py results/run_full_circ_f0.csv results/run_full_e002_f0.csv results/run_full_ecc_f0.csv | tee results/precession_check.txt

# ---- 放射あり（B3、片側）。因子は掛けない（読み (a)(b) は出力のみ）。124000 節点 = 1000 再帰 = 23000 周
python3 $P --em full --f 1.0 --delta 0.0   --nodes 124000 --out results/run_full_circ_f1.csv | tee results/run_full_circ_f1.log
python3 $P --em full --f 1.0 --delta 0.065 --nodes 124000 --out results/run_full_ecc_f1.csv  | tee results/run_full_ecc_f1.log

python3 summarize.py > results/SUMMARY.md
python3 plot.py
gzip -9 -f results/*.csv
shasum -a 256 $P kepler_reference.py precession_check.py summarize.py plot.py run_all.sh README.md results/*.csv.gz results/*.log results/*.txt results/SUMMARY.md results/*.svg results/*.png > SHA256SUMS.txt
