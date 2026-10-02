#!/bin/bash
# v7（トーラス波 Ψ_{n,m}）の動作確認走行一式。仕様: ../coefficients_relational_full_list_ja_final_v7_20261003.md
# 表の構築（初回 45 秒、以後キャッシュ）＋ 4 走行 ＋ 図 ＋ SHA256。本走行は別途指示で。
set -e
cd "$(dirname "$0")"
mkdir -p results_v7
PY="python3 -W error::RuntimeWarning em_two_body_torus_node_map.py"
echo "== selftest（Newton のみの表 vs Kepler の閉じた式）"
$PY --selftest | tee results_v7/selftest.txt
echo "== t1: δ=0, f=0, 4 呼吸（単一 (0,1)：D=0、W 不変）"
$PY --delta 0 --f 0 --breaths 4 --out results_v7/t1_delta0_f0.csv | tee results_v7/t1.txt
echo "== t2: δ=0.065, f=0, 100 呼吸（保存：W と ε̄ が厳密に不変）"
$PY --delta 0.065 --f 0 --breaths 100 --out results_v7/t2_delta065_f0.csv | tee results_v7/t2.txt
echo "== t3: δ=0.065, f=1, 1000 呼吸（放射：w(0,2) 減少、w(1,1) 不変）"
$PY --delta 0.065 --f 1 --breaths 1000 --out results_v7/t3_delta065_f1.csv | tee results_v7/t3.txt
echo "== t4: δ=0.065, f=1, 20000 呼吸（e 折り時間の推定、約 4 分）"
$PY --delta 0.065 --f 1 --breaths 20000 --log-every 62 --out results_v7/t4_delta065_f1_20000.csv | tee results_v7/t4.txt
echo "== figures"
python3 plot_readouts_v7.py results_v7/t2_delta065_f0.csv results_v7/t3_delta065_f1.csv results_v7/t4_delta065_f1_20000.csv --out results_v7/figures > /dev/null
echo "== SHA256"
shasum -a 256 em_two_body_node_map.py em_two_body_torus_node_map.py checks.py plot_readouts.py plot_readouts_v7.py run_all.sh run_all_v7.sh README.md \
  results/*.csv results/*.json results/*.txt results/figures/* results_B3/* results_v7/*.csv results_v7/*.json results_v7/*.txt results_v7/figures/* > SHA256SUMS.txt
echo done
