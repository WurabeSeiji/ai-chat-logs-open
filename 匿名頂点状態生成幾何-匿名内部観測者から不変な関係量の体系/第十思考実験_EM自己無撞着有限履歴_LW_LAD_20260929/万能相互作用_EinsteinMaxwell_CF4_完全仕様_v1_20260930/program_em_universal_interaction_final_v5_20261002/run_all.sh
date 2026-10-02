#!/bin/bash
# 動作確認走行（短時間。本走行 B3 ≈ 2×10^5 周は別途指示で）
# 仕様: ../coefficients_relational_full_list_ja_final_v5_20261002.md
set -e
cd "$(dirname "$0")"
mkdir -p results
PY="python3 -W error::RuntimeWarning em_two_body_node_map.py"

echo "== selftest（単項式表 = T1..T22、exp 級数）"
$PY --selftest | tee results/selftest.txt

echo "== test1a: δ=0, f=0, スピン 0, 4 周（§3・§4 の数値の再現：γ−1、a サイクル/周、ψ1 位相/周）"
$PY --delta 0 --f 0 --sa zero --sb zero --orbits 4 --out results/test1a_delta0_f0_spin0.csv | tee results/test1a.txt

echo "== test1b: δ=0, f=0, スピン up/down, 4 周"
$PY --delta 0 --f 0 --sa up --sb down --orbits 4 --out results/test1b_delta0_f0_updown.csv | tee results/test1b.txt

echo "== test2: δ=0.065, f=0, 呼吸 100 周（保存：同じ呼吸位相で ε̄・A・⟨m²⟩ が一定）"
$PY --delta 0.065 --f 0 --orbits 100 --log-every 62 --out results/test2_delta065_f0.csv | tee results/test2.txt

echo "== test3s: δ=0.065, f=1, 呼吸 100 周, 全節点記録（エネルギー収支：Δε̄ = ∫P_rad dθ）"
$PY --delta 0.065 --f 1 --orbits 100 --log-every 1 --out results/test3s_delta065_f1_every_node.csv | tee results/test3s.txt

echo "== test3: δ=0.065, f=1, 呼吸 1000 周（放射反作用の方向：⟨m²⟩ 減少、ε̄ 減少）"
$PY --delta 0.065 --f 1 --orbits 1000 --log-every 62 --out results/test3_delta065_f1.csv | tee results/test3.txt

echo "== checks"
python3 checks.py | tee results/checks.txt

echo "== SHA256"
shasum -a 256 em_two_body_node_map.py checks.py plot_readouts.py run_all.sh README.md results/*.csv results/*.json results/*.txt results/figures/* results_B3/* > SHA256SUMS.txt
echo done
