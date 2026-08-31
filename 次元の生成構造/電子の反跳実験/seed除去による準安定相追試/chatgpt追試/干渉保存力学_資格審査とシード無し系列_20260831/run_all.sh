#!/bin/bash
# 干渉保存力学：v2補完実験と同一プロトコル。資格審査 → 親54＋走行前予測固定 → 埋め込み → 走行54×40000 → 集計 → 図 → 星構造。
set -e
cd "$(dirname "$0")"
python3 program/pass0_qualification.py
python3 program/pass1_parents.py
python3 program/pass2_embed_random.py
for d in data/*/; do
  tag=$(basename "$d")
  [ "$tag" = "reference" ] && continue
  python3 program/pass2_run.py "$tag"
done
python3 program/pass5_analysis.py
python3 program/pass6_figures.py
python3 program/pass7_final_structure.py
# 初期拡大診断（L=12400・500 step、追加。木原指示 2026-08-31）
for d in data/*/; do
  tag=$(basename "$d")
  [ "$tag" = "reference" ] && continue
  python3 program/pass2c_run_L12400.py "$tag"
done
python3 program/pass6c_figures_L12400.py
# 初期拡大診断第 2 段（L=124000000・500 step。木原指示 2026-08-31、L=1240000 段は破棄）
for d in data/*/; do
  tag=$(basename "$d")
  [ "$tag" = "reference" ] && continue
  python3 program/pass2d_run_L124000000.py "$tag"
done
python3 program/pass6d_figures_L124000000.py
# 読出し・分析パス（パス8〜15）
python3 program/pass8_ladder_geometry.py
python3 program/pass9_composite_wave.py
python3 program/pass10_harmonic_ladder.py
python3 program/pass11_edge_wavelengths.py
python3 program/pass12_wavelength_table.py
python3 program/pass13_k_enumeration.py
python3 program/pass14_k_search_star.py
python3 program/pass15_verification_records.py
python3 program/pass16_k_search_ne_N3.py
python3 program/pass17_hm_series_k.py
shasum -a 256 program/*.py run_all.sh README.md *.md results/*.csv results/*.md results/*.json > SHA256SUMS.txt
echo ALL DONE
