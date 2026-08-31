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
shasum -a 256 program/*.py run_all.sh README.md *.md results/*.csv results/*.md results/*.json > SHA256SUMS.txt
echo ALL DONE
