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
shasum -a 256 program/*.py run_all.sh README.md *.md results/*.csv results/*.md results/*.json > SHA256SUMS.txt
echo ALL DONE
