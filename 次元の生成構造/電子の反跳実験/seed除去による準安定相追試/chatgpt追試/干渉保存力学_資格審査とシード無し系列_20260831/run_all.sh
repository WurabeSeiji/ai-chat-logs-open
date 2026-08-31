#!/bin/bash
# 干渉保存力学：資格審査 → 親生成＋走行前予測固定 → シード無し走行 40 本 → 判定 → 図。全再現 約 15 分。
set -e
cd "$(dirname "$0")"
python3 program/pass0_qualification.py
python3 program/pass1_parents.py
for d in data/*/; do
  tag=$(basename "$d")
  python3 program/pass2_run.py "$tag"
done
python3 program/pass3_analysis.py
python3 program/pass4_figures.py
shasum -a 256 program/*.py run_all.sh README.md results/*.csv results/*.md results/*.json > SHA256SUMS.txt
echo ALL DONE
