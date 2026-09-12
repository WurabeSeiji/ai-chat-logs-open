#!/usr/bin/env bash
# 巡回二乗波の因数分解定理 検証（読出しのみ・物理無変更）
set -euo pipefail
cd "$(dirname "$0")"
python3 verify_cyclic_square_wave_factorization_20260913.py
echo "--- SHA256SUMS ---"
shasum -a 256 verify_cyclic_square_wave_factorization_20260913.py \
  分析_巡回二乗波因数分解定理_20260913.md README.md run_all.sh \
  results/cyclic_factorization_summary.csv results/実行ログ_20260913.log > SHA256SUMS.txt
cat SHA256SUMS.txt
