#!/usr/bin/env bash
# 縮退真空の頂点置換軌道テスト N=6（読出しのみ・物理無変更）
set -euo pipefail
cd "$(dirname "$0")"
python3 verify_vacua_permutation_orbit_N6_20260913.py
echo "--- SHA256SUMS ---"
shasum -a 256 verify_vacua_permutation_orbit_N6_20260913.py \
  分析_縮退真空の頂点置換軌道テスト_N6_20260913.md README.md run_all.sh \
  results/実行ログ_20260913.log > SHA256SUMS.txt
cat SHA256SUMS.txt
