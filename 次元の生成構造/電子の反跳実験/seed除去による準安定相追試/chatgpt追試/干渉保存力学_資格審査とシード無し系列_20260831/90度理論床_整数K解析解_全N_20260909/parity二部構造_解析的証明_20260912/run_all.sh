#!/usr/bin/env bash
# 90度整数K床の parity 二部構造 解析的一本化 検証（読出しのみ・物理無変更）
set -euo pipefail
cd "$(dirname "$0")"
python3 verify_parity_bipartite_analytic_20260912.py
echo "--- SHA256SUMS ---"
shasum -a 256 verify_parity_bipartite_analytic_20260912.py \
  分析_parity二部構造_解析的一本化_20260912.md README.md run_all.sh \
  results/parity_bipartite_summary.csv results/実行ログ_20260912.log > SHA256SUMS.txt
cat SHA256SUMS.txt
