#!/usr/bin/env bash
# 加法法則【90度解析床】実験一式の再現（2026-09-09）。実行順: 正規化親 → N=9単独 → N=7+6混合 → 部分系別図。
set -e
cd "$(dirname "$0")"
python3 make_normalized_parents_N7N6_90degAnalytic_v1.py 2>&1 | tee 実行ログ_make_normalized_parents_90degAnalytic_20260909.log
python3 run_single_simplex_N9_90degAnalytic_v1.py        2>&1 | tee 実行ログ_single_simplex_N9_90degAnalytic_20260909.log
python3 run_mixed_simplex_N7N6_90degAnalytic_v1.py       2>&1 | tee 実行ログ_mixed_simplex_N7N6_90degAnalytic_20260909.log
python3 plot_Hperp_mixed_parts_N7N6_90degAnalytic_v1.py  2>&1 | tee 実行ログ_plot_Hperp_mixed_parts_90degAnalytic_20260909.log
echo "ALL EXPERIMENTS DONE"
