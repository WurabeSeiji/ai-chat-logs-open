#!/usr/bin/env bash
# 90度理論床（整数K σ_max 解析解）一式の再現（2026-09-09）。
set -e
cd "$(dirname "$0")"
python3 make_90deg_floor_analytic_v1.py     2>&1 | tee 実行ログ_make_90deg_floor_analytic_20260909.log
python3 closed_form_90deg_floor_verify_v1.py 2>&1 | tee 実行ログ_closed_form_verify_20260909.log
echo "ALL DONE"
