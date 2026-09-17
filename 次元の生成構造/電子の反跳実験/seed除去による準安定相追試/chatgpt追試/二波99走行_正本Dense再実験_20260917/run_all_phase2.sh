#!/bin/zsh
# Phase 2 実行手順（最小系列 L=8, den=8, 4096 steps × 3条件）— 承認後にのみ実行する。
# Phase 1（コード監査）の時点では本スクリプトは実行しない。
set -e
cd "$(dirname "$0")"

# 0) 事前登録 manifest 生成（走行行列の唯一の定義源。走行後に変更しない）
python3 gen_manifest.py

# 1) §4 構造診断（生値記録のみ。判定・ゲートなし — 監査するのはコード、判定するのは実験後のデータ）
python3 structural_audit_v1.py L8_ma1_mb2_den8 L8_ma1_mb3_den8 L8_ma2_mb4_den8

# 診断値は pre_run_structural_audit.json に記録される（停止条件ではない）

# 2) 最小系列 第1走行（奇偶 (1,2) のみ。結果と図化の確認後に残り2本）
# python3 run_dense_rerun.py L8_ma1_mb2_den8

# 3) 確認後、残り2本
# python3 run_dense_rerun.py L8_ma1_mb3_den8 L8_ma2_mb4_den8

# 4) 図化 — plot_scripts_run/ は「存在する run だけ図化（欠損パネルは空欄）」対応済み
#    （木原承認 2026-09-17、README.md 参照）。走行数に関わらずそのまま実行できる。
# cd plot_scripts_run
# python3 plot_initial_states_all99_v4_20260916.py  > ../plot_logs/initial.log 2>&1
# python3 plot_final_states_all99_20260916.py       > ../plot_logs/final.log 2>&1
# python3 plot_inflation_hperp_all99_20260916.py    > ../plot_logs/inflation.log 2>&1
