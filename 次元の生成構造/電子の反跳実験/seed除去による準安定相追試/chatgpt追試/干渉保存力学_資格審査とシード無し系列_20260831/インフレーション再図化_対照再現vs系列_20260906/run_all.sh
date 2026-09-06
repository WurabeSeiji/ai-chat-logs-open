#!/bin/zsh
# インフレーション再図化: 論文A対照再現（δ=0.01 / δ=0）vs 現行シード無し系列（読み出しのみ・新規走行なし）
cd "$(dirname "$0")"
python3 plot_inflation_control_vs_seedless_v1.py
python3 plot_inflation_delta_sweep_v1.py
