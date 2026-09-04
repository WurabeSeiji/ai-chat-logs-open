#!/bin/bash
# 静的親差し替え実験の一括ラッパー（走行→入力ゲート→図化）
set -e
cd "$(dirname "$0")"
python3 run_N40_staticparent_v1.py
python3 check_staticparent_inputs_v1.py
python3 plot_complex_plane_N40_staticparent_v1.py
echo 'STATICPARENT ALL DONE'
