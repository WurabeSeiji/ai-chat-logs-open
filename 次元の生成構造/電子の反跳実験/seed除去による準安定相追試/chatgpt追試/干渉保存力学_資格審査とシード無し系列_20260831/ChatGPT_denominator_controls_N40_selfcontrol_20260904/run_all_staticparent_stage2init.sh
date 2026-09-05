#!/bin/bash
# 段2を初期化時のみ実行する検証実験の一括ラッパー
set -e
cd "$(dirname "$0")"
python3 run_N40_staticparent_stage2init_v1.py
python3 check_staticparent_stage2init_inputs_v1.py
python3 plot_complex_plane_N40_staticparent_stage2init_v1.py
python3 plot_inflation_N40_stage2init_v1.py
echo STAGE2INIT ALL DONE
