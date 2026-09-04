#!/bin/bash
# 段2除去版（振幅込み・虚部のみ生成子）の一括ラッパー
set -e
cd "$(dirname "$0")"
python3 run_N40_staticparent_ampimK_v1.py
python3 check_staticparent_ampimK_inputs_v1.py
python3 plot_complex_plane_N40_staticparent_ampimK_v1.py
python3 plot_inflation_N40_ampimK_v1.py
echo AMPIMK ALL DONE
