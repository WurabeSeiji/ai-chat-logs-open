#!/bin/bash
# 段3（虚部のみ生成子 iK・実直交回転）の一括ラッパー
set -e
cd "$(dirname "$0")"
python3 run_N40_staticparent_imK_v1.py
python3 check_staticparent_imK_inputs_v1.py
python3 plot_complex_plane_N40_staticparent_imK_v1.py
python3 plot_inflation_N40_imK_v1.py
echo IMK ALL DONE
