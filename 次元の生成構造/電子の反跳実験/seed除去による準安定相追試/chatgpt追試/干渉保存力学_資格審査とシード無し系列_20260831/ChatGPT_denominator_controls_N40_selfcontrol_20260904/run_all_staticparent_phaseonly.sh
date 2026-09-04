#!/bin/bash
# 位相正規化版（振幅正規化の追加のみ）の一括ラッパー
set -e
cd "$(dirname "$0")"
python3 run_N40_staticparent_phaseonly_v1.py
python3 check_staticparent_phaseonly_inputs_v1.py
python3 plot_complex_plane_N40_staticparent_phaseonly_v1.py
echo PHASEONLY ALL DONE
