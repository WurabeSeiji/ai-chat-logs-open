#!/bin/bash
# 段2+段3+σ時計（段1の固定Δτを除去）の一括ラッパー
set -e
cd "$(dirname "$0")"
python3 run_N40_staticparent_sigmaclock_v1.py
python3 check_staticparent_sigmaclock_inputs_v1.py
python3 plot_complex_plane_N40_staticparent_sigmaclock_v1.py
python3 plot_inflation_N40_sigmaclock_v1.py
echo SIGMACLOCK ALL DONE
