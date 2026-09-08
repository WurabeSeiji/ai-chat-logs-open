#!/bin/zsh
# N=100・den=100・2000歩の一括再現（走行→図化3本）。実行は木原指示後。
set -e
cd "$(dirname "$0")"
python3 -u wrapper_run_N100_2000_v1.py
python3 plot_inflation_N100_v1.py
python3 plot_complex_plane_final_N100_v1.py
python3 plot3d_N100_2000_v1.py
echo 'N100 2000STEP ALL DONE'
