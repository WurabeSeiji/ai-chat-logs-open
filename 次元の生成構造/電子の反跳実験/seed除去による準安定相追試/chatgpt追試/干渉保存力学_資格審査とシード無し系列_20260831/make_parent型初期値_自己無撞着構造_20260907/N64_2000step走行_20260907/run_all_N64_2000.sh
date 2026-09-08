#!/bin/zsh
set -e
cd "$(dirname "$0")"
python3 -u wrapper_run_N64_2000_v1.py
python3 plot_inflation_N64_v1.py
python3 plot_complex_plane_final_N64_v1.py
python3 plot3d_N64_2000_v1.py
echo 'N64 2000STEP ALL DONE'
