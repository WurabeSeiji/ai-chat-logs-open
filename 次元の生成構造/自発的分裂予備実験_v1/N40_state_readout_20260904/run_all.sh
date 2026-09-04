#!/bin/bash
# N=40 状態読出し一式の再実行（対照実験→状態保存走行→図化）
set -e
cd "$(dirname "$0")"
python3 run_spontaneous_splitting_largeN_v1.py 40 1e-15 --after=1500 --tol=1e-12
diff largeN_splitting_result_v1/fcurve_N00040_delta1e-15_seed0.csv \
     ../largeN_splitting_result_v1/fcurve_N00040_delta1e-15_seed0.csv && echo FCURVE_IDENTICAL
python3 run_spontaneous_splitting_largeN_v1_savestate.py 40 1e-15 --after=1500 --tol=1e-12
diff largeN_splitting_result_v1/fcurve_N00040_delta1e-15_seed0.csv \
     ../largeN_splitting_result_v1/fcurve_N00040_delta1e-15_seed0.csv && echo FCURVE_STILL_IDENTICAL
python3 plot_complex_plane_N40_v1.py
