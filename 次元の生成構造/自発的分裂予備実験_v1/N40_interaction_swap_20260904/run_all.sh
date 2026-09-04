#!/bin/bash
# N=40 相互作用スワップ実験の再実行（走行→図化）
set -e
cd "$(dirname "$0")"
python3 run_spontaneous_splitting_largeN_v1_interactionswap.py 40 1e-15 --after=1500 --tol=1e-12
python3 plot_swap_results_v1.py
