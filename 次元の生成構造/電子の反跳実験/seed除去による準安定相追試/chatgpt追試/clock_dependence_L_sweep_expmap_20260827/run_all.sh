#!/bin/bash
# 全再現: python3 + numpy のみ。標準出力を results/*.log に保存。
set -e; cd "$(dirname "$0")"; mkdir -p results; export LC_ALL=en_US.UTF-8
python3 sweep_L_and_expmap.py 2>&1 | grep -v Warning | tee results/sweep_L_and_expmap.log
python3 ordering_vs_L.py 2>&1 | grep -v Warning | tee results/ordering_vs_L.log
python3 floquet_vs_L.py 2>&1 | grep -v Warning | tee results/floquet_vs_L.log
