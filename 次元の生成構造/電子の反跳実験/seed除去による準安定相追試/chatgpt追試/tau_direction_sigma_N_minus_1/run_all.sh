#!/bin/bash
# 全再現: python3 + numpy のみ。各スクリプトの標準出力を results/*.log に保存。
set -e; cd "$(dirname "$0")"; mkdir -p results
export LC_ALL=en_US.UTF-8
python3 analyze_n5_recurrence.py 2>&1 | grep -v Warning | tee results/analyze_n5_recurrence.log
python3 sweep_sigma2.py 5000 3,4,5,6,7,8,9,10,11,12,13,14,15,16 2>&1 | grep -v Warning | tee results/sweep_sigma2.log
python3 sweep_floquet.py 3,4,5,6,7,8,9,10,11,12,13,14,15,16 2>&1 | grep -v Warning | tee results/sweep_floquet.log
python3 verify_theorem.py 2>&1 | grep -v Warning | tee results/verify_theorem.log
