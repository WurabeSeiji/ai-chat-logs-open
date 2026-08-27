#!/bin/bash
# 全再現: python3 + numpy + matplotlib。参照パッケージ ../N5_linear124_all3fix_seedless_parentnorm_removed_20260828/data が必要。
set -e; cd "$(dirname "$0")"; mkdir -p results data figures; export LC_ALL=en_US.UTF-8
python3 program/run_amplitude_only_fix.py > results/rerun_stdout.log 2> results/rerun_progress.log   # 参照プログラムを無変更で再実行 → data/, figures/
python3 program/compare_with_reference.py 2>&1 | grep -v Warning | tee results/compare_with_reference.log
python3 program/scale_and_engine_check.py 2>&1 | grep -v Warning | tee results/scale_and_engine_check.log
python3 program/parent_symmetry_check.py 2>&1 | grep -v Warning | tee results/parent_symmetry_check.log
