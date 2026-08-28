#!/bin/bash
# 全再現: python3 + numpy + matplotlib。参照 ../N5_linear124_all3fix_seedless_parentnorm_removed_40000_20260828/data が必要。約 1 分。
set -e; cd "$(dirname "$0")"; mkdir -p results data figures; export LC_ALL=en_US.UTF-8
python3 program/run_amplitude_only_fix.py > results/rerun_stdout.log 2> results/rerun_progress.log   # 参照プログラムを無変更で再実行（40000 step）
python3 program/compare_with_reference.py 2>&1 | grep -v Warning | tee results/compare_with_reference.log
python3 program/verify_fixes.py 2>&1 | grep -v Warning | tee results/verify_fixes.log
python3 program/divergence_onset.py | tee results/divergence_onset.log
