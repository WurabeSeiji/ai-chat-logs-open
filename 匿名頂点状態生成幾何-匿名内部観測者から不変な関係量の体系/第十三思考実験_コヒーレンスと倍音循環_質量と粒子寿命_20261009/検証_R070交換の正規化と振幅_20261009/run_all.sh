#!/bin/bash
# 0713 原本の忠実コピーで対照走行（約 5.5 分）→ 正規化の記録と純粋な U_R との比較
set -e
cd "$(dirname "$0")"
( cd original_copy && python3 run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1.py > control_stdout.txt 2>&1 )
python3 measure_normalization_trace.py > results/measure_stdout.txt
find . -name ".matplotlib" -type d -prune -exec rm -rf {} +
find . -name "__pycache__" -type d -prune -exec rm -rf {} +
