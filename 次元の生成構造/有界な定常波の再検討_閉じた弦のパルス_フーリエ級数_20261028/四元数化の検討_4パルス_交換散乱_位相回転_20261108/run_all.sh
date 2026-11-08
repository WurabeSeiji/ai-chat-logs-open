#!/bin/sh
# 対照 0（原本の表の再現）→ 対照 C1/C2/C4 → 本走行（forward）→ 逆順（C3）
set -e
cd "$(dirname "$0")"
python3 control_original_v1.py > results/control_original_v1.txt 2>&1
rm -rf "原本/exchange_scattering_matrix_fermionic_localization_transfer_preliminary_result_v1"
python3 run_exchange_v10q4.py --controls > results/controls.txt 2>&1
python3 run_exchange_v10q4.py --tag main > results/run_main.txt 2>&1
python3 run_exchange_v10q4.py --tag reverse --order reverse > results/run_reverse.txt 2>&1
echo done
