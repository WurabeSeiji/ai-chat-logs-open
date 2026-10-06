#!/bin/bash
# 一式の再現：対照実験（control_report.md）、引力の走行（audit_general_attract.txt）、斥力の走行（audit_general_repel.txt）、SHA256SUMS。
# 共通部分は ../exchange_rel_20261005/exchange_cascade.py を import する（そのフォルダが必要）。
set -euo pipefail
cd "$(dirname "$0")"
python3 control_general_vs_closed_form.py > /dev/null
python3 two_body_general.py ep > /dev/null
python3 two_body_general.py ee > /dev/null
python3 two_body_general.py pp > /dev/null
python3 scattering_same_sign.py > /dev/null
python3 absorber_temperature.py > /dev/null
rm -f audit_general_attract.txt audit_general_repel.txt
shasum -a 256 two_body_general.py control_general_vs_closed_form.py scattering_same_sign.py absorber_temperature.py run_all.sh \
    control_report.md audit_general_ep.txt audit_general_ee.txt audit_general_pp.txt scattering_same_sign.md \
    absorber_temperature.md audit_temperature_T0.0000.txt audit_temperature_T2.7255.txt README.md > SHA256SUMS
echo "done: $(wc -l < SHA256SUMS) files in SHA256SUMS"
