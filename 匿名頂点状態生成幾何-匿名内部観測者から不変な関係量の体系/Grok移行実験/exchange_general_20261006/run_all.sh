#!/bin/bash
# 一式の再現：対照実験（control_report.md）、引力の走行（audit_general_attract.txt）、斥力の走行（audit_general_repel.txt）、SHA256SUMS。
# 共通部分は ../exchange_rel_20261005/exchange_cascade.py を import する（そのフォルダが必要）。
set -euo pipefail
cd "$(dirname "$0")"
python3 control_general_vs_closed_form.py > /dev/null
python3 two_body_general.py -1 > /dev/null
python3 two_body_general.py +1 > /dev/null
shasum -a 256 two_body_general.py control_general_vs_closed_form.py run_all.sh \
    control_report.md audit_general_attract.txt audit_general_repel.txt README.md > SHA256SUMS
echo "done: $(wc -l < SHA256SUMS) files in SHA256SUMS"
