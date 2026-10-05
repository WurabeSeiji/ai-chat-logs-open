#!/bin/bash
# 一式の再現：本体（audit.txt、cascade_result.png）、検算 3 本（*_check.txt）、差分、SHA256SUMS。
# 本体はフォントを OS に合わせて落とすので、どの OS でも直接実行できる。
set -euo pipefail
cd "$(dirname "$0")"
python3 exchange_cascade.py
python3 coulomb_levels_check.py > /dev/null
python3 gravity_levels_check.py > /dev/null
python3 rates_check.py > /dev/null
diff -u ../exchange_latest_20261005/exchange_cascade.py exchange_cascade.py > CHANGES_vs_original.diff || true
shasum -a 256 exchange_cascade.py coulomb_levels_check.py gravity_levels_check.py rates_check.py run_all.sh \
    audit.txt coulomb_levels_check.txt gravity_levels_check.txt rates_check.txt cascade_result.png \
    CHANGES_vs_original.diff README.md 重力項目_調査と修正_ja_20261005.md > SHA256SUMS
echo "done: $(wc -l < SHA256SUMS) files in SHA256SUMS"
