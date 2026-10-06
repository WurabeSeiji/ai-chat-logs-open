#!/bin/bash
# 一式の再現：水素の回帰テスト（regression_hydrogen.md）、H₂ の検証（validation_h2.md）、H₂ の走行（audit_h2_v1J3.txt、audit_h2_v3J5.txt）、SHA256SUMS。
# 水素の部品は ../exchange_rel_20261005/exchange_cascade.py を import する（そのフォルダが必要）。H₂ の部品は h2_data/ の表だけで閉じる。
set -euo pipefail
cd "$(dirname "$0")"
python3 regression_hydrogen.py > /dev/null
python3 validation_h2.py > /dev/null
python3 run_h2.py 1 3 > /dev/null
python3 run_h2.py 3 5 > /dev/null
python3 run_h2.py 14 0 > /dev/null
python3 run_h2.py 14 1 > /dev/null
python3 radiative_association_h2.py > /dev/null
shasum -a 256 engine.py system_hydrogen.py regression_hydrogen.py system_h2.py validation_h2.py run_h2.py radiative_association_h2.py run_all.sh \
    regression_hydrogen.md validation_h2.md audit_h2_v1J3.txt audit_h2_v3J5.txt audit_h2_v14J0.txt audit_h2_v14J1.txt radiative_association_h2.md README.md \
    validation_h2.json audit_h2_v1J3.json audit_h2_v3J5.json audit_h2_v14J0.json audit_h2_v14J1.json radiative_association_h2.json \
    h2_data/h2_V_BO_pachucki2010.dat h2_data/h2_Q_WSD1998.dat h2_data/h2_g_PK2011.dat h2_data/SOURCES.md \
    h2_data/validation/roueff2019_table2.dat h2_data/validation/roueff2019_ReadMe.txt h2_data/validation/h2spectre74_H2_E2full.dat > SHA256SUMS
echo "done: $(wc -l < SHA256SUMS) files in SHA256SUMS"
