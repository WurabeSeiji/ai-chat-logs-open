#!/bin/sh
# 確認計算を実行して results を更新し、SHA256SUMS を再生成する
set -e
cd "$(dirname "$0")"
python3 -X utf8 check_values_v2c2.py > results_check_values_v2c2.txt
python3 -X utf8 check_values_two_pulses.py > results_check_values_two_pulses.txt
python3 -X utf8 check_values_v2c3.py > results_check_values_v2c3.txt
shasum -a 256 index_v2c2.html index_v2c3.html check_values_v2c2.py results_check_values_v2c2.txt check_values_v2c3.py results_check_values_v2c3.txt check_values_two_pulses.py results_check_values_two_pulses.txt \
  "論文原案_Rの軌道_同じ弦の上の二つの複素パルス_20261029.md" README.md run_all.sh figs/*.png > SHA256SUMS
echo "done: $(wc -l < SHA256SUMS) files"
