#!/bin/sh
# 確認計算を実行して results を更新し、SHA256SUMS を再生成する
set -e
cd "$(dirname "$0")"
python3 -X utf8 check_values_v2c.py > results_check_values_v2c.txt
python3 -X utf8 check_square_sum_parity.py > results_check_square_sum_parity.txt
shasum -a 256 index_v2c.html check_values_v2c.py results_check_values_v2c.txt check_square_sum_parity.py results_check_square_sum_parity.txt \
  "論文原案_Rの検討_複素化した閉じた弦のパルス_20261029.md" README.md run_all.sh figs/*.png > SHA256SUMS
echo "done: $(wc -l < SHA256SUMS) files"
