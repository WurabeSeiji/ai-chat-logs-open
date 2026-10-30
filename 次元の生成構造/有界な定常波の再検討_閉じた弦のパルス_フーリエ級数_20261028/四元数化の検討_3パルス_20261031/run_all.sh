#!/bin/sh
# 確認計算と原本との対照を実行して results を更新し、SHA256SUMS を再生成する。
#   sh run_all.sh           … Python の確認計算・差分・HTML 対照の集計（既にある DOM を使う）
#   sh run_all.sh --shoot   … 加えて描画記録と DOM の取り直し（Chrome と通信が要る）
set -e
cd "$(dirname "$0")"
python3 -X utf8 check_values_v2q3.py > results/results_check_values_v2q3.txt
diff -u 原本/index_v2c3.html index_v2q3.html > results/diff_index_v2c3_v2q3.txt || true
diff -u 原本/check_values_v2c3.py check_values_v2q3.py > results/diff_check_values_v2c3_v2q3.txt || true
if [ "$1" = "--shoot" ]; then sh shoot.sh; fi
python3 -X utf8 control_html.py > results/results_control_html.txt
shasum -a 256 index_v2q3.html check_values_v2q3.py control_html.py shoot.sh run_all.sh README.md \
  "分析_四元数化した閉じた弦の3パルス_複素数との対応表_20261031.md" 原本/* results/*.txt results/*.html figs/*.png > SHA256SUMS
echo "done: $(wc -l < SHA256SUMS) files"
