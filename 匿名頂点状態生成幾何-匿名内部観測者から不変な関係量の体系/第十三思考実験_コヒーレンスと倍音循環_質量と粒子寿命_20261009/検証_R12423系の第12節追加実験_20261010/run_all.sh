#!/bin/sh
# 第十三思考実験 §12.2〜12.6 の追加実験（R_{124,23} 系）：記録 → 図化。
# 原本は同じ親フォルダの 検証_R070交換の正規化と振幅_20261009/original_copy/ から取り、SHA256 を照合する。
set -e
cd "$(dirname "$0")"
ORIG=original_copy/run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1.py
if [ ! -f "$ORIG" ]; then
  mkdir -p original_copy
  cp "../検証_R070交換の正規化と振幅_20261009/$ORIG" "$ORIG"
fi
echo "f815320f5632ae1b23ccade3a53b01e9110ad770da8407e2b10fa6065ef1695c  $ORIG" | sha256sum -c -
python3 run_section12_R12423.py > results/run_stdout.txt
python3 plot_section12_R12423.py > figures/plot_stdout.txt
echo "done (about 2 min). check: sha256sum -c SHA256SUMS"
