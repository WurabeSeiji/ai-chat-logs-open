#!/bin/bash
# 図の再現（日本語：v1 のプログラムで再生成し、フォルダ直下の SVG とバイト比較）と英語の図の生成
set -e
cd "$(dirname "$0")/.."
python3 make_figures_quaternion_crystallization_20261021.py 図表/ja > 図表/ja_run_log.txt 2>&1
for f in fig0*.svg; do cmp -s "$f" "図表/ja/$f" && echo "一致 $f" || echo "不一致 $f"; done > 図表/ja_compare.txt
python3 make_figures_quaternion_crystallization_20261021_v1.1.py 図表/ja_v1.1check ja > /dev/null 2>&1
for f in fig0*.svg; do cmp -s "$f" "図表/ja_v1.1check/$f" && echo "一致 $f" || echo "不一致 $f"; done > 図表/ja_v1.1_compare.txt
rm -rf 図表/ja_v1.1check
python3 make_figures_quaternion_crystallization_20261021_v1.1.py 図表/en en > 図表/en_run_log.txt 2>&1
cd 図表 && shasum -a 256 ja/* en/* *.txt > SHA256SUMS.txt
