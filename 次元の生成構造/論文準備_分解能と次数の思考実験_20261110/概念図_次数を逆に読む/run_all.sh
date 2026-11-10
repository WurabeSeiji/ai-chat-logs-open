#!/bin/sh
# 概念図 v0.1 の再生成
cd "$(dirname "$0")" || exit 1
mkdir -p figs
python3 fig_reverse_reading.py 2>&1 | grep -v findfont > figs/check_v0.1.txt
cat figs/check_v0.1.txt
