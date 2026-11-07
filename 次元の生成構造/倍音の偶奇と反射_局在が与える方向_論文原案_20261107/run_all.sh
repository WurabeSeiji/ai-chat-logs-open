#!/bin/sh
# 図（SVG）と数値（results/values.json）の再生成
set -e
cd "$(dirname "$0")"
python3 make_figures.py > results/make_figures_log.txt 2>&1
echo done
