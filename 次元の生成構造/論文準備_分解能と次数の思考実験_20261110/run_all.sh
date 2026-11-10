#!/bin/bash
# 図 1〜3 の再生成
set -e
cd "$(dirname "$0")"
python3 図表/make_crops.py
bash 概念図_次数を逆に読む/run_all.sh
