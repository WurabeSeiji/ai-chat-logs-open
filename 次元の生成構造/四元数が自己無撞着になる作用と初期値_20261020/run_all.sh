#!/bin/bash
set -e
cd "$(dirname "$0")"
python3 quaternion_selfconsistent_actions_20261020.py > results/run_log.txt 2> results/run_stderr.txt
rmdir figures 2>/dev/null || true   # 読み込み元の閉塞プログラムが作る空のフォルダ
shasum -a 256 *.py results/* > SHA256SUMS.txt
