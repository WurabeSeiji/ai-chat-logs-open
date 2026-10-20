#!/bin/bash
set -e
cd "$(dirname "$0")"
python3 octonion_selfconsistent_actions_20261020.py > results/run_log.txt 2> results/run_stderr.txt
python3 check_second_layer_identity_20261020.py > results/second_layer_identity_log.txt 2>/dev/null
rmdir figures 2>/dev/null || true   # 読み込み元の閉塞プログラムが作る空のフォルダ
shasum -a 256 octonion_selfconsistent_actions_20261020.py check_second_layer_identity_20261020.py results/* > SHA256SUMS.txt
