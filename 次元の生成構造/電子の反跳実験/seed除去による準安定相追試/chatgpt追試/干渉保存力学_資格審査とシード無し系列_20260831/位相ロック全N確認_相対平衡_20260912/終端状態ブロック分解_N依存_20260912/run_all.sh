#!/bin/bash
# 追試B: 終端状態ブロック分解のN依存（読出しのみ）
set -e
cd "$(dirname "$0")"
python3 analyze_terminal_block_decomposition_20260912.py
