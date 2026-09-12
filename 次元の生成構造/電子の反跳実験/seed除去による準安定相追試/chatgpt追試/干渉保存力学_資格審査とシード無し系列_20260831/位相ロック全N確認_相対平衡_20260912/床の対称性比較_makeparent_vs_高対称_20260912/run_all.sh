#!/bin/bash
# 床の対称性比較（読出しのみ）
set -e
cd "$(dirname "$0")"
python3 compare_floor_symmetry_20260912.py
