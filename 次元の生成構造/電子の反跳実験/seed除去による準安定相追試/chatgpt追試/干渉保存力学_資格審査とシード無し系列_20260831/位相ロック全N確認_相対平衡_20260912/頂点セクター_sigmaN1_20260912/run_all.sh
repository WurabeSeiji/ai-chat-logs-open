#!/bin/bash
# 追試C: 頂点セクターとσ=N-1の対応（読出しのみ）
set -e
cd "$(dirname "$0")"
python3 analyze_vertex_sector_sigmaN1_20260912.py
