#!/bin/bash
# 90°骨格 D^{-1}KD=iW 同値とPerron検証（読出しのみ）
set -e
cd "$(dirname "$0")"
python3 verify_DKD_iW_equivalence_20260912.py
