#!/bin/bash
# 追試A: ランプ超過＝時変非正規ヤコビアン積 の直接Lyapunov照合（読出しのみ）
set -e
cd "$(dirname "$0")"
python3 analyze_ramp_lyapunov_20260912.py
