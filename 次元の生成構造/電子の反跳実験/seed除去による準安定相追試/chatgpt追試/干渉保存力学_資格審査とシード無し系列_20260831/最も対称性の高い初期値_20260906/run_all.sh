#!/bin/zsh
# 最も対称性の高い初期値 N=3..40 の生成・検査・一覧・図・文書（新規走行なし）
cd "$(dirname "$0")"
python3 make_symmetric_parents_N3_N40_v1.py
