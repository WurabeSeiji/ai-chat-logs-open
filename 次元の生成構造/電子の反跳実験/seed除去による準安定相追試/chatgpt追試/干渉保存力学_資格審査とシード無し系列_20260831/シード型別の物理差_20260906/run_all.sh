#!/bin/zsh
# シード型別の物理差: 論文A正本 δ=1e-2 T42000 N=12 の5系＋真空を横並び比較（読み出しのみ）
cd "$(dirname "$0")"
python3 check_mode_physics_diff_v1.py
