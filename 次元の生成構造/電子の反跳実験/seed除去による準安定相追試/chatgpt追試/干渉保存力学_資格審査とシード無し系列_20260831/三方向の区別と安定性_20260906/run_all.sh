#!/bin/zsh
# 三方向の区別と安定性: 3検査を順に実行（読み出しのみ・新規走行なし）
cd "$(dirname "$0")"
python3 check_three_direction_frame_stability_v1.py
python3 check_inplane_direction_conservation_v1.py
python3 check_final_inplane_ellipse_v1.py
python3 check_frame_coordinate_trajectory_v1.py
