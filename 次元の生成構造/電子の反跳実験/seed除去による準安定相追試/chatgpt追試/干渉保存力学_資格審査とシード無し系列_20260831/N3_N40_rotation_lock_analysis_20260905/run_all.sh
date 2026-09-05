#!/bin/bash
# 回転数ロック全系列調査（読み出しのみ・入力SHAゲート＋対照テスト内蔵）
set -euo pipefail
cd "$(dirname "$0")"
python3 analyze_rotation_lock_N3_N40_v1.py
python3 check_N5_convergence_v1.py
python3 check_spectrum_harmonics_v1.py
echo "run_all done"
