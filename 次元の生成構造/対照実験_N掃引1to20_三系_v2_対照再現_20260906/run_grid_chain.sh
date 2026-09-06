#!/bin/zsh
# 第1段階格子拡張: 1つのδについて control → phase-balanced を逐次実行
D="$1"; TAG="$2"; REF="$3"
export PB_DELTA="$D"
[ -n "$REF" ] && export PB_REF_NPZ="$REF"
echo "=== chain $TAG : delta=$D ref=${REF:-default} start $(date '+%H:%M:%S') ==="
python3 -u run_phase_balanced_mixed_grid_v1.py run --arm control --replicate ctl1 --allow-science-run || { echo "CONTROL FAILED"; exit 1; }
echo "=== control done $(date '+%H:%M:%S') ==="
python3 -u run_phase_balanced_mixed_grid_v1.py run --arm phase-balanced --replicate r1 --allow-science-run || { echo "BALANCED FAILED"; exit 1; }
echo "=== chain $TAG done $(date '+%H:%M:%S') ==="
