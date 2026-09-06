#!/bin/zsh
# A_standard の完全対照が通るまで待ってから残りの腕を実行
echo "=== A_standard の完了を待機中 $(date '+%H:%M:%S') ==="
while [ ! -f electron_affine16_A_standard_T42000_v1.json ]; do sleep 20; done
echo "=== A_standard 完了を検知 $(date '+%H:%M:%S') ==="
for A in "$@"; do
  echo "=== arm $A start $(date '+%H:%M:%S') ==="
  python3 -u run_electron_affine16_instability_diagnostic_v1.py --arm "$A" || { echo "FAILED $A"; exit 1; }
  echo "=== arm $A done $(date '+%H:%M:%S') ==="
done
echo "=== chain complete $(date '+%H:%M:%S') ==="
