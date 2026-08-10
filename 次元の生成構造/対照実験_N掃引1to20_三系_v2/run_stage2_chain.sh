#!/bin/zsh
for A in "$@"; do
  echo "=== arm $A start $(date '+%H:%M:%S') ==="
  python3 -u run_electron_affine16_instability_diagnostic_v1.py --arm "$A" || { echo "FAILED $A"; exit 1; }
  echo "=== arm $A done $(date '+%H:%M:%S') ==="
done
echo "=== chain complete $(date '+%H:%M:%S') ==="
