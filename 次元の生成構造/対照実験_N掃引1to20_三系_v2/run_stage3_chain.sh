#!/bin/zsh
# 第3段階 root probe: (mode, delta, deltatag) の並びを逐次実行
while [ $# -gt 0 ]; do
  M="$1"; D="$2"; DT="$3"; shift 3
  SUF="rootprobe-${M}-d${DT//./}"
  echo "=== probe $M delta=$D suffix=$SUF start $(date '+%H:%M:%S') ==="
  python3 -u run_finite_order_root_probe_stage3_v1.py \
    --mode "$M" --delta "$D" --output-suffix "$SUF" \
    --reference-npz "nsweep_${M}_T42000_d${DT}_N12_v2.npz" || { echo "FAILED $M $D"; exit 1; }
  echo "=== probe $M delta=$D done $(date '+%H:%M:%S') ==="
done
echo "=== chain complete $(date '+%H:%M:%S') ==="
