#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
DATA_ROOT="$ROOT/data"
OUT="$HERE/../results/recomputed"
python3 "$HERE/reproduce_hm_rank3_audit.py" \
  --data-root "$DATA_ROOT" \
  --out "$OUT" \
  --n-min 3 --n-max 16 \
  --windows 4096 8192 16384 \
  --baseline-cmax 100
