#!/bin/bash
# 凝縮体構造の状態空間内調査（読み出しのみ）
set -euo pipefail
cd "$(dirname "$0")"
python3 check_vertex_star_closure_v1.py
python3 check_transition_cowinding_v1.py
python3 check_anonymity_entropy_v1.py
python3 check_temperature_readout_v1.py
echo "run_all done"
