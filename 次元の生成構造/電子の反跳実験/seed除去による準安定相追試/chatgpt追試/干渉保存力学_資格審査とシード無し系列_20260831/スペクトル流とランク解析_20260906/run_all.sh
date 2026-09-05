#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
python3 analyze_spectral_flow_rank_v1.py
echo "run_all done"
