#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
python3 make_fig3_5color_stage123_v1.py
python3 check_direction_axes_mapping_v1.py
echo "run_all done"
