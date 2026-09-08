#!/bin/sh
set -eu
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
python3 "$HERE/audit_simplex_dimension_evolution_v1.py" "$HERE/inputs" "$HERE"
