#!/bin/sh
set -eu
cd "$(dirname "$0")"
python3 audit_highsym_even_photon_pairs_v2.py
python3 verify_results.py
