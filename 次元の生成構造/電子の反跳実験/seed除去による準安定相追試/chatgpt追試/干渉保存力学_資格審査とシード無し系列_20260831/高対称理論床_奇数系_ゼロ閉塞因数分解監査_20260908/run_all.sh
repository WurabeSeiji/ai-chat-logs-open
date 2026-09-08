#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python analyze_highsym_odd_closure_factorization_v1.py
python verify_results_v1.py
sha256sum analyze_highsym_odd_closure_factorization_v1.py verify_results_v1.py odd_N_summary.csv exact_2wave_pair_candidates.csv constructive_primitive_zero_closure_groups.csv small_closure_examples.csv distance_class_zero_closure.csv RUN_METADATA.json README.md 高対称理論床_奇数系_ゼロ閉塞因数分解監査_20260908.md > SHA256SUMS.txt
