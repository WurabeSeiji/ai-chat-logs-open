#!/usr/bin/env bash
set -euo pipefail
python search_dual_balance_strict_v2.py --q 0.3  --out dual_q0p3_v2.json  > dual_q0p3_v2_stdout.txt
python search_dual_balance_strict_v2.py --q 0.6  --out dual_q0p6_v2.json  > dual_q0p6_v2_stdout.txt
python search_dual_balance_strict_v2.py --q 0.9  --out dual_q0p9_v2.json  > dual_q0p9_v2_stdout.txt
python search_dual_balance_strict_v2.py --q 0.99 --out dual_q0p99_v2.json > dual_q0p99_v2_stdout.txt
python verify_dual_balance_repro_v1.py > dual_balance_repro_verification_v1_stdout.txt
