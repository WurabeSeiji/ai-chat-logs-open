#!/bin/sh
set -e
cd "$(dirname "$0")"
python3 check1_true_square_and_norm.py
python3 check2_triality_d4.py
python3 check2b_triality_as_quaternion_product.py
