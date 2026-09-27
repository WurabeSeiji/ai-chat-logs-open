#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations
import csv
from pathlib import Path
import numpy as np
from paper6_init_G_q0_c1 import CASES, init_case, P, N, C, D, K1, K2, K3, K4, DV, Q, QB

HERE = Path(__file__).resolve().parent
OUT = HERE / "generated_initial_rows_G_q0_c1.csv"

header = ["case_id","micro","q_index","P","N","C","D","k1p","k2p","k3p","k4p","dP","Q"]
with OUT.open("w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(header)
    for case_id, *_ in CASES:
        z, meta = init_case(case_id)
        w.writerow([
            case_id, 0, int(np.argmax(z[QB:QB+11])), z[P], z[N], z[C], z[D],
            z[K1], z[K2], z[K3], z[K4], z[DV], z[Q]
        ])
print(OUT)
