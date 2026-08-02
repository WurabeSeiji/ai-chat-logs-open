#!/usr/bin/env python3
"""1/3 タング精密スキャン v1（インライン実行の正式ランナー化）

README の表「精密スキャン（800衝突・tail 300）」を生成した実験の再現コード。
振幅 2.4..3.6 の13点で、V1（coherent X power）読出しの下の後期回転数を測り、
タング内で tail 平均が厳密に 1/3 に一致することを確認する。
"""

from __future__ import annotations

import csv
import importlib.util
import math
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
PROBE_PATH = HERE / "run_mode_locking_probe_pre_v1.py"

spec = importlib.util.spec_from_file_location("probe_for_fine_scan_v1", PROBE_PATH)
probe = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = probe
probe.__name__ = "probe_for_fine_scan_v1"
spec.loader.exec_module(probe)

toy, base = probe.toy, probe.base

COLLISIONS = 800
TAIL = 300
AMPLITUDES = np.linspace(2.4, 3.6, 13)
EXACT_TOL = 1.0e-11


def main() -> None:
    params = base.Params(high_n=63, recursive_collision_count=COLLISIONS)
    sp = base.build_source_params(params)
    a_t, b_t = probe.make_templates(sp)

    rows = []
    n_exact = 0
    for amp in AMPLITUDES:
        a = a_t.copy()
        b = float(amp) * b_t
        hist = np.zeros(COLLISIONS)
        for j in range(COLLISIONS):
            th = probe.theta_coherent(a, b, sp)
            hist[j] = th / math.pi
            a, b = toy.rotate_ab(a, b, th)
        tail = hist[-TAIL:]
        mean = float(np.mean(tail))
        std = float(np.std(tail))
        dist = abs(mean - 1.0 / 3.0)
        exact = dist <= EXACT_TOL
        n_exact += int(exact)
        rows.append(
            {
                "amplitude": float(amp),
                "theta_over_pi_tail_mean": mean,
                "tail_std": std,
                "dist_to_one_third": dist,
                "locked": exact,
            }
        )
        print(
            f"amp={amp:.3f} rho={mean:.12f} std={std:.2e}"
            f" |rho-1/3|={dist:.2e} {'LOCK' if exact else ''}"
        )

    csv_path = HERE / "one_third_tongue_fine_scan_rows_v1.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"locked points: {n_exact}/{len(rows)} (expected: interior 2.6..3.3 exact)")


if __name__ == "__main__":
    main()
