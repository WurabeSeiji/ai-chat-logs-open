#!/usr/bin/env python3
"""39/124 マイクロズーム v1（インライン実行の正式ランナー化）

README「マイクロズーム（27点×1500衝突、amp 2.470-2.522）」を生成した実験の
再現コード。1/3 タング下端直下で、目標 rho = 39/124 = 0.3145161290...への
最接近距離とプラトー不在（幅 < 0.002 の上界）を確認する。
"""

from __future__ import annotations

import csv
import importlib.util
import math
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
TS_PATH = HERE / "run_tongue_spectroscopy_pre_v1.py"

spec = importlib.util.spec_from_file_location("ts_for_micro_zoom_v1", TS_PATH)
ts = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = ts
ts.__name__ = "ts_for_micro_zoom_v1"
spec.loader.exec_module(ts)

toy, base = ts.toy, ts.base

COLLISIONS = 1500
TAIL = 500
AMPLITUDES = np.linspace(2.470, 2.522, 27)
TARGET = 39.0 / 124.0


def main() -> None:
    params = base.Params(high_n=63, recursive_collision_count=COLLISIONS)
    sp = base.build_source_params(params)
    a_t, b_t = ts.make_templates(sp)
    v1 = ts.READOUTS["R2_X_power_V1"]

    rows = []
    for amp in AMPLITUDES:
        a = a_t.copy()
        b = float(amp) * b_t
        hist = np.zeros(COLLISIONS)
        for j in range(COLLISIONS):
            th = v1(a, b, sp)
            hist[j] = th / math.pi
            a, b = toy.rotate_ab(a, b, th)
        tail = hist[-TAIL:]
        mean = float(np.mean(tail))
        frac = Fraction(mean).limit_denominator(200)
        rows.append(
            {
                "amplitude": float(amp),
                "rho_tail_mean": mean,
                "tail_std": float(np.std(tail)),
                "dist_to_39_124": abs(mean - TARGET),
                "nearest_rational": f"{frac.numerator}/{frac.denominator}",
            }
        )
        print(
            f"amp={amp:.4f} rho={mean:.10f} d(39/124)={abs(mean-TARGET):.2e}"
        )

    best = min(rows, key=lambda r: r["dist_to_39_124"])
    print(
        f"closest approach: amp={best['amplitude']:.4f}"
        f" d={best['dist_to_39_124']:.2e} (plateau at 39/124: not detected)"
    )
    csv_path = HERE / "micro_zoom_39_124_rows_v1.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
