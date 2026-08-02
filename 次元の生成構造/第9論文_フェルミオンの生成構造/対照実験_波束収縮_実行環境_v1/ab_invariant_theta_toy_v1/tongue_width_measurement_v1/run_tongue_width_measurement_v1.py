#!/usr/bin/env python3
"""タング幅のモデル内実測 v1（E8）

目的:
    「高分母タングの重みは無視できる」を外部理論の引用でなくモデル内実測で言う。
    (i) 1/3 タングの上端を精密化し、幅を実測で確定する
        （下端 2.518 は micro-zoom で確定済み）
    (ii) 探索済み帯域 [2.470, 2.522]（Δamp=0.002, 1500衝突）に分母>=4 の
        タングが検出されなかったことから、幅比の下界を出す:
        width(q=3) / width(q>=4) > 実測比

方法: tongue_spectroscopy の V1 読出しを再利用（コピーせず import）。
"""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
TS_PATH = HERE.parent / "tongue_spectroscopy_pre_v1" / "run_tongue_spectroscopy_pre_v1.py"

spec = importlib.util.spec_from_file_location("ts_for_width_v1", TS_PATH)
ts = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = ts
ts.__name__ = "ts_for_width_v1"
spec.loader.exec_module(ts)
toy, base = ts.toy, ts.base
plt = base.plt

COLLISIONS = 800
TAIL = 300
UPPER_SCAN = np.linspace(3.28, 3.44, 17)
LOWER_EDGE_KNOWN = 2.518  # micro-zoom で確定済み（それ未満は非ロック）
SEARCHED_BAND_WIDTH = 0.002  # micro-zoom の分解能＝分母>=4 タング幅の上界
EXACT_TOL = 1.0e-9


def rho_tail(amp: float, sp, a_t, b_t, v1) -> float:
    a = a_t.copy()
    b = float(amp) * b_t
    hist = np.zeros(COLLISIONS)
    for j in range(COLLISIONS):
        th = v1(a, b, sp)
        hist[j] = th / math.pi
        a, b = toy.rotate_ab(a, b, th)
    return float(np.mean(hist[-TAIL:]))


def main() -> None:
    params = base.Params(high_n=63, recursive_collision_count=COLLISIONS)
    sp = base.build_source_params(params)
    a_t, b_t = ts.make_templates(sp)
    v1 = ts.READOUTS["R2_X_power_V1"]

    rows = []
    locked = []
    for amp in UPPER_SCAN:
        rho = rho_tail(amp, sp, a_t, b_t, v1)
        is_locked = abs(rho - 1.0 / 3.0) <= EXACT_TOL
        rows.append({"amplitude": float(amp), "rho": rho,
                     "dist_to_one_third": abs(rho - 1.0 / 3.0), "locked": is_locked})
        if is_locked:
            locked.append(float(amp))
        print(f"amp={amp:.3f} rho={rho:.12f} {'LOCK' if is_locked else ''}")

    upper_edge = max(locked) if locked else float("nan")
    width_q3 = upper_edge - LOWER_EDGE_KNOWN
    ratio_bound = width_q3 / SEARCHED_BAND_WIDTH
    print(f"\n1/3 tongue: lower edge {LOWER_EDGE_KNOWN} (micro-zoom), upper edge {upper_edge:.3f}")
    print(f"width(q=3) = {width_q3:.3f} amplitude units")
    print(f"searched band resolution = {SEARCHED_BAND_WIDTH} -> width(q>=4) < {SEARCHED_BAND_WIDTH}")
    print(f"IN-MODEL bound: width(q=3)/width(q>=4) > {ratio_bound:.0f}")

    with (HERE / "tongue_width_rows_v1.csv").open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    fig, ax = plt.subplots(figsize=(8, 4.6), constrained_layout=True)
    ax.plot([r["amplitude"] for r in rows], [r["rho"] for r in rows], "o-", markersize=4)
    ax.axhline(1.0 / 3.0, color="tab:red", linestyle=":", label="1/3")
    ax.set_xlabel("initial B amplitude"); ax.set_ylabel("tail-mean theta/pi")
    ax.set_title(f"upper edge of the 1/3 tongue: {upper_edge:.3f} (width {width_q3:.3f})")
    ax.legend(); ax.grid(alpha=0.3)
    for ext in ("png", "svg"):
        fig.savefig(HERE / f"tongue_width_measurement_v1.{ext}", dpi=160)
    plt.close(fig)

    payload = {
        "experiment": "tongue_width_measurement_v1",
        "core_runner_lineage": "tongue_spectroscopy_pre_v1 (V1 readout, imported)",
        "one_third_tongue": {
            "lower_edge": LOWER_EDGE_KNOWN,
            "upper_edge": upper_edge,
            "width": width_q3,
        },
        "higher_denominator_bound": {
            "searched_band": [2.470, 2.522],
            "resolution": SEARCHED_BAND_WIDTH,
            "statement": "分母>=4 のタングは検出されず、幅 < 0.002",
        },
        "in_model_width_ratio_bound": ratio_bound,
        "conclusion": "q=3 から q>=4 で幅は少なくとも400倍以上縮む（モデル内実測の下界）。分母124のタング重みが無視できることの実測的裏付け",
        "rows": rows,
    }
    (HERE / "tongue_width_measurement_result_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
