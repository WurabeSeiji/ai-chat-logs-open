#!/usr/bin/env python3
"""混合結合・予備テスト v1 —— クーロン前係数 q_Aq_B の合成則判定

目的（2026-07-30 メモリ「残る決定打」項目2、提案のまま未実行だった実験）:
    二閉包が固有の透過分率 f_A, f_B を持つとき、実効交換強度
    T_eff(f_A, f_B) がどう合成されるかを測る。
      - 積型   T_eff ∝ f_A f_B      → クーロン前係数 q_Aq_B/4π がそのまま立つ
      - 平方根型 ∝ sqrt(f_A f_B)     → 頂点量は sqrt(T)
      - 平均型  ∝ (f_A+f_B)/2       → QED 型頂点合成はゼロ次に存在しない

方法（合成則を仮定しない状態構成）:
    各側の固有分率 f を、全ボゾン的単一倍音（マスク外, f=0）と
    全フェルミオン的単一倍音（マスク内, f=1）の二倍音混合
        state(f) = sqrt(1-f)·s_bos + sqrt(f)·s_fer
    で連続に調整する（パワー直交性より分率は厳密に f）。
    A側: s_bos = A(1)（ビン0,2）, s_fer = A(5)（ビン4,6）
    B側: s_bos = B(1)（ビン0,-2）, s_fer = B(7)（ビン6,-8）
    ビン表で各単一状態の分率が厳密に 0 / 1 であることを先に検証する。

測定:
    Part 1: 実効交換強度 R_joint = sin^2(theta_from_ab) を 6x6 グリッドで測定し、
            積・平均・平方根の三仮説と比較（ゼロ次合成則の確定）
    Part 2: 実際のチャネルノルム移乗（16衝突のスイング）と局在移乗 L の
            スイングを同グリッドで測定し、同じ三仮説へ回帰
"""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
TOY_RUNNER_PATH = HERE.parent / "run_ab_invariant_theta_toy_v1.py"

F_GRID = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
SWING_COLLISIONS = 16
FRACTION_TOLERANCE = 1.0e-12
LAW_TOLERANCE = 1.0e-12


def load_toy_module() -> Any:
    spec = importlib.util.spec_from_file_location(
        "ab_theta_toy_for_mixed_coupling_v1", TOY_RUNNER_PATH
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load toy runner: {TOY_RUNNER_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


toy = load_toy_module()
base = toy.base
plt = base.plt


def unit_norm(v: np.ndarray) -> np.ndarray:
    n = math.sqrt(float(np.vdot(v, v).real))
    return v / n


def make_single(k: int, which: str, source_params: Any) -> np.ndarray:
    case = base.explicit_packet_case(
        mode=f"mix_{which}_{k}", packet_a=(k,), packet_b=(k,)
    )
    return unit_norm(base.make_case_state(source_params, case, which, hair_enabled=True))


def masked_fraction(state: np.ndarray, source_params: Any) -> float:
    zero = np.zeros_like(state)
    frequencies, power = toy.combined_chi_power(state, zero, source_params)
    f = np.abs(frequencies)
    mask = (f >= 4) & (f % 2 == 0)
    total = float(np.sum(power))
    return float(np.sum(power[mask])) / total


def mixed_state(f: float, s_bos: np.ndarray, s_fer: np.ndarray) -> np.ndarray:
    return unit_norm(math.sqrt(1.0 - f) * s_bos + math.sqrt(f) * s_fer)


def channel_norms(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    return float(np.vdot(a, a).real), float(np.vdot(b, b).real)


def main() -> None:
    params = base.Params(high_n=63, recursive_collision_count=SWING_COLLISIONS)
    source_params = base.build_source_params(params)
    metric_context = base.MetricContext(source_params)

    a_bos = make_single(1, "A", source_params)
    a_fer = make_single(5, "A", source_params)
    b_bos = make_single(1, "B", source_params)
    b_fer = make_single(7, "B", source_params)

    prep = {
        "A_bosonic_k1": masked_fraction(a_bos, source_params),
        "A_fermionic_k5": masked_fraction(a_fer, source_params),
        "B_bosonic_k1": masked_fraction(b_bos, source_params),
        "B_fermionic_k7": masked_fraction(b_fer, source_params),
        "cross_overlap_abs_max": max(
            abs(complex(np.vdot(a_bos, a_fer))),
            abs(complex(np.vdot(b_bos, b_fer))),
        ),
    }
    prep_ok = (
        prep["A_bosonic_k1"] <= FRACTION_TOLERANCE
        and abs(prep["A_fermionic_k5"] - 1.0) <= FRACTION_TOLERANCE
        and prep["B_bosonic_k1"] <= FRACTION_TOLERANCE
        and abs(prep["B_fermionic_k7"] - 1.0) <= FRACTION_TOLERANCE
        and prep["cross_overlap_abs_max"] <= 1.0e-12
    )
    print("state preparation:", "PASS" if prep_ok else "CHECK", prep)

    rows: list[dict[str, Any]] = []
    for f_a in F_GRID:
        a0 = mixed_state(f_a, a_bos, a_fer)
        for f_b in F_GRID:
            b0 = mixed_state(f_b, b_bos, b_fer)
            readout = toy.theta_from_ab(a0, b0, source_params)
            r_joint = readout.reflection_rate

            a, b = a0.copy(), b0.copy()
            nb0 = channel_norms(a, b)[1]
            lb0 = toy.state_metrics(b, metric_context)["L"]
            norm_swing = 0.0
            loc_swing = 0.0
            for _ in range(SWING_COLLISIONS):
                r = toy.theta_from_ab(a, b, source_params)
                a, b = toy.rotate_ab(a, b, r.theta)
                nb = channel_norms(a, b)[1]
                lb = toy.state_metrics(b, metric_context)["L"]
                norm_swing = max(norm_swing, abs(nb - nb0))
                loc_swing = max(loc_swing, abs(lb - lb0))

            rows.append(
                {
                    "f_A": f_a,
                    "f_B": f_b,
                    "R_joint": r_joint,
                    "hyp_mean": (f_a + f_b) / 2.0,
                    "hyp_product": f_a * f_b,
                    "hyp_sqrt_product": math.sqrt(f_a * f_b),
                    "norm_swing": norm_swing,
                    "loc_swing": loc_swing,
                }
            )

    # ---- Part 1 判定: R_joint がどの合成則に一致するか ----
    err = {
        law: max(abs(row["R_joint"] - row[f"hyp_{law}"]) for row in rows)
        for law in ("mean", "product", "sqrt_product")
    }
    matched = [law for law, e in err.items() if e <= LAW_TOLERANCE]
    print("Part1 composition law for R_joint:")
    for law, e in err.items():
        print(f"   {law:<13} max|err|={e:.3e}" + ("  <== MATCH" if e <= LAW_TOLERANCE else ""))

    # ---- Part 2 回帰: スイングは f_A f_B 積構造を持つか ----
    def correlate(target_key: str) -> dict[str, float]:
        y = np.asarray([row[target_key] for row in rows])
        result = {}
        for law in ("mean", "product", "sqrt_product"):
            x = np.asarray([row[f"hyp_{law}"] for row in rows])
            if np.std(x) == 0 or np.std(y) == 0:
                result[law] = float("nan")
                continue
            result[law] = float(np.corrcoef(x, y)[0, 1])
        return result

    corr_norm = correlate("norm_swing")
    corr_loc = correlate("loc_swing")
    print("Part2 correlation of transfer swings vs laws:")
    print(f"   norm_swing: {corr_norm}")
    print(f"   loc_swing : {corr_loc}")

    csv_path = HERE / "mixed_coupling_pre_rows_v1.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
    n = len(F_GRID)
    for ax, key, title in (
        (axes[0], "R_joint", "R_joint (exchange strength)"),
        (axes[1], "norm_swing", "channel-norm swing (16 collisions)"),
        (axes[2], "loc_swing", "localization swing"),
    ):
        grid = np.asarray([row[key] for row in rows]).reshape(n, n)
        im = ax.imshow(grid, origin="lower", extent=(0, 1, 0, 1), aspect="auto")
        ax.set_xlabel("f_B")
        ax.set_ylabel("f_A")
        ax.set_title(title)
        fig.colorbar(im, ax=ax, shrink=0.85)
    fig.suptitle("Mixed coupling: composition law of exchange strength and transfer")
    figure_names = []
    for ext in ("png", "svg"):
        path = HERE / f"mixed_coupling_pre_v1.{ext}"
        fig.savefig(path, dpi=160)
        figure_names.append(path.name)
    plt.close(fig)

    payload = {
        "experiment": "mixed_coupling_pre_v1",
        "purpose": "Coulomb prefactor composition test (memory 2026-07-30 item 2, previously proposed but never run)",
        "core_runner": {"path": TOY_RUNNER_PATH.name, "sha256": toy.sha256(TOY_RUNNER_PATH)},
        "state_preparation": {**prep, "verdict": "PASS" if prep_ok else "CHECK"},
        "part1_law_errors": err,
        "part1_matched_laws": matched,
        "part2_correlations": {"norm_swing": corr_norm, "loc_swing": corr_loc},
        "rows": rows,
        "figures": figure_names,
    }
    (HERE / "mixed_coupling_pre_result_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
