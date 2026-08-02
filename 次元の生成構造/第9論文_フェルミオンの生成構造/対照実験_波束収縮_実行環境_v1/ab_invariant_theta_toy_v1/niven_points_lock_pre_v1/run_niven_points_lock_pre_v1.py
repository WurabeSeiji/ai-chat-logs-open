#!/usr/bin/env python3
"""Niven 交点ロック・予備テスト v1

背景:
    数え上げ有理数 R 族と有理角ロック族の交点は、Niven の定理により
    R ∈ {0, 1/4, 1/2, 3/4, 1} の5点に限られる。R=1/2（周期8, (m,n)=(1,4)）は
    counting_rule_lock_bridge_pre_v1 で実証済み。本テストは残る非自明点を検証する。

予言（測定前にコード内定数として固定）:
    P1: 混合パリティ支持（偶数 n 個 + 奇数>=5 を n 個）は正準マスク下で
        R = 1/4 に厳密着地し、theta = pi/6、実回転周期 12、(m,n) = (1,3)
        （奇数位数）で厳密ロックする。
    P2: even>=2 変種マスク下の odds 5..63 等重み状態は R = 3/4 に厳密着地し
        （規約監査で測定済み）、theta = pi/3、周期 6、(m,n) = (1,6) でロックする。
    P3: R = 0 状態（B が基本波のみ）は不変（theta = 0）。
    A:  アンカー = R = 1/2（{1,3}-free 支持、周期 8、(1,4)）の再現。

設計境界:
    - 散乱本体・回転は無変更。振幅探索は一切使用しない（全状態が等重み数え上げ構成）
    - P2 の変種マスクは規約監査 v1 と同一定義の再利用であり、その旨を明示
"""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Callable

import numpy as np


HERE = Path(__file__).resolve().parent
TOY_RUNNER_PATH = HERE.parent / "run_ab_invariant_theta_toy_v1.py"

HIGH_N = 63
LOCK_COLLISIONS = 60
RECURRENCE_TOLERANCE = 1.0e-8
R_TOLERANCE = 1.0e-12
RATIONAL_MAX_DEN = 200


def load_toy() -> Any:
    spec = importlib.util.spec_from_file_location("toy_for_niven_v1", TOY_RUNNER_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


toy = load_toy()
base = toy.base
plt = base.plt


def canonical_theta(a: np.ndarray, b: np.ndarray, sp: Any) -> float:
    return toy.theta_from_ab(a, b, sp).theta


def m2_theta(a: np.ndarray, b: np.ndarray, sp: Any) -> float:
    """規約監査 v1 の M2_even_ge2 変種（同一定義の再利用）。"""
    shape = (sp.chi_grid_n, sp.eta_grid_n)
    af = np.fft.fft(a.reshape(shape), axis=0, norm="ortho")
    bf = np.fft.fft(b.reshape(shape), axis=0, norm="ortho")
    freqs = np.rint(np.fft.fftfreq(sp.chi_grid_n, d=1.0 / sp.chi_grid_n)).astype(int)
    f = np.abs(freqs)
    mask = (f >= 2) & (f % 2 == 0)
    power = np.sum(np.abs(af) ** 2 + np.abs(bf) ** 2, axis=1)
    pf = float(np.sum(power[mask]))
    pb = float(np.sum(power[~mask]))
    return math.atan2(math.sqrt(max(pf, 0.0)), math.sqrt(max(pb, 0.0)))


def unit_norm(v: np.ndarray) -> np.ndarray:
    return v / math.sqrt(float(np.vdot(v, v).real))


def make_pair(packet: tuple[int, ...], sp: Any) -> tuple[np.ndarray, np.ndarray]:
    case = base.explicit_packet_case(
        mode=f"niven_len{len(packet)}_{packet[0]}", packet_a=(1,), packet_b=packet
    )
    a = unit_norm(base.make_case_state(sp, case, "A", hair_enabled=True))
    b = unit_norm(base.make_case_state(sp, case, "B", hair_enabled=True))
    return a, b


CASES = (
    # (case_id, support, theta_fn_name, predicted_R, predicted_period, predicted_m_n)
    ("A_anchor_R12_odds5_63", tuple(range(5, HIGH_N + 1, 2)), "canonical",
     Fraction(1, 2), 8, (1, 4)),
    ("P1a_R14_mixed_2_5", (2, 5), "canonical", Fraction(1, 4), 12, (1, 3)),
    ("P1b_R14_mixed_2_4_5_7", (2, 4, 5, 7), "canonical", Fraction(1, 4), 12, (1, 3)),
    ("P1c_R14_mixed_6elem", (2, 4, 6, 5, 7, 9), "canonical", Fraction(1, 4), 12, (1, 3)),
    ("P2_R34_M2_odds5_63", tuple(range(5, HIGH_N + 1, 2)), "m2",
     Fraction(3, 4), 6, (1, 6)),
    ("P3_R0_fundamental", (1,), "canonical", Fraction(0, 1), None, None),
)

THETA_FNS: dict[str, Callable] = {"canonical": canonical_theta, "m2": m2_theta}


def run_case(case_id, support, fn_name, pred_r, pred_period, pred_mn, sp) -> dict[str, Any]:
    theta_fn = THETA_FNS[fn_name]
    a, b = make_pair(support, sp)
    initial_a, initial_b = a.copy(), b.copy()
    initial_norm = toy.pair_hermitian_norm(a, b)
    theta0 = theta_fn(a, b, sp)
    r0 = math.sin(theta0) ** 2

    residuals = []
    accumulated = 0.0
    first_return = None
    acc_at_return = None
    for collision in range(1, LOCK_COLLISIONS + 1):
        th = theta_fn(a, b, sp)
        a, b = toy.rotate_ab(a, b, th)
        accumulated += th
        diff = float(
            np.vdot(a - initial_a, a - initial_a).real
            + np.vdot(b - initial_b, b - initial_b).real
        )
        residual = math.sqrt(max(diff, 0.0) / initial_norm)
        residuals.append(residual)
        if first_return is None and residual <= RECURRENCE_TOLERANCE:
            first_return = collision
            acc_at_return = accumulated

    row: dict[str, Any] = {
        "case": case_id,
        "readout": fn_name,
        "support_size": len(support),
        "predicted_R": str(pred_r),
        "measured_R": r0,
        "R_err": abs(r0 - float(pred_r)),
        "R_pass": abs(r0 - float(pred_r)) <= R_TOLERANCE,
        "predicted_period": pred_period,
        "first_return": first_return,
        "period_pass": (first_return == pred_period) if pred_period else None,
        "max_residual": max(residuals),
    }
    if pred_period is None:
        row["invariant_pass"] = max(residuals) <= 1.0e-12
    if first_return and acc_at_return is not None:
        theta_mean = acc_at_return / first_return
        x = 0.5 - theta_mean / math.pi
        frac = Fraction(x).limit_denominator(RATIONAL_MAX_DEN)
        row.update(
            {
                "reconstructed_m_n": [frac.numerator, frac.denominator],
                "predicted_m_n": list(pred_mn) if pred_mn else None,
                "m_n_pass": (frac.numerator, frac.denominator) == pred_mn if pred_mn else None,
                "wind_per_period": acc_at_return / (2.0 * math.pi),
            }
        )
    row["_residuals"] = residuals
    return row


def main() -> None:
    params = base.Params(high_n=HIGH_N, recursive_collision_count=LOCK_COLLISIONS)
    sp = base.build_source_params(params)

    rows = []
    residual_store = {}
    for case in CASES:
        row = run_case(*case, sp)
        residual_store[row["case"]] = row.pop("_residuals")
        rows.append(row)
        print(
            f"{row['case']:>24}: R={row['measured_R']:.12f}"
            f" (pred {row['predicted_R']}, err {row['R_err']:.1e},"
            f" {'PASS' if row['R_pass'] else 'FAIL'})"
            f" period={row['first_return']} (pred {row['predicted_period']},"
            f" {row.get('period_pass')})"
            f" (m,n)={row.get('reconstructed_m_n')} pred {row.get('predicted_m_n')}"
        )

    all_pass = all(
        row["R_pass"]
        and (row["period_pass"] is not False)
        and (row.get("m_n_pass") is not False)
        and (row.get("invariant_pass") is not False)
        for row in rows
    )
    print(f"\nALL NIVEN POINT PREDICTIONS: {'PASS' if all_pass else 'FAIL'}")

    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    for case_id, residuals in residual_store.items():
        if case_id.startswith("P3"):
            continue
        ax.semilogy(range(1, len(residuals) + 1), residuals, label=case_id, linewidth=1.0)
    ax.axhline(RECURRENCE_TOLERANCE, color="0.5", linestyle=":", label="tolerance")
    for p in (6, 8, 12):
        ax.axvline(p, color="0.85", linewidth=0.6)
    ax.set_xlabel("collision")
    ax.set_ylabel("pair return residual")
    ax.set_title(
        "Niven-point counting locks: periods 12 (R=1/4), 6 (R=3/4), 8 (R=1/2 anchor)"
    )
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)
    figure_names = []
    for ext in ("png", "svg"):
        path = HERE / f"niven_points_lock_pre_v1.{ext}"
        fig.savefig(path, dpi=160)
        figure_names.append(path.name)
    plt.close(fig)

    csv_path = HERE / "niven_points_lock_pre_rows_v1.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as h:
        fieldnames = sorted({k for row in rows for k in row})
        w = csv.DictWriter(h, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    payload = {
        "experiment": "niven_points_lock_pre_v1",
        "core_runner": {"path": TOY_RUNNER_PATH.name, "sha256": toy.sha256(TOY_RUNNER_PATH)},
        "predictions_fixed_before_measurement": True,
        "niven_intersection": ["0", "1/4", "1/2", "3/4", "1"],
        "note_R1": "R=1 は A 基本波が全マスクでボゾン的なため数え上げ構成では到達不能（振幅調整族での周期4回帰は R sweep 実験で確認済み）",
        "rows": rows,
        "all_pass": all_pass,
        "figures": figure_names,
    }
    (HERE / "niven_points_lock_pre_result_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
