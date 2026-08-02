#!/usr/bin/env python3
"""勘定則検証と数え上げ→ロック直結・予備テスト v1

三部構成:
  Part 1: 二記憶バッテリーの5点から解読した勘定則
              R_eq(S) = (2|S| - c(S)) / (4|S|)
              c(S) = sum over k in S of { k=1: 2, k=3: 1, even k: 2, odd k>=5: 0 }
          を、バッテリーで未使用の留め置き支持で予言検証する。
          予言値はコード内に定数として測定前に書き込む（仮説→反証の作法）。

  Part 2: {1,3}-free 支持（c=0）では勘定則が R = 1/2 を予言する。
          R = 1/2 = cos^2(pi/4) は「有理数 R 族」と「有理角 theta 族」の
          唯一の非自明な交点であり、位数 n=4（実回転周期 8）の厳密ロック。
          等重み・振幅探索ゼロの純粋な数え上げ構成だけで、
          8衝突厳密回帰が出ることを確認する（数え上げ→ロックの初の直結）。

  Part 3: ロック状態に対し ord/wind 読出しを実装検証する。
          ord  = 最初の厳密回帰衝突数 P（予言: 8）
          wind = 一周期の累積回転角 / 2pi（予言: 1）
          theta_meas = 累積角/P から x = 1/2 - theta_meas/pi を作り、
          有理数再構成で (m, n) を得る（予言: (1, 4)）。
          この読出し機構が、後に B63 で (23, 124) を狙う本丸実験の道具になる。

設計境界:
  - theta_from_ab・rotate_ab は無変更で呼ぶだけ
  - 状態は explicit_packet_case / make_case_state の等重みパケットのみ
  - **振幅探索（search_initial_b_amplitude）は一切使用しない**
    （数え上げ構成だけでロックに到達することが Part 2 の主張だから）
"""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
TOY_RUNNER_PATH = HERE.parent / "run_ab_invariant_theta_toy_v1.py"

HIGH_N = 63
FORMULA_TOLERANCE = 1.0e-12
RECURRENCE_TOLERANCE = 1.0e-8
LOCK_COLLISIONS = 100
RATIONAL_MAX_DENOMINATOR = 2000


def load_toy_module() -> Any:
    spec = importlib.util.spec_from_file_location(
        "ab_theta_toy_for_counting_rule_bridge_v1",
        TOY_RUNNER_PATH,
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


def excluded_half_bins(support: tuple[int, ...]) -> int:
    """勘定則の c(S)。バッテリー5点から解読した仮説（本実験の検証対象）。"""
    c = 0
    for k in support:
        if k == 1:
            c += 2
        elif k == 3:
            c += 1
        elif k % 2 == 0:
            c += 2
    return c


def predicted_r(support: tuple[int, ...]) -> Fraction:
    size = len(support)
    return Fraction(2 * size - excluded_half_bins(support), 4 * size)


# ---- Part 1 の留め置き支持（バッテリー未使用）と、測定前の予言値 ----
HELD_OUT_SUPPORTS: dict[str, tuple[int, ...]] = {
    "H_pair_1_3": (1, 3),
    "H_single_5": (5,),
    "H_odds_5_to_63": tuple(range(5, HIGH_N + 1, 2)),
    "H_triple_1_5_7": (1, 5, 7),
    "H_pair_3_7": (3, 7),
    "H_odds_1_to_15": tuple(range(1, 16, 2)),
    "H_mixed_1_2_5": (1, 2, 5),
}
PREDICTIONS: dict[str, Fraction] = {
    name: predicted_r(support) for name, support in HELD_OUT_SUPPORTS.items()
}

# ---- Part 2 のロック候補（{1,3}-free、c=0 → R=1/2 予言）----
LOCK_SUPPORTS: dict[str, tuple[int, ...]] = {
    "L_odds_5_to_63": tuple(range(5, HIGH_N + 1, 2)),
    "L_single_5": (5,),
    "L_step4_from_5": tuple(range(5, HIGH_N + 1, 4)),
}
PREDICTED_LOCK = {
    "R": Fraction(1, 2),
    "theta_over_pi": Fraction(1, 4),
    "real_rotation_period": 8,
    "wind_per_period": 1,
    "root_convention_m_n": (1, 4),
}


def unit_norm(vector: np.ndarray) -> np.ndarray:
    norm = math.sqrt(float(np.vdot(vector, vector).real))
    if norm <= 0.0:
        raise ValueError("zero-norm state")
    return vector / norm


def make_states(
    support: tuple[int, ...], source_params: Any
) -> tuple[np.ndarray, np.ndarray]:
    case = base.explicit_packet_case(
        mode="counting_" + "_".join(str(k) for k in support[:4]) + f"_len{len(support)}",
        packet_a=(1,),
        packet_b=support,
    )
    a = unit_norm(base.make_case_state(source_params, case, "A", hair_enabled=True))
    b = unit_norm(base.make_case_state(source_params, case, "B", hair_enabled=True))
    return a, b


def measure_r(support: tuple[int, ...], source_params: Any) -> float:
    a, b = make_states(support, source_params)
    return toy.theta_from_ab(a, b, source_params).reflection_rate


def lock_run(
    label: str,
    support: tuple[int, ...],
    source_params: Any,
) -> dict[str, Any]:
    a, b = make_states(support, source_params)
    initial_a = a.copy()
    initial_b = b.copy()
    initial_norm = toy.pair_hermitian_norm(a, b)
    initial_readout = toy.theta_from_ab(a, b, source_params)

    residuals: list[float] = []
    accumulated_theta = 0.0
    accumulated_theta_first_period: float | None = None
    first_return: int | None = None
    for collision in range(1, LOCK_COLLISIONS + 1):
        readout = toy.theta_from_ab(a, b, source_params)
        a, b = toy.rotate_ab(a, b, readout.theta)
        accumulated_theta += readout.theta
        diff = float(
            np.vdot(a - initial_a, a - initial_a).real
            + np.vdot(b - initial_b, b - initial_b).real
        )
        residual = math.sqrt(max(diff, 0.0) / initial_norm)
        residuals.append(residual)
        if first_return is None and residual <= RECURRENCE_TOLERANCE:
            first_return = collision
            accumulated_theta_first_period = accumulated_theta

    period_multiples_ok = None
    if first_return is not None:
        period_multiples_ok = all(
            residuals[j - 1] <= RECURRENCE_TOLERANCE
            for j in range(first_return, LOCK_COLLISIONS + 1, first_return)
        )

    result: dict[str, Any] = {
        "label": label,
        "support_size": len(support),
        "measured_R": initial_readout.reflection_rate,
        "abs_R_minus_half": abs(initial_readout.reflection_rate - 0.5),
        "first_exact_return_collision": first_return,
        "predicted_period": PREDICTED_LOCK["real_rotation_period"],
        "period_matches_prediction": first_return
        == PREDICTED_LOCK["real_rotation_period"],
        "all_period_multiples_return": period_multiples_ok,
        "min_residual": min(residuals),
        "residual_at_8": residuals[7] if len(residuals) >= 8 else None,
        "amplitude_search_used": False,
    }

    if first_return is not None and accumulated_theta_first_period is not None:
        theta_meas = accumulated_theta_first_period / first_return
        wind = accumulated_theta_first_period / (2.0 * math.pi)
        x = 0.5 - theta_meas / math.pi
        frac = Fraction(x).limit_denominator(RATIONAL_MAX_DENOMINATOR)
        result.update(
            {
                "theta_measured": theta_meas,
                "wind_per_period_measured": wind,
                "wind_matches_prediction": abs(
                    wind - PREDICTED_LOCK["wind_per_period"]
                )
                <= 1.0e-9,
                "x_half_minus_theta_over_pi": x,
                "reconstructed_m": frac.numerator,
                "reconstructed_n": frac.denominator,
                "m_n_matches_prediction": (frac.numerator, frac.denominator)
                == PREDICTED_LOCK["root_convention_m_n"],
                "rational_reconstruction_error": abs(x - float(frac)),
            }
        )
    result["_residuals"] = residuals
    return result


def main() -> None:
    params = base.Params(high_n=HIGH_N, recursive_collision_count=LOCK_COLLISIONS)
    source_params = base.build_source_params(params)

    # ---- Part 1: 留め置き予言検証 ----
    part1_rows: list[dict[str, Any]] = []
    for name, support in HELD_OUT_SUPPORTS.items():
        prediction = PREDICTIONS[name]
        measured = measure_r(support, source_params)
        error = abs(measured - float(prediction))
        part1_rows.append(
            {
                "label": name,
                "support": " ".join(str(k) for k in support[:8])
                + ("..." if len(support) > 8 else ""),
                "support_size": len(support),
                "c_excluded_half_bins": excluded_half_bins(support),
                "predicted_R": str(prediction),
                "predicted_R_float": float(prediction),
                "measured_R": measured,
                "abs_error": error,
                "verdict": "PASS" if error <= FORMULA_TOLERANCE else "FAIL",
            }
        )
    part1_pass = all(row["verdict"] == "PASS" for row in part1_rows)

    # ---- Part 2 & 3: 数え上げ→ロック直結と ord/wind 読出し ----
    lock_rows: list[dict[str, Any]] = []
    residual_series: dict[str, list[float]] = {}
    for label, support in LOCK_SUPPORTS.items():
        row = lock_run(label, support, source_params)
        residual_series[label] = row.pop("_residuals")
        lock_rows.append(row)
    part2_pass = all(
        row["abs_R_minus_half"] <= FORMULA_TOLERANCE
        and row["period_matches_prediction"]
        and row["all_period_multiples_return"]
        for row in lock_rows
    )
    part3_pass = all(
        row.get("m_n_matches_prediction", False)
        and row.get("wind_matches_prediction", False)
        for row in lock_rows
    )

    # ---- 図 ----
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    labels = [row["label"] for row in part1_rows]
    ax1.scatter(
        range(len(part1_rows)),
        [row["predicted_R_float"] for row in part1_rows],
        marker="o",
        s=70,
        facecolors="none",
        edgecolors="tab:blue",
        label="predicted (written before measurement)",
    )
    ax1.scatter(
        range(len(part1_rows)),
        [row["measured_R"] for row in part1_rows],
        marker="x",
        s=45,
        color="tab:red",
        label="measured",
    )
    ax1.set_xticks(range(len(part1_rows)))
    ax1.set_xticklabels(labels, rotation=30, ha="right", fontsize=7)
    ax1.set_ylabel("R")
    ax1.set_title("Part 1: held-out counting-rule predictions")
    ax1.legend(fontsize=7)
    ax1.grid(alpha=0.3)

    for label, residuals in residual_series.items():
        ax2.semilogy(range(1, len(residuals) + 1), residuals, label=label, linewidth=1.0)
    ax2.axhline(RECURRENCE_TOLERANCE, color="0.5", linestyle=":", label="tolerance")
    for j in range(8, LOCK_COLLISIONS + 1, 8):
        ax2.axvline(j, color="0.85", linewidth=0.6)
    ax2.set_xlabel("collision")
    ax2.set_ylabel("pair return residual")
    ax2.set_title("Part 2: counting-built {1,3}-free states, period-8 exact returns")
    ax2.legend(fontsize=7)
    ax2.grid(alpha=0.3)
    fig.suptitle("Counting rule verification and counting-to-lock bridge (pre v1)")
    figure_names = []
    for ext in ("png", "svg"):
        path = HERE / f"counting_rule_lock_bridge_pre_v1.{ext}"
        fig.savefig(path, dpi=160)
        figure_names.append(path.name)
    plt.close(fig)

    # ---- 保存 ----
    csv_path = HERE / "counting_rule_lock_bridge_pre_rows_v1.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(part1_rows[0]))
        writer.writeheader()
        writer.writerows(part1_rows)

    payload = {
        "experiment": "counting_rule_lock_bridge_pre_v1",
        "design_boundary": {
            "theta_readout_modified": False,
            "forward_scattering_modified": False,
            "amplitude_search_used": False,
            "states": "equal-weight packet constructions only (pure counting)",
        },
        "counting_rule_hypothesis": {
            "formula": "R_eq(S) = (2|S| - c(S)) / (4|S|)",
            "c_rule": "k=1 -> 2, k=3 -> 1, even k -> 2, odd k>=5 -> 0",
            "decoded_from": "two_memory_readout_battery_pre_v1 (5 support classes)",
        },
        "core_runner": {
            "path": TOY_RUNNER_PATH.name,
            "sha256": toy.sha256(TOY_RUNNER_PATH),
        },
        "part1_held_out_verification": {
            "rows": part1_rows,
            "tolerance": FORMULA_TOLERANCE,
            "verdict": "PASS" if part1_pass else "FAIL",
        },
        "part2_counting_to_lock": {
            "predicted": {
                "R": str(PREDICTED_LOCK["R"]),
                "real_rotation_period": PREDICTED_LOCK["real_rotation_period"],
            },
            "rows": lock_rows,
            "verdict": "PASS" if part2_pass else "FAIL",
        },
        "part3_ord_wind_readout": {
            "predicted_m_n": list(PREDICTED_LOCK["root_convention_m_n"]),
            "predicted_wind": PREDICTED_LOCK["wind_per_period"],
            "verdict": "PASS" if part3_pass else "FAIL",
        },
        "figures": figure_names,
    }
    (HERE / "counting_rule_lock_bridge_pre_result_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print(f"Part1 held-out formula: {'PASS' if part1_pass else 'FAIL'}")
    for row in part1_rows:
        print(
            f"  {row['label']:>18}: pred={row['predicted_R']:>6}"
            f" meas={row['measured_R']:.12f} err={row['abs_error']:.3e} {row['verdict']}"
        )
    print(f"Part2 counting-to-lock: {'PASS' if part2_pass else 'FAIL'}")
    for row in lock_rows:
        print(
            f"  {row['label']:>16}: |R-1/2|={row['abs_R_minus_half']:.3e}"
            f" first_return={row['first_exact_return_collision']}"
            f" residual@8={row['residual_at_8']:.3e}"
        )
    print(f"Part3 ord/wind readout: {'PASS' if part3_pass else 'FAIL'}")
    for row in lock_rows:
        print(
            f"  {row['label']:>16}: wind={row.get('wind_per_period_measured')}"
            f" (m,n)=({row.get('reconstructed_m')},{row.get('reconstructed_n')})"
        )


if __name__ == "__main__":
    main()
