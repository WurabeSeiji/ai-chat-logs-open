#!/usr/bin/env python3
"""規約監査・予備テスト v1

目的:
    本日の主張群を、トイの実装規約（フェルミオンマスクの定義）を変えて
    再測定し、規約非依存（物理級）と規約依存（簿記級）に分類する。

監査対象の主張:
    C1 等重み点の R は数え上げ（ビン勘定）で予言できる  [メタ勘定則]
    C2 等重み点の R は小分母有理数になる                [有理数性]
    C3 振幅変形で R が動く                              [可視層汚染の普遍性]
    C4 「全ビン通過」倍音だけの支持は振幅探索ゼロで R=1/2 に着地し
       周期8で厳密ロックする                            [ロック橋頭堡]
    C5 パリティ反転マスクではロック支持が偶数倍音側へ移る
       （ボゾン/フェルミオンの役割は規約ラベル）        [役割の規約相対性]

方法:
    - theta 読出しをマスク関数でパラメタ化した変種を実装し、
      正準マスクで toy.theta_from_ab と厳密一致することを対照テストで確認
      （コピー→対照テスト→拡張の規約に従う）
    - 各倍音の「ビン表」（どの周波数ビンにどれだけの割合で乗るか）を
      単一倍音状態から実測し、予言はビン表×マスクの勘定だけで立てる
    - 前進動力学は rotate_ab 無変更、theta のみ変種読出し
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
PREDICTION_TOLERANCE = 1.0e-12
RATIONAL_TOLERANCE = 1.0e-12
RATIONAL_MAX_DEN = 1000
AMPLITUDE_SPREAD_MIN = 1.0e-3
LOCK_COLLISIONS = 24
RECURRENCE_TOLERANCE = 1.0e-8
BIN_POWER_FLOOR = 1.0e-12


def load_toy_module() -> Any:
    spec = importlib.util.spec_from_file_location(
        "ab_theta_toy_for_convention_audit_v1", TOY_RUNNER_PATH
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

ODD_FULL = tuple(range(1, HIGH_N + 1, 2))
EVEN_FULL = tuple(range(2, HIGH_N, 2))


MASKS: dict[str, Callable[[np.ndarray], np.ndarray]] = {
    "M0_canonical_even_ge4": lambda f: (f >= 4) & (f % 2 == 0),
    "M1_even_ge6": lambda f: (f >= 6) & (f % 2 == 0),
    "M2_even_ge2": lambda f: (f >= 2) & (f % 2 == 0),
    "M3_any_ge4": lambda f: f >= 4,
    "M4_parity_flip_odd_ge3": lambda f: (f >= 3) & (f % 2 == 1),
}


def variant_theta(
    a: np.ndarray,
    b: np.ndarray,
    source_params: Any,
    mask_fn: Callable[[np.ndarray], np.ndarray],
) -> tuple[float, float]:
    """toy.theta_from_ab のマスク部だけをパラメタ化した変種。

    正準マスクでは元実装と厳密一致することを anchor_test で確認する。
    """
    frequencies, pair_power = toy.combined_chi_power(a, b, source_params)
    abs_frequency = np.abs(frequencies)
    fermionic_mask = mask_fn(abs_frequency)
    total_power = float(np.sum(pair_power))
    fermionic_power = float(np.sum(pair_power[fermionic_mask]))
    numerical_floor = 1024.0 * np.finfo(float).eps * max(total_power, 1.0)
    if abs(fermionic_power) <= numerical_floor:
        fermionic_power = 0.0
    fermionic_power = min(max(fermionic_power, 0.0), total_power)
    bosonic_power = max(total_power - fermionic_power, 0.0)
    theta = math.atan2(math.sqrt(fermionic_power), math.sqrt(bosonic_power))
    return theta, math.sin(theta) ** 2


def unit_norm(vector: np.ndarray) -> np.ndarray:
    norm = math.sqrt(float(np.vdot(vector, vector).real))
    if norm <= 0.0:
        raise ValueError("zero-norm state")
    return vector / norm


def make_state(packet: tuple[int, ...], which: str, source_params: Any) -> np.ndarray:
    case = base.explicit_packet_case(
        mode=f"audit_{which}_len{len(packet)}_{packet[0]}",
        packet_a=(1,),
        packet_b=packet,
    )
    return unit_norm(
        base.make_case_state(source_params, case, which, hair_enabled=True)
    )


def bin_profile(state: np.ndarray, source_params: Any) -> dict[int, float]:
    """状態単独の周波数ビン別パワー割合（合計1）。"""
    zero = np.zeros_like(state)
    frequencies, power = toy.combined_chi_power(state, zero, source_params)
    total = float(np.sum(power))
    profile: dict[int, float] = {}
    for f, p in zip(frequencies, power):
        if p / total > BIN_POWER_FLOOR:
            profile[int(f)] = profile.get(int(f), 0.0) + float(p / total)
    return profile


def predicted_r_from_bins(
    support: tuple[int, ...],
    weights: dict[int, float] | None,
    a_profile: dict[int, float],
    b_profiles: dict[int, dict[int, float]],
    mask_fn: Callable[[np.ndarray], np.ndarray],
) -> float:
    """ビン表×マスクだけから R を予言（状態の再FFTはしない）。"""
    def masked_fraction(profile: dict[int, float]) -> float:
        keys = np.abs(np.asarray(list(profile), dtype=int))
        vals = np.asarray(list(profile.values()))
        return float(np.sum(vals[mask_fn(keys)]))

    if weights is None:
        w2 = {h: 1.0 / len(support) for h in support}
    else:
        norm = sum(w * w for w in weights.values())
        w2 = {h: weights[h] * weights[h] / norm for h in support}
    p_f = masked_fraction(a_profile)  # A unit power
    for h in support:
        p_f += w2[h] * masked_fraction(b_profiles[h])
    return p_f / 2.0


def lock_supports_from_bins(
    pool: tuple[int, ...],
    a_profile: dict[int, float],
    b_profiles: dict[int, dict[int, float]],
    mask_fn: Callable[[np.ndarray], np.ndarray],
) -> tuple[tuple[int, ...], bool]:
    """A が全ボゾンかつ「全ビン通過」倍音の集合（あれば）を返す。"""
    def masked_fraction(profile: dict[int, float]) -> float:
        keys = np.abs(np.asarray(list(profile), dtype=int))
        vals = np.asarray(list(profile.values()))
        return float(np.sum(vals[mask_fn(keys)]))

    a_bosonic = masked_fraction(a_profile) <= 1.0e-12
    fully_fermionic = tuple(
        h for h in pool if abs(masked_fraction(b_profiles[h]) - 1.0) <= 1.0e-12
    )
    return fully_fermionic, a_bosonic


def lock_run(
    support: tuple[int, ...],
    source_params: Any,
    mask_fn: Callable[[np.ndarray], np.ndarray],
) -> dict[str, Any]:
    a = make_state(support, "A", source_params)
    b = make_state(support, "B", source_params)
    initial_a, initial_b = a.copy(), b.copy()
    initial_norm = toy.pair_hermitian_norm(a, b)
    _, r0 = variant_theta(a, b, source_params, mask_fn)
    accumulated = 0.0
    first_return = None
    acc_at_return = None
    for collision in range(1, LOCK_COLLISIONS + 1):
        theta, _ = variant_theta(a, b, source_params, mask_fn)
        a, b = toy.rotate_ab(a, b, theta)
        accumulated += theta
        diff = float(
            np.vdot(a - initial_a, a - initial_a).real
            + np.vdot(b - initial_b, b - initial_b).real
        )
        residual = math.sqrt(max(diff, 0.0) / initial_norm)
        if first_return is None and residual <= RECURRENCE_TOLERANCE:
            first_return = collision
            acc_at_return = accumulated
    row: dict[str, Any] = {
        "support_size": len(support),
        "support_head": " ".join(str(k) for k in support[:5]),
        "measured_R": r0,
        "abs_R_minus_half": abs(r0 - 0.5),
        "first_return": first_return,
        "period8": first_return == 8,
    }
    if first_return and acc_at_return is not None:
        theta_mean = acc_at_return / first_return
        x = 0.5 - theta_mean / math.pi
        frac = Fraction(x).limit_denominator(RATIONAL_MAX_DEN)
        row.update(
            {
                "wind": acc_at_return / (2.0 * math.pi),
                "m_n": [frac.numerator, frac.denominator],
                "m_n_is_1_4": (frac.numerator, frac.denominator) == (1, 4),
            }
        )
    return row


def main() -> None:
    params = base.Params(high_n=HIGH_N, recursive_collision_count=LOCK_COLLISIONS)
    source_params = base.build_source_params(params)

    # ---- 対照テスト: 正準マスク変種 = 元実装 ----
    anchor_states = [
        (make_state(ODD_FULL, "A", source_params), make_state(ODD_FULL, "B", source_params)),
        (make_state((1, 5, 7), "A", source_params), make_state((1, 5, 7), "B", source_params)),
    ]
    anchor_max = 0.0
    for a, b in anchor_states:
        _, r_variant = variant_theta(
            a, b, source_params, MASKS["M0_canonical_even_ge4"]
        )
        r_original = toy.theta_from_ab(a, b, source_params).reflection_rate
        anchor_max = max(anchor_max, abs(r_variant - r_original))
    anchor_pass = anchor_max <= 1.0e-15
    print(f"anchor (variant==original under M0): {'PASS' if anchor_pass else 'FAIL'} ({anchor_max:.2e})")
    if not anchor_pass:
        raise SystemExit("anchor failed")

    # ---- ビン表の実測 ----
    a_profile = bin_profile(make_state((1,), "A", source_params), source_params)
    pool = tuple(sorted((*ODD_FULL, *EVEN_FULL)))
    b_profiles = {
        h: bin_profile(make_state((h,), "B", source_params), source_params)
        for h in pool
    }

    inverse_weights = {h: 1.0 / h for h in ODD_FULL}
    test_states: dict[str, tuple[tuple[int, ...], dict[int, float] | None]] = {
        "B63_equal": (ODD_FULL, None),
        "B63_inverse_k": (ODD_FULL, inverse_weights),
        "odds_5_to_63_equal": (tuple(range(5, HIGH_N + 1, 2)), None),
        "evens_equal": (EVEN_FULL, None),
    }

    prebuilt: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for label, (support, weights) in test_states.items():
        a = make_state(support, "A", source_params)
        if weights is None:
            b = make_state(support, "B", source_params)
        else:
            singles = {h: make_state((h,), "B", source_params) for h in support}
            b = unit_norm(sum(weights[h] * singles[h] for h in support))
        prebuilt[label] = (a, b)

    rows: list[dict[str, Any]] = []
    lock_rows: list[dict[str, Any]] = []
    for mask_name, mask_fn in MASKS.items():
        r_by_state: dict[str, float] = {}
        for label, (support, weights) in test_states.items():
            a, b = prebuilt[label]
            _, r_meas = variant_theta(a, b, source_params, mask_fn)
            r_pred = predicted_r_from_bins(
                support, weights, a_profile, b_profiles, mask_fn
            )
            frac = Fraction(r_meas).limit_denominator(RATIONAL_MAX_DEN)
            is_rational = abs(r_meas - float(frac)) <= RATIONAL_TOLERANCE
            rows.append(
                {
                    "mask": mask_name,
                    "state": label,
                    "R_measured": r_meas,
                    "R_predicted_bin_counting": r_pred,
                    "prediction_error": abs(r_meas - r_pred),
                    "C1_counting_prediction": "PASS"
                    if abs(r_meas - r_pred) <= PREDICTION_TOLERANCE
                    else "FAIL",
                    "nearest_rational": f"{frac.numerator}/{frac.denominator}",
                    "C2_rational_at_equal_weight": (
                        "PASS" if is_rational else ("FAIL" if weights is None else "n/a")
                    )
                    if weights is None
                    else ("expected_irrational_" + ("yes" if not is_rational else "NO")),
                }
            )
            r_by_state[label] = r_meas
        spread = abs(r_by_state["B63_equal"] - r_by_state["B63_inverse_k"])
        rows.append(
            {
                "mask": mask_name,
                "state": "C3_amplitude_spread(B63 equal vs 1/k)",
                "R_measured": spread,
                "R_predicted_bin_counting": None,
                "prediction_error": None,
                "C1_counting_prediction": "",
                "nearest_rational": "",
                "C2_rational_at_equal_weight": "PASS" if spread >= AMPLITUDE_SPREAD_MIN else "FAIL",
            }
        )

        # ---- ロック橋頭堡（C4/C5）----
        fully_f, a_bosonic = lock_supports_from_bins(
            pool, a_profile, b_profiles, mask_fn
        )
        entry: dict[str, Any] = {
            "mask": mask_name,
            "A_fully_bosonic": a_bosonic,
            "fully_fermionic_harmonics_head": " ".join(str(k) for k in fully_f[:6]),
            "fully_fermionic_count": len(fully_f),
            "parity_of_lock_harmonics": (
                "odd" if fully_f and fully_f[0] % 2 == 1 else ("even" if fully_f else "none")
            ),
        }
        if fully_f and a_bosonic:
            result = lock_run(fully_f, source_params, mask_fn)
            entry.update(result)
        else:
            entry["note"] = "no counting-lock support under this mask"
        lock_rows.append(entry)

    csv_path = HERE / "convention_audit_pre_rows_v1.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    payload = {
        "experiment": "convention_audit_pre_v1",
        "audited_claims": {
            "C1": "equal-weight R predictable from measured bin table + mask (meta counting rule)",
            "C2": "equal-weight R is small-denominator rational under every mask",
            "C3": "amplitude deformation moves R under every mask",
            "C4": "fully-transmitting supports lock at R=1/2, period 8, no tuning",
            "C5": "parity-flipped mask moves the lock to even harmonics",
        },
        "anchor_variant_equals_original": anchor_max,
        "core_runner": {
            "path": TOY_RUNNER_PATH.name,
            "sha256": toy.sha256(TOY_RUNNER_PATH),
        },
        "readout_rows": rows,
        "lock_rows": lock_rows,
    }
    (HERE / "convention_audit_pre_result_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    for row in rows:
        if "spread" in str(row["state"]):
            print(f"{row['mask']:<26} C3 spread={row['R_measured']:.6f} {row['C2_rational_at_equal_weight']}")
        else:
            print(
                f"{row['mask']:<26} {row['state']:<20}"
                f" R={row['R_measured']:.12f} pred_err={row['prediction_error']:.2e}"
                f" {row['C1_counting_prediction']:<4} rat={row['nearest_rational']:<8}"
                f" {row['C2_rational_at_equal_weight']}"
            )
    print("---- locks ----")
    for entry in lock_rows:
        print(
            f"{entry['mask']:<26} parity={entry.get('parity_of_lock_harmonics'):<5}"
            f" n={entry.get('fully_fermionic_count')}"
            f" R-1/2={entry.get('abs_R_minus_half', 'n/a')}"
            f" period8={entry.get('period8', 'n/a')}"
            f" (m,n)={entry.get('m_n', 'n/a')}"
        )


if __name__ == "__main__":
    main()
