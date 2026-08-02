#!/usr/bin/env python3
"""二記憶判別バッテリー・予備テスト v1

目的:
    波形の二種類の記憶
        (1) 大きさの記憶 = ノルム（振幅二乗総和）
        (2) 形の記憶     = 倍音の存在集合（support）
    に対し、現行の各読出しがどちらに依存するかを分類する。

    条件A: 存在集合を固定し、振幅分布だけを変える（同ペアノルム）
        -> 構造層読出しなら不変のはず。現行のパワー比 theta 読出しが
           どれだけ変動するかを定量化する（可視層汚染の測定）
    条件B: ペアノルムを固定し、存在集合だけを変える
        -> 構造層読出しは離散的に変わるはず。数え上げ仮説
           n_counting = 4 * N_rel（N_rel = Aと非共通の倍音数）の予言値も併記

設計境界:
    - theta_from_ab・rotate_ab は無変更で呼ぶだけ（複製・改変しない）
    - 状態構築は既存の explicit_packet_case / make_case_state のみを使用
      （単一倍音状態の重み付き和で振幅分布を構成。散乱本体に手を入れない）

対照テスト（シリーズ内完結の再現性規約）:
    等重み単一倍音和で構成した B63 状態の R 読出しが、既存テンプレート
    （make_case_state による等振幅 B63）の R 読出しと一致するかを最初に
    確認する。不一致の場合はその値を記録し、等重み和を条件Aの基線とする。
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

HIGH_N = 63
EVOLUTION_COLLISIONS = 100
INVARIANT_TOLERANCE = 1.0e-10
CONSTRUCTION_MATCH_TOLERANCE = 1.0e-12
SUPPORT_POWER_FLOOR_RATIO = 1.0e-12
RANDOM_SEED = 7


def load_toy_module() -> Any:
    spec = importlib.util.spec_from_file_location(
        "ab_theta_toy_for_two_memory_battery_v1",
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

ODD_FULL = tuple(range(1, HIGH_N + 1, 2))


def make_single_harmonic_b(k: int, source_params: Any) -> np.ndarray:
    case = base.explicit_packet_case(
        mode=f"single_harmonic_{k}",
        packet_a=(1,),
        packet_b=(k,),
    )
    return base.make_case_state(source_params, case, "B", hair_enabled=True)


def make_packet_b(packet: tuple[int, ...], source_params: Any) -> np.ndarray:
    case = base.explicit_packet_case(
        mode="packet_" + "_".join(str(k) for k in packet[:4]) + f"_len{len(packet)}",
        packet_a=(1,),
        packet_b=packet,
    )
    return base.make_case_state(source_params, case, "B", hair_enabled=True)


def make_a_state(source_params: Any) -> np.ndarray:
    case = base.explicit_packet_case(
        mode="fundamental_a",
        packet_a=(1,),
        packet_b=(1,),
    )
    return base.make_case_state(source_params, case, "A", hair_enabled=True)


def unit_norm(vector: np.ndarray) -> np.ndarray:
    norm = math.sqrt(float(np.vdot(vector, vector).real))
    if norm <= 0.0:
        raise ValueError("zero-norm state")
    return vector / norm


def weighted_sum_b(
    weights: dict[int, float],
    singles: dict[int, np.ndarray],
) -> np.ndarray:
    total = None
    for k, w in weights.items():
        term = w * singles[k]
        total = term if total is None else total + term
    return unit_norm(total)


def readout_row(
    label: str,
    condition: str,
    support: tuple[int, ...],
    a: np.ndarray,
    b: np.ndarray,
    source_params: Any,
    metric_context: Any,
) -> dict[str, Any]:
    readout = toy.theta_from_ab(a, b, source_params)
    pair_norm = toy.pair_hermitian_norm(a, b)
    metrics_b = toy.state_metrics(b, metric_context)
    n_rel = len([k for k in support if k != 1])
    return {
        "label": label,
        "condition": condition,
        "support_size": len(support),
        "support_min": min(support),
        "support_max": max(support),
        "N_rel_noncommon_with_A": n_rel,
        "n_counting_hypothesis_4Nrel": 4 * n_rel,
        "pair_norm": pair_norm,
        "R_power_ratio_readout": readout.reflection_rate,
        "theta": readout.theta,
        "P_fermionic": readout.fermionic_relation_power,
        "P_bosonic": readout.bosonic_relation_power,
        "L_B": metrics_b["L"],
        "N_eff_B": metrics_b["N_eff"],
    }


def evolution_invariance_check(
    a: np.ndarray,
    b: np.ndarray,
    source_params: Any,
) -> dict[str, float]:
    a_run = a.copy()
    b_run = b.copy()
    initial = toy.theta_from_ab(a_run, b_run, source_params).reflection_rate
    max_drift = 0.0
    for _ in range(EVOLUTION_COLLISIONS):
        readout = toy.theta_from_ab(a_run, b_run, source_params)
        max_drift = max(max_drift, abs(readout.reflection_rate - initial))
        a_run, b_run = toy.rotate_ab(a_run, b_run, readout.theta)
    return {
        "collisions": EVOLUTION_COLLISIONS,
        "max_R_drift": max_drift,
        "verdict": "PASS" if max_drift <= INVARIANT_TOLERANCE else "CHECK",
    }


def main() -> None:
    params = base.Params(high_n=HIGH_N, recursive_collision_count=EVOLUTION_COLLISIONS)
    source_params = base.build_source_params(params)
    metric_context = base.MetricContext(source_params)

    a_state = unit_norm(make_a_state(source_params))
    singles = {k: unit_norm(make_single_harmonic_b(k, source_params)) for k in ODD_FULL}

    rows: list[dict[str, Any]] = []

    # ---- 対照テスト: 等重み和 vs 既存テンプレート（等振幅 B63）----
    b_template = unit_norm(make_packet_b(ODD_FULL, source_params))
    b_equal_sum = weighted_sum_b({k: 1.0 for k in ODD_FULL}, singles)
    r_template = toy.theta_from_ab(a_state, b_template, source_params).reflection_rate
    r_equal_sum = toy.theta_from_ab(a_state, b_equal_sum, source_params).reflection_rate
    construction_diff = abs(r_template - r_equal_sum)
    construction_match = construction_diff <= CONSTRUCTION_MATCH_TOLERANCE
    overlap = abs(complex(np.vdot(b_template, b_equal_sum)))
    anchor = {
        "anchor": "equal-weight single-harmonic sum vs make_case_state B63 template",
        "R_template": r_template,
        "R_equal_weight_sum": r_equal_sum,
        "abs_R_difference": construction_diff,
        "state_overlap_abs": overlap,
        "verdict": "MATCH" if construction_match else "RECORDED_MISMATCH",
        "note": (
            "一致すれば構成法は等価。僅差なら等重み和を条件A基線として使用"
        ),
    }

    # ---- 条件A: 同存在集合・異振幅分布（Bは単位ノルム、Aは共通）----
    rng = np.random.default_rng(RANDOM_SEED)
    count = len(ODD_FULL)
    profiles: dict[str, dict[int, float]] = {
        "A_equal": {k: 1.0 for k in ODD_FULL},
        "A_linear_decay": {
            k: 1.0 - 0.8 * (idx / (count - 1)) for idx, k in enumerate(ODD_FULL)
        },
        "A_inverse_k": {k: 1.0 / k for k in ODD_FULL},
        "A_exp_decay": {k: math.exp(-k / 21.0) for k in ODD_FULL},
        "A_random_seed7": {
            k: float(w) for k, w in zip(ODD_FULL, rng.uniform(0.2, 1.0, count))
        },
    }
    for label, weights in profiles.items():
        b_state = weighted_sum_b(weights, singles)
        rows.append(
            readout_row(
                label,
                "A_same_support_diff_amplitude",
                ODD_FULL,
                a_state,
                b_state,
                source_params,
                metric_context,
            )
        )

    # ---- 条件B: 同ノルム・異存在集合（等重み、Bは単位ノルム）----
    supports: dict[str, tuple[int, ...]] = {
        "B_odds_full_63": ODD_FULL,
        "B_odds_upto_31": tuple(range(1, 32, 2)),
        "B_odds_sparse_step4": tuple(range(1, HIGH_N + 1, 4)),
        "B_odds_full_minus_3": tuple(k for k in ODD_FULL if k != 3),
        "B_evens_control": (1, *tuple(range(2, HIGH_N, 2))),
    }
    for label, support in supports.items():
        needed = {k: 1.0 for k in support}
        extra_singles = {
            k: (singles[k] if k in singles else unit_norm(make_single_harmonic_b(k, source_params)))
            for k in support
        }
        b_state = weighted_sum_b(needed, extra_singles)
        rows.append(
            readout_row(
                label,
                "B_same_norm_diff_support",
                support,
                a_state,
                b_state,
                source_params,
                metric_context,
            )
        )

    # ---- 発展下の不変性（基線＋各条件から1本ずつ）----
    evolution_checks: dict[str, dict[str, float]] = {}
    for label, weights_or_support in (
        ("A_equal", profiles["A_equal"]),
        ("A_inverse_k", profiles["A_inverse_k"]),
        ("B_odds_upto_31", supports["B_odds_upto_31"]),
    ):
        if label.startswith("A"):
            b_state = weighted_sum_b(weights_or_support, singles)
        else:
            ss = {
                k: (singles[k] if k in singles else unit_norm(make_single_harmonic_b(k, source_params)))
                for k in weights_or_support
            }
            b_state = weighted_sum_b({k: 1.0 for k in weights_or_support}, ss)
        evolution_checks[label] = evolution_invariance_check(
            a_state, b_state, source_params
        )

    # ---- 集計: 条件Aの R 変動幅（可視層汚染の定量）----
    a_rows = [r for r in rows if r["condition"].startswith("A_")]
    b_rows = [r for r in rows if r["condition"].startswith("B_")]
    a_r_values = [r["R_power_ratio_readout"] for r in a_rows]
    baseline_r = next(
        r["R_power_ratio_readout"] for r in a_rows if r["label"] == "A_equal"
    )
    summary = {
        "condition_A_R_min": min(a_r_values),
        "condition_A_R_max": max(a_r_values),
        "condition_A_R_spread": max(a_r_values) - min(a_r_values),
        "condition_A_max_abs_deviation_from_equal": max(
            abs(v - baseline_r) for v in a_r_values
        ),
        "condition_A_structural_readout_requirement": (
            "spread は構造層読出しなら 0 のはず。非零なら現行パワー比読出しは"
            "可視層（振幅）に依存しており、数え上げ型読出しへの置換が必要"
        ),
        "condition_B_R_by_support": {
            r["label"]: r["R_power_ratio_readout"] for r in b_rows
        },
        "condition_B_counting_hypothesis": {
            r["label"]: {
                "N_rel": r["N_rel_noncommon_with_A"],
                "n_counting_4Nrel": r["n_counting_hypothesis_4Nrel"],
            }
            for r in b_rows
        },
    }

    csv_path = HERE / "two_memory_readout_battery_pre_rows_v1.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    a_labels = [r["label"] for r in a_rows]
    ax_a.bar(range(len(a_rows)), a_r_values, color="tab:blue")
    ax_a.axhline(baseline_r, color="0.4", linestyle=":", label="A_equal baseline")
    ax_a.set_xticks(range(len(a_rows)))
    ax_a.set_xticklabels(a_labels, rotation=30, ha="right", fontsize=7)
    ax_a.set_ylabel("R (power-ratio readout)")
    ax_a.set_title("Condition A: same support, different amplitude profile")
    ax_a.legend(fontsize=7)
    ax_a.grid(alpha=0.3)

    b_labels = [r["label"] for r in b_rows]
    b_values = [r["R_power_ratio_readout"] for r in b_rows]
    ax_b.bar(range(len(b_rows)), b_values, color="tab:orange")
    ax_b.set_xticks(range(len(b_rows)))
    ax_b.set_xticklabels(b_labels, rotation=30, ha="right", fontsize=7)
    ax_b.set_title("Condition B: same norm, different support")
    ax_b.grid(alpha=0.3)
    fig.suptitle("Two-memory readout battery (pre-test v1)", fontsize=13)
    figure_names = []
    for ext in ("png", "svg"):
        path = HERE / f"two_memory_readout_battery_pre_v1.{ext}"
        fig.savefig(path, dpi=160)
        figure_names.append(path.name)
    plt.close(fig)

    payload = {
        "experiment": "two_memory_readout_battery_pre_v1",
        "design_boundary": {
            "theta_readout_modified": False,
            "forward_scattering_modified": False,
            "state_construction": (
                "existing explicit_packet_case / make_case_state only; "
                "amplitude profiles via weighted sums of single-harmonic states"
            ),
        },
        "conditions": {
            "high_n": HIGH_N,
            "evolution_collisions": EVOLUTION_COLLISIONS,
            "invariant_tolerance": INVARIANT_TOLERANCE,
            "random_seed": RANDOM_SEED,
            "A_state": "fundamental (1,) with hair, unit norm",
            "B_states": "unit norm; pair norm equal across all states",
        },
        "core_runner": {
            "path": TOY_RUNNER_PATH.name,
            "sha256": toy.sha256(TOY_RUNNER_PATH),
        },
        "anchor_construction_test": anchor,
        "rows": rows,
        "evolution_invariance": evolution_checks,
        "summary": summary,
        "figures": figure_names,
    }
    (HERE / "two_memory_readout_battery_pre_result_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print("anchor:", anchor["verdict"], f"dR={anchor['abs_R_difference']:.3e}")
    print(
        "condition A spread:",
        f"{summary['condition_A_R_spread']:.6e}",
        f"(baseline R={baseline_r:.12f})",
    )
    for r in a_rows:
        print(f"  {r['label']:>16}: R={r['R_power_ratio_readout']:.12f}")
    print("condition B:")
    for r in b_rows:
        print(
            f"  {r['label']:>22}: R={r['R_power_ratio_readout']:.12f}",
            f"N_rel={r['N_rel_noncommon_with_A']}",
            f"4*N_rel={r['n_counting_hypothesis_4Nrel']}",
        )
    for label, check in evolution_checks.items():
        print(f"evolution {label}: {check['verdict']} drift={check['max_R_drift']:.3e}")


if __name__ == "__main__":
    main()
