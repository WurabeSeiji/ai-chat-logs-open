#!/usr/bin/env python3
"""R質量うなり二腕干渉計・予備テスト v1

目的:
    R_read の異なる二状態（二腕）を同一の無変更前進処理で別々に発展させ、
    腕間の重なり積分 O_ij(j) = <psi_i(j)|psi_j(j)> の位相が衝突番号 j とともに
    蓄積するか（うなり）を測る配管検証。うなり傾き Omega が Delta R に
    比例するなら質量はスケール円の巻き面（第4共役対）。Omega=0 が機械精度で
    出るなら零次 no-go の候補（逆二乗論文 v1→v2 と同じ手順で次段へ）。

設計境界（既存実験と同一）:
    - 散乱本体・theta 読出しは無変更（theta は毎衝突、現在のABから再導出）
    - 目標 R は初期状態探索にのみ使用し、前進処理へは渡さない
    - 重ね合わせは物理的に走らせない（theta 読出しが状態依存で非線形のため）。
      二腕を別々に発展させ、重なり積分を後処理で追う二腕干渉計方式。

対照テスト（シリーズ内完結の再現性規約）:
    実行冒頭で R124_23 の 248 衝突厳密回帰（residual <= 1e-8）を再現し、
    既存 exact_finite_order_roots_longrun_v1 の公表値と一致することを確認
    してから本測定に進む。失敗した場合は本測定を実行しない。
"""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
INITIAL_SEARCH_PATH = HERE.parent / "search_initial_conditions_and_plot_v1.py"

MAX_COLLISION_PRE = 300
SNAPSHOT_COLLISIONS = (0, 1, 2, 10, 100, 248, 300)
SEARCH_TOLERANCE = 1.0e-15
INVARIANT_TOLERANCE = 1.0e-10
ANCHOR_RECURRENCE_TOLERANCE = 1.0e-8
OVERLAP_FIT_MIN_ABS = 1.0e-3


def load_initial_search_module() -> Any:
    spec = importlib.util.spec_from_file_location(
        "initial_state_search_for_r_mass_beat_pre_v1",
        INITIAL_SEARCH_PATH,
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load initial-state search: {INITIAL_SEARCH_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


initial_search = load_initial_search_module()
toy = initial_search.toy
plt = initial_search.plt


ARM_SPECS = (
    # (arm_id, target_R, note)
    ("R030", 0.30, "generic mid-low R"),
    ("R050", 0.50, "exact 8-collision closure point"),
    ("R050_dup", 0.50, "duplicate of R050: null-test arm (expects phi=0, |O|=1)"),
    ("R070", 0.70, "generic mid-high R"),
    ("R124_23", math.cos(math.pi * 23.0 / 124.0) ** 2, "exact root, period 248"),
    ("R620_117", math.cos(math.pi * 117.0 / 620.0) ** 2, "exact root, period 1240"),
)


def prepare_arm(
    target_r: float,
    a_template: np.ndarray,
    b_template: np.ndarray,
    source_params: Any,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """目標Rを初期探索にだけ使い、単位ペアノルムに規格化した腕状態を返す。"""
    search_result = initial_search.search_initial_b_amplitude(
        target_r,
        a_template,
        b_template,
        source_params,
        tolerance=SEARCH_TOLERANCE,
    )
    a = a_template.copy()
    b = search_result.initial_b_amplitude * b_template
    norm = toy.pair_hermitian_norm(a, b)
    scale = 1.0 / math.sqrt(norm)
    a = a * scale
    b = b * scale
    readout = toy.theta_from_ab(a, b, source_params)
    info = {
        "search": asdict(search_result),
        "normalized_pair_norm": toy.pair_hermitian_norm(a, b),
        "initial_R_after_normalization": readout.reflection_rate,
        "initial_theta_after_normalization": readout.theta,
    }
    return a, b, info


def evolve_one_collision(
    a: np.ndarray,
    b: np.ndarray,
    source_params: Any,
) -> tuple[np.ndarray, np.ndarray, Any]:
    """無変更の前進処理: theta を現在のABから再導出して1回転する。"""
    readout = toy.theta_from_ab(a, b, source_params)
    a2, b2 = toy.rotate_ab(a, b, readout.theta)
    return a2, b2, readout


def pair_overlap(
    a1: np.ndarray,
    b1: np.ndarray,
    a2: np.ndarray,
    b2: np.ndarray,
) -> complex:
    """ペアエルミート内積 <psi_1|psi_2> = <a1|a2> + <b1|b2>。"""
    return complex(np.vdot(a1, a2) + np.vdot(b1, b2))


def anchor_control_test(
    a_template: np.ndarray,
    b_template: np.ndarray,
    source_params: Any,
) -> dict[str, Any]:
    """R124_23 の248衝突厳密回帰を再現する対照テスト。"""
    target_r = math.cos(math.pi * 23.0 / 124.0) ** 2
    a, b, info = prepare_arm(target_r, a_template, b_template, source_params)
    initial_a = a.copy()
    initial_b = b.copy()
    initial_norm = toy.pair_hermitian_norm(a, b)
    initial_r = toy.theta_from_ab(a, b, source_params).reflection_rate
    max_r_drift = 0.0
    for _ in range(248):
        a, b, readout = evolve_one_collision(a, b, source_params)
        max_r_drift = max(max_r_drift, abs(readout.reflection_rate - initial_r))
    diff = float(
        np.vdot(a - initial_a, a - initial_a).real
        + np.vdot(b - initial_b, b - initial_b).real
    )
    residual = math.sqrt(max(diff, 0.0) / initial_norm)
    passed = (
        residual <= ANCHOR_RECURRENCE_TOLERANCE
        and max_r_drift <= INVARIANT_TOLERANCE
    )
    return {
        "anchor": "R124_23 period-248 exact return",
        "target_R": target_r,
        "achieved_initial_R": info["initial_R_after_normalization"],
        "return_residual_at_248": residual,
        "max_R_drift_through_248": max_r_drift,
        "recurrence_tolerance": ANCHOR_RECURRENCE_TOLERANCE,
        "invariant_tolerance": INVARIANT_TOLERANCE,
        "verdict": "PASS" if passed else "FAIL",
    }


def fit_phase_slope(
    collisions: np.ndarray,
    unwrapped_phase: np.ndarray,
    overlap_abs: np.ndarray,
) -> dict[str, Any]:
    """|O| が閾値以上の点だけで unwrap 済み位相を一次フィットする。"""
    mask = overlap_abs >= OVERLAP_FIT_MIN_ABS
    n_used = int(np.count_nonzero(mask))
    if n_used < 3:
        return {
            "slope_rad_per_collision": None,
            "intercept_rad": None,
            "linear_fit_r_squared": None,
            "points_used": n_used,
            "note": "insufficient points above |O| threshold",
        }
    x = collisions[mask].astype(float)
    y = unwrapped_phase[mask]
    coeffs = np.polyfit(x, y, 1)
    predicted = np.polyval(coeffs, x)
    ss_res = float(np.sum((y - predicted) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else 1.0
    return {
        "slope_rad_per_collision": float(coeffs[0]),
        "intercept_rad": float(coeffs[1]),
        "linear_fit_r_squared": r_squared,
        "points_used": n_used,
        "note": "",
    }


def make_figures(
    pair_ids: list[tuple[str, str]],
    collisions: np.ndarray,
    phase_by_pair: dict[tuple[str, str], np.ndarray],
    abs_by_pair: dict[tuple[str, str], np.ndarray],
    fits: dict[tuple[str, str], dict[str, Any]],
    delta_r_by_pair: dict[tuple[str, str], float],
) -> list[str]:
    figure_names: list[str] = []

    fig, (ax_phase, ax_abs) = plt.subplots(
        2, 1, figsize=(12, 9), sharex=True, constrained_layout=True
    )
    for pair in pair_ids:
        label = f"{pair[0]}-{pair[1]} (dR={delta_r_by_pair[pair]:.3f})"
        ax_phase.plot(collisions, phase_by_pair[pair], label=label, linewidth=1.0)
        ax_abs.plot(collisions, abs_by_pair[pair], label=label, linewidth=1.0)
    ax_phase.set_ylabel("unwrapped arg <psi_i|psi_j> [rad]")
    ax_phase.set_title("Two-arm overlap phase vs collision (pre-test)")
    ax_phase.grid(alpha=0.3)
    ax_phase.legend(fontsize=6)
    ax_abs.set_ylabel("|<psi_i|psi_j>|")
    ax_abs.set_xlabel("collision")
    ax_abs.set_yscale("log")
    ax_abs.grid(alpha=0.3)
    fig.suptitle("R mass-beat two-arm interferometer, preliminary v1", fontsize=13)
    for stem in ("overlap_phase_and_abs_pre_v1",):
        for ext in ("png", "svg"):
            path = HERE / f"{stem}.{ext}"
            fig.savefig(path, dpi=160)
            figure_names.append(path.name)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    xs, ys = [], []
    for pair in pair_ids:
        fit = fits[pair]
        if fit["slope_rad_per_collision"] is None:
            continue
        xs.append(delta_r_by_pair[pair])
        ys.append(fit["slope_rad_per_collision"])
        ax.annotate(
            f"{pair[0]}-{pair[1]}",
            (delta_r_by_pair[pair], fit["slope_rad_per_collision"]),
            fontsize=6,
        )
    ax.scatter(xs, ys, color="tab:purple")
    ax.set_xlabel("|Delta R| between arms")
    ax.set_ylabel("fitted phase slope Omega [rad/collision]")
    ax.set_title("Beat slope vs Delta R (pre-test)")
    ax.grid(alpha=0.3)
    for ext in ("png", "svg"):
        path = HERE / f"beat_slope_vs_delta_R_pre_v1.{ext}"
        fig.savefig(path, dpi=160)
        figure_names.append(path.name)
    plt.close(fig)
    return figure_names


def main() -> None:
    params = toy.base.Params(
        high_n=63,
        recursive_collision_count=MAX_COLLISION_PRE,
    )
    source_params = toy.base.build_source_params(params)
    a_template, b_template, case = initial_search.make_unit_templates(source_params)

    anchor = anchor_control_test(a_template, b_template, source_params)
    print(
        "anchor control:",
        anchor["verdict"],
        f"residual={anchor['return_residual_at_248']:.3e}",
        f"R_drift={anchor['max_R_drift_through_248']:.3e}",
    )
    if anchor["verdict"] != "PASS":
        payload = {"experiment": "R_mass_beat_two_arm_pre_v1", "anchor_control": anchor}
        (HERE / "R_mass_beat_two_arm_pre_result_v1.json").write_text(
            json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        raise SystemExit("anchor control FAILED: not proceeding to measurement")

    arms: dict[str, dict[str, Any]] = {}
    states: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for arm_id, target_r, note in ARM_SPECS:
        if arm_id == "R050_dup":
            base_a, base_b = states["R050"]
            states[arm_id] = (base_a.copy(), base_b.copy())
            arms[arm_id] = dict(arms["R050"])
            arms[arm_id]["note"] = note
            continue
        a, b, info = prepare_arm(target_r, a_template, b_template, source_params)
        states[arm_id] = (a, b)
        arms[arm_id] = {"target_R": target_r, "note": note, **info}

    arm_ids = [spec[0] for spec in ARM_SPECS]
    pair_ids = [
        (arm_ids[i], arm_ids[j])
        for i in range(len(arm_ids))
        for j in range(i + 1, len(arm_ids))
    ]
    delta_r_by_pair = {
        (p, q): abs(
            arms[p]["initial_R_after_normalization"]
            - arms[q]["initial_R_after_normalization"]
        )
        for (p, q) in pair_ids
    }

    collisions = np.arange(MAX_COLLISION_PRE + 1)
    raw_overlap: dict[tuple[str, str], list[complex]] = {p: [] for p in pair_ids}
    r_drift: dict[str, float] = {arm_id: 0.0 for arm_id in arm_ids}
    initial_r = {
        arm_id: arms[arm_id]["initial_R_after_normalization"] for arm_id in arm_ids
    }
    snapshots: dict[str, dict[int, np.ndarray]] = {"collisions": {}}

    for j in range(MAX_COLLISION_PRE + 1):
        for p, q in pair_ids:
            a1, b1 = states[p]
            a2, b2 = states[q]
            raw_overlap[(p, q)].append(pair_overlap(a1, b1, a2, b2))
        if j in SNAPSHOT_COLLISIONS:
            for arm_id in arm_ids:
                a, b = states[arm_id]
                snapshots.setdefault(arm_id, {})[j] = np.stack([a, b])
        if j < MAX_COLLISION_PRE:
            for arm_id in arm_ids:
                a, b = states[arm_id]
                a2, b2, readout = evolve_one_collision(a, b, source_params)
                states[arm_id] = (a2, b2)
                r_drift[arm_id] = max(
                    r_drift[arm_id],
                    abs(readout.reflection_rate - initial_r[arm_id]),
                )

    phase_by_pair: dict[tuple[str, str], np.ndarray] = {}
    abs_by_pair: dict[tuple[str, str], np.ndarray] = {}
    fits: dict[tuple[str, str], dict[str, Any]] = {}
    for pair in pair_ids:
        series = np.asarray(raw_overlap[pair])
        abs_by_pair[pair] = np.abs(series)
        phase_by_pair[pair] = np.unwrap(np.angle(series))
        fits[pair] = fit_phase_slope(
            collisions, phase_by_pair[pair], abs_by_pair[pair]
        )

    null_pair = ("R050", "R050_dup")
    null_max_phase = float(np.max(np.abs(phase_by_pair[null_pair])))
    null_min_abs = float(np.min(abs_by_pair[null_pair]))
    null_pass = null_max_phase <= 1.0e-12 and abs(null_min_abs - 1.0) <= 1.0e-12

    rows_path = HERE / "R_mass_beat_two_arm_pre_rows_v1.csv"
    with rows_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["arm_i", "arm_j", "collision", "overlap_abs", "overlap_phase_unwrapped"]
        )
        for pair in pair_ids:
            for j in range(MAX_COLLISION_PRE + 1):
                writer.writerow(
                    [
                        pair[0],
                        pair[1],
                        j,
                        f"{abs_by_pair[pair][j]:.17g}",
                        f"{phase_by_pair[pair][j]:.17g}",
                    ]
                )

    npz_payload: dict[str, np.ndarray] = {}
    for arm_id in arm_ids:
        for j, stacked in snapshots.get(arm_id, {}).items():
            npz_payload[f"{arm_id}_collision_{j}"] = stacked
    np.savez_compressed(HERE / "arm_state_snapshots_pre_v1.npz", **npz_payload)

    figure_names = make_figures(
        pair_ids, collisions, phase_by_pair, abs_by_pair, fits, delta_r_by_pair
    )

    payload = {
        "experiment": "R_mass_beat_two_arm_pre_v1",
        "design_boundary": {
            "target_R_used_only_in": "standalone initial-state error minimization",
            "core_theta_readout_modified": False,
            "forward_scattering_external_R_or_theta": False,
            "superposition_evolved_physically": False,
            "protocol": "two separately evolved arms, post-hoc overlap phase",
        },
        "case": asdict(case),
        "conditions": {
            "max_collision": MAX_COLLISION_PRE,
            "snapshot_collisions": list(SNAPSHOT_COLLISIONS),
            "search_tolerance": SEARCH_TOLERANCE,
            "invariant_tolerance": INVARIANT_TOLERANCE,
            "overlap_fit_min_abs": OVERLAP_FIT_MIN_ABS,
            "arms_normalized_to_unit_pair_norm": True,
        },
        "core_runner": {
            "path": str(
                initial_search.TOY_RUNNER_PATH.relative_to(initial_search.TOY_DIR)
            ),
            "sha256": toy.sha256(initial_search.TOY_RUNNER_PATH),
        },
        "anchor_control": anchor,
        "arms": arms,
        "max_R_drift_by_arm": r_drift,
        "null_pair_test": {
            "pair": list(null_pair),
            "max_abs_phase": null_max_phase,
            "min_overlap_abs": null_min_abs,
            "verdict": "PASS" if null_pass else "FAIL",
        },
        "pair_results": [
            {
                "arm_i": pair[0],
                "arm_j": pair[1],
                "delta_R": delta_r_by_pair[pair],
                "min_overlap_abs": float(np.min(abs_by_pair[pair])),
                "mean_overlap_abs": float(np.mean(abs_by_pair[pair])),
                "final_unwrapped_phase": float(phase_by_pair[pair][-1]),
                **fits[pair],
            }
            for pair in pair_ids
        ],
        "figures": figure_names,
    }
    (HERE / "R_mass_beat_two_arm_pre_result_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print(f"null pair {null_pair}:", "PASS" if null_pass else "FAIL")
    for pair in pair_ids:
        fit = fits[pair]
        slope = fit["slope_rad_per_collision"]
        slope_text = "n/a" if slope is None else f"{slope:.6e}"
        print(
            f"pair={pair[0]}-{pair[1]}",
            f"dR={delta_r_by_pair[pair]:.6f}",
            f"slope={slope_text}",
            f"R2={fit['linear_fit_r_squared']}",
            f"min|O|={np.min(abs_by_pair[pair]):.3e}",
        )


if __name__ == "__main__":
    main()
