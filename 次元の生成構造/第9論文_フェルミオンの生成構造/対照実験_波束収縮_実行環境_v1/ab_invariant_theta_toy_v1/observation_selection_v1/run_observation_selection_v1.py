#!/usr/bin/env python3
"""観測選択実験 v1（E6）

仮説（主張4の実験的裏付け・初の直接検証）:
    「観測（初期状態への回帰の確認）は、厳密回帰する有限位数状態に
     アンサンブルを指数集中させる」

方法:
    - 状態族: 従来の振幅族（A固定・B振幅スケール）61点＋厳密根状態2点
      （R=1/2、R_{124,23}。逆算探索で構成——探索は初期状態のみ、前進は無変更）
    - 各状態を内生 theta 力学で J 衝突発展させ、回帰忠実度
        F_J = |<psi_0|psi_J>|^2 / |<psi_0|psi_0>|^2
      を測る（psi = (a,b) 対、ペアエルミート内積）
    - 観測を n 回繰り返した後の生存重み F_J^n が、厳密回帰状態に
      指数集中することを確認する
    - 観測時計 J = 8 と J = 248 の二種を測り、どの状態族が選ばれるかが
      J の約数構造に依存すること（選択が観測時計に継承されること）を記録する

予言（測定前に固定）:
    P1: F_8 = 1（機械精度）は R=1/2 状態のみ。F_248 = 1 は R=1/2（8|248）と
        R_{124,23}（周期248）の両方。
    P2: 非回帰状態の F は 1 より真に小さく、n=20 観測後の生存重みは
        回帰状態に対し指数的に劣る（集中率 > 10^3）。
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
TOY_RUNNER_PATH = HERE.parent / "run_ab_invariant_theta_toy_v1.py"
SEARCH_PATH = HERE.parent / "inverse_initial_conditions_v1" / "search_initial_conditions_and_plot_v1.py"

spec = importlib.util.spec_from_file_location("search_for_obs_sel_v1", SEARCH_PATH)
search = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = search
spec.loader.exec_module(search)
toy = search.toy
base = toy.base
plt = base.plt

J_CLOCKS = (8, 248)
N_OBS = 20
AMPLITUDES = np.geomspace(0.05, 5.0, 61)
R_ROOT = math.cos(math.pi * 23.0 / 124.0) ** 2


def pair_fidelity(a0, b0, a, b) -> float:
    ov = complex(np.vdot(a0, a) + np.vdot(b0, b))
    n0 = float(np.vdot(a0, a0).real + np.vdot(b0, b0).real)
    return abs(ov) ** 2 / n0**2


def evolve_fidelities(a0, b0, sp, j_max: int, checkpoints) -> dict[int, float]:
    a, b = a0.copy(), b0.copy()
    out = {}
    for j in range(1, j_max + 1):
        readout = toy.theta_from_ab(a, b, sp)
        a, b = toy.rotate_ab(a, b, readout.theta)
        if j in checkpoints:
            out[j] = pair_fidelity(a0, b0, a, b)
    return out


def main() -> None:
    params = base.Params(high_n=63, recursive_collision_count=max(J_CLOCKS))
    sp = base.build_source_params(params)
    a_template, b_template, _ = search.make_unit_templates(sp)

    states = []
    for amp in AMPLITUDES:
        states.append((f"amp_{amp:.4f}", a_template.copy(), float(amp) * b_template, "family"))
    for label, target in (("root_R050", 0.5), ("root_R124_23", R_ROOT)):
        result = search.search_initial_b_amplitude(
            target, a_template, b_template, sp, tolerance=1.0e-15
        )
        states.append(
            (label, a_template.copy(), result.initial_b_amplitude * b_template, "exact_root")
        )

    rows = []
    for label, a0, b0, kind in states:
        fids = evolve_fidelities(a0, b0, sp, max(J_CLOCKS), set(J_CLOCKS))
        r0 = toy.theta_from_ab(a0, b0, sp).reflection_rate
        rows.append(
            {
                "state": label, "kind": kind, "R_initial": r0,
                "F_8": fids[8], "F_248": fids[248],
                "survival_8_n20": fids[8] ** N_OBS,
                "survival_248_n20": fids[248] ** N_OBS,
            }
        )

    # ---- 判定 ----
    root50 = next(r for r in rows if r["state"] == "root_R050")
    root124 = next(r for r in rows if r["state"] == "root_R124_23")
    family = [r for r in rows if r["kind"] == "family"]

    p1a = abs(root50["F_8"] - 1.0) <= 1e-10 and abs(root50["F_248"] - 1.0) <= 1e-10
    p1b = abs(root124["F_248"] - 1.0) <= 1e-10
    p1c = root124["F_8"] < 1.0 - 1e-6  # 124根は周期248なのでJ=8では戻らない
    best_family_248 = max(r["survival_248_n20"] for r in family)
    concentration = root124["survival_248_n20"] / max(best_family_248, 1e-300)
    p2 = concentration > 1e3

    print(f"P1: F_8(R=1/2)={root50['F_8']:.15f}, F_248(R=1/2)={root50['F_248']:.15f} -> {p1a}")
    print(f"P1: F_248(root124)={root124['F_248']:.15f}, F_8(root124)={root124['F_8']:.6f} -> {p1b and p1c}")
    print(f"P2: best family survival(J=248,n=20)={best_family_248:.3e}")
    print(f"P2: concentration ratio root124/family = {concentration:.3e} -> {p2}")
    all_pass = p1a and p1b and p1c and p2
    print(f"OBSERVATION-SELECTION MECHANISM: {'PASS' if all_pass else 'FAIL'}")
    print("note: J=248 は 8|248 により R=1/2 族も選ぶ＝選択は観測時計の約数構造に継承される")

    with (HERE / "observation_selection_rows_v1.csv").open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)
    fam_r = [r["R_initial"] for r in family]
    ax1.semilogy(fam_r, [max(r["F_248"], 1e-16) for r in family], ".", markersize=4, label="family F_248")
    ax1.semilogy([root124["R_initial"]], [root124["F_248"]], "r*", markersize=14, label="root 124 (F=1)")
    ax1.semilogy([root50["R_initial"]], [root50["F_248"]], "g*", markersize=14, label="root 1/2 (F=1)")
    ax1.set_xlabel("initial R"); ax1.set_ylabel("return fidelity F_248")
    ax1.set_title("return fidelity landscape (J=248)")
    ax1.legend(fontsize=8); ax1.grid(alpha=0.3)
    ax2.semilogy(fam_r, [max(r["survival_248_n20"], 1e-300) for r in family], ".", markersize=4)
    ax2.semilogy([root124["R_initial"]], [root124["survival_248_n20"]], "r*", markersize=14)
    ax2.semilogy([root50["R_initial"]], [root50["survival_248_n20"]], "g*", markersize=14)
    ax2.set_xlabel("initial R"); ax2.set_ylabel(f"survival after {N_OBS} observations")
    ax2.set_title("iterated observation concentrates on exact-recurrence states")
    ax2.grid(alpha=0.3)
    for ext in ("png", "svg"):
        fig.savefig(HERE / f"observation_selection_v1.{ext}", dpi=160)
    plt.close(fig)

    payload = {
        "experiment": "observation_selection_v1",
        "core_runner": {"path": TOY_RUNNER_PATH.name, "sha256": toy.sha256(TOY_RUNNER_PATH)},
        "predictions_fixed_before_measurement": True,
        "verdicts": {
            "P1_root50_both_clocks": p1a,
            "P1_root124_J248_only": p1b and p1c,
            "P2_concentration_ratio": concentration,
            "all_pass": all_pass,
        },
        "note": "選択は観測時計 J の約数構造に継承される（J=248 は周期8族も選ぶ）。「なぜ124か」は「観測時計がなぜ248可換か」に変換される——蒸し返しでなく問題の座標変換",
        "rows": rows,
    }
    (HERE / "observation_selection_result_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
