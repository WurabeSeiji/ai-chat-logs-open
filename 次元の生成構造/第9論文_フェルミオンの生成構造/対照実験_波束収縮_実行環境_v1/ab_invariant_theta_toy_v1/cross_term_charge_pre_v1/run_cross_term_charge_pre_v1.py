#!/usr/bin/env python3
"""交差項＝電荷文法・予備テスト v1

前実験（mixed_coupling_pre_v1）の帰結:
    対角セクターは平均合成（重力的：加法・符号なし）。積構造と符号は
    干渉交差項 2Re(A*B) に住むはず、という予言を動力学で検証する。

予言（反証可能）:
    B のフェルミオン成分だけに位相 phi を付けると、衝突1回のチャネル
    ノルム流 dN_B は
        dN_B(phi) = offset + M cos(phi + delta)
    と変調し、変調振幅は
        M = kappa * sin(2 theta) * sqrt(f_A f_B)      [振幅積 = 頂点積]
    で、kappa と delta はグリッド全体で一定。さらに phi の半回転で
    流れの向きが反転する（符号＝相対位相：引力/斥力の類似）。
    f_A=0 または f_B=0 では M=0（中性閉包は符号を感じない＝遮蔽類似）。

方法:
    A(f_A) = sqrt(1-f_A) A1 + sqrt(f_A) A5
    B(f_B, phi) = sqrt(1-f_B) B1 + e^{i phi} sqrt(f_B) B7
    （A5 と B7 はビン+6 を共有 → マスク内交差項の担い手）
    グリッド f in {0, 0.2, ..., 1.0}^2 × phi 8点。衝突は無変更の
    theta_from_ab + rotate_ab を1回適用し、dN_B を測る。
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
PHASES = tuple(2.0 * math.pi * k / 8 for k in range(8))
KAPPA_CONSTANCY_TOL = 1.0e-10
ZERO_MOD_TOL = 1.0e-13


def load_toy_module() -> Any:
    spec = importlib.util.spec_from_file_location(
        "ab_theta_toy_for_cross_term_charge_v1", TOY_RUNNER_PATH
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
    return v / math.sqrt(float(np.vdot(v, v).real))


def make_single(k: int, which: str, source_params: Any) -> np.ndarray:
    case = base.explicit_packet_case(
        mode=f"xterm_{which}_{k}", packet_a=(k,), packet_b=(k,)
    )
    return unit_norm(base.make_case_state(source_params, case, which, hair_enabled=True))


def fit_cosine(phases: np.ndarray, values: np.ndarray) -> tuple[float, float, float]:
    """values ~ offset + C cos(phi) + S sin(phi) を最小二乗で解き、
    (offset, M=sqrt(C^2+S^2), delta) を返す。"""
    design = np.column_stack(
        [np.ones_like(phases), np.cos(phases), np.sin(phases)]
    )
    coef, *_ = np.linalg.lstsq(design, values, rcond=None)
    offset, c, s = coef
    return float(offset), float(math.hypot(c, s)), float(math.atan2(-s, c))


def main() -> None:
    params = base.Params(high_n=63, recursive_collision_count=1)
    source_params = base.build_source_params(params)

    a1 = make_single(1, "A", source_params)
    a5 = make_single(5, "A", source_params)
    b1 = make_single(1, "B", source_params)
    b7 = make_single(7, "B", source_params)

    rows: list[dict[str, Any]] = []
    fits: list[dict[str, Any]] = []
    for f_a in F_GRID:
        a0 = unit_norm(math.sqrt(1.0 - f_a) * a1 + math.sqrt(f_a) * a5)
        for f_b in F_GRID:
            flows = []
            thetas = []
            for phi in PHASES:
                b0 = unit_norm(
                    math.sqrt(1.0 - f_b) * b1
                    + np.exp(1j * phi) * math.sqrt(f_b) * b7
                )
                readout = toy.theta_from_ab(a0, b0, source_params)
                a2, b2 = toy.rotate_ab(a0.copy(), b0.copy(), readout.theta)
                dn_b = float(np.vdot(b2, b2).real) - 1.0
                flows.append(dn_b)
                thetas.append(readout.theta)
                rows.append(
                    {
                        "f_A": f_a,
                        "f_B": f_b,
                        "phi": phi,
                        "theta": readout.theta,
                        "dN_B_one_collision": dn_b,
                    }
                )
            offset, m_mod, delta = fit_cosine(
                np.asarray(PHASES), np.asarray(flows)
            )
            theta_mean = float(np.mean(thetas))
            sin2t = math.sin(2.0 * theta_mean)
            product_amp = math.sqrt(f_a * f_b)
            kappa = (
                m_mod / (abs(sin2t) * product_amp)
                if abs(sin2t) * product_amp > 1e-15
                else None
            )
            fits.append(
                {
                    "f_A": f_a,
                    "f_B": f_b,
                    "offset": offset,
                    "modulation_M": m_mod,
                    "delta": delta,
                    "sin_2theta": sin2t,
                    "sqrt_fA_fB": product_amp,
                    "kappa": kappa,
                    "flow_at_delta": offset + m_mod,
                    "flow_at_delta_plus_pi": offset - m_mod,
                    "sign_flip": (offset + m_mod) * (offset - m_mod) < 0.0,
                }
            )

    # ---- 判定 ----
    charged = [f for f in fits if f["kappa"] is not None]
    neutral = [f for f in fits if f["kappa"] is None]
    kappas = np.asarray([f["kappa"] for f in charged])
    deltas = np.asarray([f["delta"] for f in charged])
    kappa_spread = float(np.max(kappas) - np.min(kappas))
    delta_spread = float(np.max(deltas) - np.min(deltas))
    neutral_max_mod = max(f["modulation_M"] for f in neutral) if neutral else 0.0
    n_sign_flips = sum(1 for f in charged if f["sign_flip"])

    product_law_pass = kappa_spread <= KAPPA_CONSTANCY_TOL
    neutrality_pass = neutral_max_mod <= ZERO_MOD_TOL

    print(f"charged cells: {len(charged)}, neutral cells: {len(neutral)}")
    print(
        f"vertex-product law: kappa mean={float(np.mean(kappas)):.15g}"
        f" spread={kappa_spread:.3e}"
        f" -> {'PASS' if product_law_pass else 'FAIL'}"
    )
    print(f"delta spread={delta_spread:.3e}")
    print(
        f"neutrality (f=0 cells): max modulation={neutral_max_mod:.3e}"
        f" -> {'PASS' if neutrality_pass else 'FAIL'}"
    )
    print(f"sign flip across half-turn: {n_sign_flips}/{len(charged)} charged cells")

    csv_path = HERE / "cross_term_charge_pre_rows_v1.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    for f in charged:
        if f["f_A"] in (0.2, 0.6, 1.0) and f["f_B"] in (0.2, 0.6, 1.0):
            sub = [
                r for r in rows if r["f_A"] == f["f_A"] and r["f_B"] == f["f_B"]
            ]
            ax1.plot(
                [r["phi"] for r in sub],
                [r["dN_B_one_collision"] for r in sub],
                ".-",
                label=f"fA={f['f_A']}, fB={f['f_B']}",
            )
    ax1.axhline(0.0, color="0.5", linewidth=0.6)
    ax1.set_xlabel("phi (phase of B fermionic component)")
    ax1.set_ylabel("dN_B after one collision")
    ax1.set_title("Flow direction flips with relative phase (sign = charge)")
    ax1.legend(fontsize=6)
    ax1.grid(alpha=0.3)

    xs = [f["sqrt_fA_fB"] * abs(f["sin_2theta"]) for f in charged]
    ys = [f["modulation_M"] for f in charged]
    ax2.plot(xs, ys, "o", markersize=4)
    ax2.set_xlabel("sqrt(f_A f_B) * |sin 2theta|")
    ax2.set_ylabel("modulation amplitude M")
    ax2.set_title(
        f"Vertex-product law: M = kappa * sin2theta * sqrt(fA fB)"
        f" (spread {kappa_spread:.1e})"
    )
    ax2.grid(alpha=0.3)
    figure_names = []
    for ext in ("png", "svg"):
        path = HERE / f"cross_term_charge_pre_v1.{ext}"
        fig.savefig(path, dpi=160)
        figure_names.append(path.name)
    plt.close(fig)

    payload = {
        "experiment": "cross_term_charge_pre_v1",
        "prediction": "dN_B = offset + kappa*sin(2theta)*sqrt(fA fB)*cos(phi+delta); neutral cells unmodulated; sign flips across half-turn",
        "core_runner": {"path": TOY_RUNNER_PATH.name, "sha256": toy.sha256(TOY_RUNNER_PATH)},
        "verdicts": {
            "vertex_product_law": "PASS" if product_law_pass else "FAIL",
            "kappa_mean": float(np.mean(kappas)),
            "kappa_spread": kappa_spread,
            "delta_spread": delta_spread,
            "neutrality": "PASS" if neutrality_pass else "FAIL",
            "neutral_max_modulation": neutral_max_mod,
            "sign_flips": f"{n_sign_flips}/{len(charged)}",
        },
        "fits": fits,
        "figures": figure_names,
    }
    (HERE / "cross_term_charge_pre_result_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
