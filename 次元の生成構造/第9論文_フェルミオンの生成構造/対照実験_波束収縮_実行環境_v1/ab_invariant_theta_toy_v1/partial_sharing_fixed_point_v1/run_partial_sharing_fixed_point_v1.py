#!/usr/bin/env python3
"""部分共有状態と固定点地図 v1（E1/E7）

E1 部分共有状態（二文法同時成立の設計要件の解消）:
    状態を「搬送波あり成分（対角読出しを駆動）＋搬送波なし共有成分
    （コヒーレント流路を開く）」の混合として構成する:
        a = sqrt(1-w) A_on(1) + sqrt(w) A_off(5)
        b = sqrt(1-w) B_on(7) + sqrt(w) e^{i phi} B_off(7)
    予言:
      (i)  正準 theta 読出しは非自明（対角セクター動作）——R はビン加法から厳密予言可能
      (ii) コヒーレント流の phi 変調 M = sin(2 theta) * w * |<A_off5|B_off7>| が非零
      両文法が単一状態対の上で同時に成立する。

E7 R4（交差スペクトル読出し）の固定点地図（見落とされた第二の動的選択チャネル）:
    R4 は自己駆動下で固定点に収束する（分類実験）。その吸引先 theta* が
    初期振幅にどう依存するかの地図を測り、特別な値（有理角・有限位数根）へ
    集中するかを判定する。
"""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
TOY_RUNNER_PATH = HERE.parent / "run_ab_invariant_theta_toy_v1.py"

spec = importlib.util.spec_from_file_location("toy_for_ps_fp_v1", TOY_RUNNER_PATH)
toy = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = toy
spec.loader.exec_module(toy)
base = toy.base
plt = base.plt

PRED_TOL = 1.0e-12


def unit_norm(v):
    return v / math.sqrt(float(np.vdot(v, v).real))


def mk(k, which, hair, sp):
    case = base.explicit_packet_case(mode=f"psfp_{which}_{k}_{int(hair)}", packet_a=(k,), packet_b=(k,))
    return unit_norm(base.make_case_state(sp, case, which, hair_enabled=hair))


def fit_cosine(phases, values):
    d = np.column_stack([np.ones_like(phases), np.cos(phases), np.sin(phases)])
    coef, *_ = np.linalg.lstsq(d, values, rcond=None)
    off, c, s = coef
    return float(off), float(math.hypot(c, s))


def masked_fraction(state, sp):
    shape = (sp.chi_grid_n, sp.eta_grid_n)
    f = np.fft.fft(state.reshape(shape), axis=0, norm="ortho")
    freqs = np.rint(np.fft.fftfreq(sp.chi_grid_n, d=1.0 / sp.chi_grid_n)).astype(int)
    fa = np.abs(freqs)
    mask = (fa >= 4) & (fa % 2 == 0)
    p = np.sum(np.abs(f) ** 2, axis=1)
    return float(np.sum(p[mask])) / float(np.sum(p))


def theta_relational(a, b, sp):
    shape = (sp.chi_grid_n, sp.eta_grid_n)
    af = np.fft.fft(a.reshape(shape), axis=0, norm="ortho")
    bf = np.fft.fft(b.reshape(shape), axis=0, norm="ortho")
    freqs = np.rint(np.fft.fftfreq(sp.chi_grid_n, d=1.0 / sp.chi_grid_n)).astype(int)
    f = np.abs(freqs)
    mask = (f >= 4) & (f % 2 == 0)
    cross = np.abs(np.sum(af * np.conj(bf), axis=1))
    pf = float(np.sum(cross[mask])); pb = float(np.sum(cross[~mask]))
    if pf + pb <= 0:
        return 0.0
    return math.atan2(math.sqrt(pf), math.sqrt(pb))


def main() -> None:
    params = base.Params(high_n=63, recursive_collision_count=1)
    sp = base.build_source_params(params)

    # ================= E1 部分共有 =================
    a_on1 = mk(1, "A", True, sp)
    b_on7 = mk(7, "B", True, sp)
    a_off5 = mk(5, "A", False, sp)
    b_off7 = mk(7, "B", False, sp)
    c_sh = complex(np.vdot(a_off5, b_off7))
    # 直交性検査（on/off 成分間の漏れ）
    leaks = [abs(complex(np.vdot(x, y))) for x, y in
             ((a_on1, a_off5), (b_on7, b_off7), (a_on1, b_off7), (a_off5, b_on7))]
    print(f"E1: shared coupling |c|={abs(c_sh):.6f}, component leaks max={max(leaks):.1e}")

    f_on1 = masked_fraction(a_on1, sp)
    f_on7 = masked_fraction(b_on7, sp)
    f_off5 = masked_fraction(a_off5, sp)
    f_off7 = masked_fraction(b_off7, sp)

    rows_e1 = []
    e1_r_err = 0.0
    e1_all_two_grammar = True
    phases = [2 * math.pi * k / 8 for k in range(8)]
    for w in (0.2, 0.4, 0.6):
        a0 = unit_norm(math.sqrt(1 - w) * a_on1 + math.sqrt(w) * a_off5)
        # R 予言（ビン加法）: f(a)= (1-w) f_on1 + w f_off5 等
        f_a = (1 - w) * f_on1 + w * f_off5
        f_b = (1 - w) * f_on7 + w * f_off7
        r_pred = (f_a + f_b) / 2.0
        flows = []
        r_meas = None
        for phi in phases:
            b0 = unit_norm(math.sqrt(1 - w) * b_on7 + math.sqrt(w) * np.exp(1j * phi) * b_off7)
            readout = toy.theta_from_ab(a0, b0, sp)
            if r_meas is None:
                r_meas = readout.reflection_rate
            a2, b2 = toy.rotate_ab(a0.copy(), b0.copy(), readout.theta)
            flows.append(float(np.vdot(b2, b2).real) - 1.0)
        off, m_mod = fit_cosine(np.asarray(phases), np.asarray(flows))
        theta_meas = math.asin(math.sqrt(r_meas))
        m_pred = abs(math.sin(2 * theta_meas)) * w * abs(c_sh)
        r_err = abs(r_meas - r_pred)
        m_err = abs(m_mod - m_pred)
        e1_r_err = max(e1_r_err, r_err, m_err)
        two_grammar = (0.01 < r_meas < 0.99) and (m_mod > 1e-3)
        e1_all_two_grammar &= two_grammar
        rows_e1.append(
            {"w": w, "R_measured": r_meas, "R_pred": r_pred, "R_err": r_err,
             "M_measured": m_mod, "M_pred": m_pred, "M_err": m_err,
             "both_grammars_active": two_grammar}
        )
        print(
            f"E1 w={w}: R={r_meas:.6f} (err {r_err:.1e})"
            f" M={m_mod:.6f} (err {m_err:.1e}) both_active={two_grammar}"
        )
    e1_pass = e1_r_err <= 1.0e-10 and e1_all_two_grammar
    print(f"E1: partial-sharing two-grammar coexistence -> {'PASS' if e1_pass else 'FAIL'}")

    # ================= E7 固定点地図 =================
    case = base.explicit_packet_case(
        mode="fp_map_b63", packet_a=(1,), packet_b=tuple(range(1, 64, 2))
    )
    a_t = unit_norm(base.make_case_state(sp, case, "A", hair_enabled=True))
    b_t = base.make_case_state(sp, case, "B", hair_enabled=True)

    amps = np.geomspace(0.05, 5.0, 41)
    rows_e7 = []
    for amp in amps:
        a = a_t.copy()
        b = float(amp) * unit_norm(b_t)
        th_prev = None
        for j in range(400):
            th = theta_relational(a, b, sp)
            a, b = toy.rotate_ab(a, b, th)
            th_prev = th
        th_star = th_prev
        rho = th_star / math.pi
        frac = Fraction(rho).limit_denominator(200)
        rows_e7.append(
            {"amplitude": float(amp), "theta_star_over_pi": rho,
             "nearest_rational": f"{frac.numerator}/{frac.denominator}",
             "dist": abs(rho - float(frac))}
        )
    # 固定点の集中判定: ヒストグラム的に見て特別値への集中があるか
    rhos = np.asarray([r["theta_star_over_pi"] for r in rows_e7])
    exact_hits = [r for r in rows_e7 if r["dist"] <= 1e-9 and r["nearest_rational"] not in ("0/1", "1/1")]
    print(f"E7: fixed points range [{rhos.min():.6f}, {rhos.max():.6f}],"
          f" exact rational hits: {len(exact_hits)}/{len(rows_e7)}")
    for r in exact_hits[:8]:
        print(f"   amp={r['amplitude']:.3f} theta*/pi={r['theta_star_over_pi']:.9f} = {r['nearest_rational']}")

    # ---- 保存 ----
    with (HERE / "partial_sharing_e1_rows_v1.csv").open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows_e1[0])); w.writeheader(); w.writerows(rows_e1)
    with (HERE / "fixed_point_map_e7_rows_v1.csv").open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows_e7[0])); w.writeheader(); w.writerows(rows_e7)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)
    ax1.errorbar([r["w"] for r in rows_e1], [r["M_measured"] for r in rows_e1],
                 fmt="o-", label="M measured")
    ax1.plot([r["w"] for r in rows_e1], [r["M_pred"] for r in rows_e1], "x--",
             label="M predicted")
    ax1b = ax1.twinx()
    ax1b.plot([r["w"] for r in rows_e1], [r["R_measured"] for r in rows_e1], "s-",
              color="tab:green", label="R (diagonal)")
    ax1.set_xlabel("sharing weight w"); ax1.set_ylabel("coherent modulation M")
    ax1b.set_ylabel("R (diagonal readout)")
    ax1.set_title("E1: both grammars active on one state pair")
    ax1.legend(loc="upper left", fontsize=8)
    ax2.semilogx(amps, rhos, ".-")
    for r in exact_hits:
        ax2.axhline(r["theta_star_over_pi"], color="tab:red", linewidth=0.5, alpha=0.4)
    ax2.set_xlabel("initial B amplitude"); ax2.set_ylabel("fixed point theta*/pi")
    ax2.set_title("E7: R4 fixed-point map")
    ax2.grid(alpha=0.3)
    for ext in ("png", "svg"):
        fig.savefig(HERE / f"partial_sharing_fixed_point_v1.{ext}", dpi=160)
    plt.close(fig)

    payload = {
        "experiment": "partial_sharing_fixed_point_v1",
        "core_runner": {"path": TOY_RUNNER_PATH.name, "sha256": toy.sha256(TOY_RUNNER_PATH)},
        "E1": {"rows": rows_e1, "max_err": e1_r_err,
               "verdict": "PASS" if e1_pass else "FAIL",
               "claim": "部分共有状態の上で対角読出し（重力文法）とコヒーレント変調（電荷文法）が同時成立"},
        "E7": {"n_points": len(rows_e7), "exact_rational_fixed_points": len(exact_hits),
               "range": [float(rhos.min()), float(rhos.max())]},
    }
    (HERE / "partial_sharing_fixed_point_result_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
