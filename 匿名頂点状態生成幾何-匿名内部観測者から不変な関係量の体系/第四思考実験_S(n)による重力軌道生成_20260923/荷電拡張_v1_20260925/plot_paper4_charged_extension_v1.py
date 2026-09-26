#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plotting / reproducibility package for Paper-4 charged extension v1.

Default input is the first stable Coulomb-dominated charged experiment:
  q_m = 1, (n_A,n_B)=(+1,-1), |F_C/F_G|=100,
  orbit_dilation=10, Newtonian conservative angular mode,
  GW RR + EM RR enabled, 3 radial orbits, 4000 steps/orbit.

Outputs
-------
1. xy orbit plot.
2. Additive/log-side five-state plots: (chi,p,e,phi,t), individually and 5-panel.
3. Multiplicative/exponential five-state plots:
     U=exp(i chi), P=p, E=e, H=exp(i phi), Q=exp(t/tau0),
   represented as arg(U), P, E, arg(H), Q versus processing step,
   individually and 5-panel.
4. First-experiment overview: orbit, normalized p/e drift, GW waveform,
   EM waveform, radiated powers, and a text diagnostics panel.
5. Derived state CSV containing U,H,Q and tau0.

The plotting code does not alter the numerical experiment results.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager


def set_japanese_font() -> None:
    for name in [
        "Noto Sans CJK JP", "Noto Sans JP", "IPAexGothic", "IPAGothic",
        "Yu Gothic", "Hiragino Sans",
    ]:
        try:
            font_manager.findfont(name, fallback_to_default=False)
            plt.rcParams["font.family"] = name
            break
        except Exception:
            pass
    plt.rcParams["axes.unicode_minus"] = False


def save_fig(fig, output_dir: Path, stem: str) -> None:
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(output_dir / f"{stem}.{ext}", dpi=180, bbox_inches="tight")
    plt.close(fig)


def line_figure(k, y, title, ylabel, output_dir, stem, xlabel="処理ステップ k"):
    fig, ax = plt.subplots(figsize=(10, 5.8))
    ax.plot(k, y, linewidth=1.6)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    save_fig(fig, output_dir, stem)


def main() -> None:
    ap = argparse.ArgumentParser()
    here = Path(__file__).resolve().parent
    default_results = here / "results_rho100_d10_v2"
    ap.add_argument("--input-csv", default=str(default_results / "charged_v1_qm1_n1_-1_rho100.csv"))
    ap.add_argument("--summary-json", default=str(default_results / "charged_v1_qm1_n1_-1_rho100_summary.json"))
    ap.add_argument("--output-dir", default=str(here / "figures_rho100_d10_v2"))
    args = ap.parse_args()

    set_japanese_font()
    input_csv = Path(args.input_csv)
    summary_json = Path(args.summary_json)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_csv)
    with summary_json.open("r", encoding="utf-8") as f:
        summary = json.load(f)

    k = df["k"].to_numpy()
    chi = df["chi"].to_numpy()
    p = df["p"].to_numpy()
    e = df["e"].to_numpy()
    phi = df["phi"].to_numpy()
    t = df["t"].to_numpy()
    x = df["x"].to_numpy()
    y = df["y"].to_numpy()
    r = df["r"].to_numpy()

    U = np.exp(1j * chi)
    H = np.exp(1j * phi)

    # Natural charged-orbit clock scale: initial Kepler period for alpha.
    a0 = float(summary["initial_orbit"]["a0"])
    alpha = float(summary["derived_couplings"]["alpha"])
    tau0 = 2.0 * math.pi * math.sqrt(a0 ** 3 / alpha)
    Q = np.exp(t / tau0)
    argU = np.mod(np.angle(U), 2.0 * math.pi)
    argH = np.mod(np.angle(H), 2.0 * math.pi)

    # Save derived state table.
    derived = df.copy()
    derived["U_re"] = U.real
    derived["U_im"] = U.imag
    derived["argU_mod_2pi"] = argU
    derived["P"] = p
    derived["E"] = e
    derived["H_re"] = H.real
    derived["H_im"] = H.imag
    derived["argH_mod_2pi"] = argH
    derived["Q"] = Q
    derived["tau0"] = tau0
    derived.to_csv(output_dir / "charged_v1_derived_five_states.csv", index=False)

    # Figure 01: orbit.
    fig, ax = plt.subplots(figsize=(7.4, 7.0))
    ax.plot(x, y, linewidth=1.4)
    ax.plot(x[0], y[0], "o", label="開始")
    ax.plot(x[-1], y[-1], "s", label="終了")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("荷電拡張 v1：クーロン優勢二体軌道\n"
                 r"$m_A/m_B=1$, $(n_A,n_B)=(+1,-1)$, $|F_C/F_G|=100$")
    ax.grid(True, alpha=0.3)
    ax.legend()
    save_fig(fig, output_dir, "figure01_charged_orbit")

    # Individual additive/log-side states.
    additive = [
        (chi, "状態1：χ（動径位相）", "χ [rad]", "figure02a_log_state01_chi"),
        (p,   "状態2：p（半直弦）", "p", "figure02b_log_state02_p"),
        (e,   "状態3：e（離心率）", "e", "figure02c_log_state03_e"),
        (phi, "状態4：φ（方位角）", "φ [rad]", "figure02d_log_state04_phi"),
        (t,   "状態5：t（物理時間）", "t", "figure02e_log_state05_t"),
    ]
    for yy, title, ylabel, stem in additive:
        line_figure(k, yy, "加法 / log側表現 — " + title, ylabel, output_dir, stem)

    # Combined additive/log-side states.
    fig, axes = plt.subplots(5, 1, figsize=(11, 15), sharex=True)
    for ax, (yy, title, ylabel, _) in zip(axes, additive):
        ax.plot(k, yy, linewidth=1.35)
        ax.set_ylabel(ylabel)
        ax.set_title(title, loc="left")
        ax.grid(True, alpha=0.25)
    axes[-1].set_xlabel("処理ステップ k")
    fig.suptitle("荷電拡張 v1：5状態の加法 / log側表現", y=0.995, fontsize=16)
    save_fig(fig, output_dir, "figure02_log_five_states_panel")

    # Individual multiplicative/exponential states.
    exponential = [
        (argU, r"状態1：$U=e^{i\chi}$", r"arg(U) mod $2\pi$ [rad]", "figure03a_exp_state01_argU"),
        (p,    r"状態2：$P=p$", "P", "figure03b_exp_state02_P"),
        (e,    r"状態3：$E=e$", "E", "figure03c_exp_state03_E"),
        (argH, r"状態4：$H=e^{i\phi}$", r"arg(H) mod $2\pi$ [rad]", "figure03d_exp_state04_argH"),
        (Q,    r"状態5：$Q=e^{t/\tau_0}$", "Q", "figure03e_exp_state05_Q"),
    ]
    for yy, title, ylabel, stem in exponential:
        line_figure(k, yy, "指数 / 乗法表現 — " + title, ylabel, output_dir, stem)

    # Combined multiplicative/exponential states.
    fig, axes = plt.subplots(5, 1, figsize=(11, 15), sharex=True)
    for ax, (yy, title, ylabel, _) in zip(axes, exponential):
        ax.plot(k, yy, linewidth=1.35)
        ax.set_ylabel(ylabel)
        ax.set_title(title, loc="left")
        ax.grid(True, alpha=0.25)
    axes[-1].set_xlabel("処理ステップ k")
    fig.suptitle(
        "荷電拡張 v1：5状態の指数 / 乗法表現\n"
        + rf"$\tau_0=2\pi\sqrt{{a_0^3/\alpha}}={tau0:.6g}$",
        y=0.995, fontsize=16,
    )
    save_fig(fig, output_dir, "figure03_exponential_five_states_panel")

    # First-orbit selection for waveform panels.
    one_orbit = chi <= 2.0 * math.pi + 1e-12
    chi1 = chi[one_orbit]

    # Figure 04: first experiment overview.
    fig, axes = plt.subplots(2, 3, figsize=(17, 10))

    ax = axes[0, 0]
    ax.plot(x, y, linewidth=1.2)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title("軌道（3動径周期）")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True, alpha=0.25)

    ax = axes[0, 1]
    ax.plot(k, p / p[0], label="p/p0")
    ax.plot(k, e / e[0], label="e/e0")
    ax.set_title("軌道要素の緩慢変化")
    ax.set_xlabel("処理ステップ k")
    ax.set_ylabel("初期値に対する比")
    ax.grid(True, alpha=0.25)
    ax.legend()

    ax = axes[0, 2]
    ax.plot(chi1, df.loc[one_orbit, "gw_plus"], label="h+")
    ax.plot(chi1, df.loc[one_orbit, "gw_cross"], label="hx")
    ax.set_title("GW 読み出し（第1周）")
    ax.set_xlabel("χ [rad]")
    ax.set_ylabel("normalized strain readout")
    ax.grid(True, alpha=0.25)
    ax.legend()

    ax = axes[1, 0]
    ax.plot(chi1, df.loc[one_orbit, "em_x"], label="EM x")
    ax.plot(chi1, df.loc[one_orbit, "em_y"], label="EM y")
    ax.set_title("EM 双極子波 読み出し（第1周）")
    ax.set_xlabel("χ [rad]")
    ax.set_ylabel("normalized field readout")
    ax.grid(True, alpha=0.25)
    ax.legend()

    ax = axes[1, 1]
    positive_gw = np.maximum(df["gw_power"].to_numpy(), np.finfo(float).tiny)
    positive_em = np.maximum(df["em_power"].to_numpy(), np.finfo(float).tiny)
    ax.plot(k, positive_gw, label="GW power")
    ax.plot(k, positive_em, label="EM power")
    ax.set_yscale("log")
    ax.set_title("放射パワー")
    ax.set_xlabel("処理ステップ k")
    ax.set_ylabel("normalized power (log)")
    ax.grid(True, alpha=0.25)
    ax.legend()

    ax = axes[1, 2]
    ax.axis("off")
    d = summary["diagnostics"]
    ro = summary["readout"]
    inp = summary["inputs"]
    dc = summary["derived_couplings"]
    text = (
        "第一荷電実験（安定化条件）\n\n"
        f"mass ratio mA/mB = {inp['mass_ratio_mA_over_mB']}\n"
        f"charges (nA,nB) = ({inp['nA']:+d},{inp['nB']:+d})\n"
        f"Gamma = {inp['Gamma']:.6g}\n"
        f"|FC/FG| = {dc['rho_coulomb_to_gravity']:.6g}\n"
        f"alpha = {dc['alpha']:.6g}\n"
        f"orbit dilation = {inp['orbit_dilation']:.6g}\n\n"
        f"max v/c = {d['max_speed_over_c']:.6g}\n"
        f"EM/GW avg power = {d['first_orbit_em_to_gw_power_ratio']:.6g}\n"
        f"max |det S - 1| = {d['max_det_S_minus_1']:.3e}\n\n"
        f"nu fit = {ro['nu_fit']:.12g}\n"
        f"Gamma beta_q^2 fit = {ro['gamma_beta2_fit']:.12g}\n"
        f"Gamma nA nB local fit = {ro['gamma_nprod_local_fit']:.12g}\n"
        f"tau0 = {tau0:.12g}"
    )
    ax.text(0.02, 0.98, text, va="top", ha="left", fontsize=11)

    fig.suptitle("第四実験 荷電拡張 v1：第一数値実験の図化", fontsize=17)
    save_fig(fig, output_dir, "figure04_first_charged_experiment_overview")

    # Figure 05: normalized geometry comparison with uncharged baseline if available.
    baseline = here / "results_baseline_newtonian" / "charged_v1_qm1_n1_-1_rho0.csv"
    if baseline.exists():
        db = pd.read_csv(baseline)
        xb = db["x"].to_numpy()
        yb = db["y"].to_numpy()
        ab = 30.0
        fig, ax = plt.subplots(figsize=(7.4, 7.0))
        ax.plot(xb / ab, yb / ab, label="uncharged baseline, normalized")
        ax.plot(x / a0, y / a0, label="charged |FC/FG|=100, normalized")
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("x / initial semi-major axis")
        ax.set_ylabel("y / initial semi-major axis")
        ax.set_title("無荷電基準と荷電系の規格化軌道形状")
        ax.grid(True, alpha=0.25)
        ax.legend()
        save_fig(fig, output_dir, "figure05_baseline_vs_charged_normalized_orbit")

    print(f"Wrote plots to: {output_dir}")
    print(f"tau0 = {tau0:.15g}")


if __name__ == "__main__":
    main()
