#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v2論文の新図生成（日英）: 保存済みJSONのみから描画（走行なし・決定論）
fig_v2a 偶奇非対称 / fig_v2b 線幅寿命 / fig_v2c 分裂読出し / fig_v2d 定常性処方"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

JA_FONT = None
for f in ("/System/Library/Fonts/ヒラギノ角ゴシック W3.ttc",
          "/System/Library/Fonts/Hiragino Sans GB.ttc"):
    if Path(f).exists():
        font_manager.fontManager.addfont(f)
        JA_FONT = font_manager.FontProperties(fname=f).get_name()
        break
plt.rcParams["axes.unicode_minus"] = False
HERE = Path(__file__).resolve().parent

T = {
"ja": dict(
 a_title="図V2a  偶奇結合の非対称性: 等振幅 f=0.500 は不動点でない",
 a_x="衝突数 j", a_y="フェルミオン分率 f",
 a_l1="f=0.5（対称なら不動点）", a_l2="実測（f0=0.500000から進化）",
 a_note="半並進射影 P±=(1±T)/2 の分離リーク: 偶側 {le:.0e} / 奇側 {lo:.0e}（厳密分離）",
 b_title="図V2b  線幅・寿命関係: τ_coh·σ_ω ≈ 一定（8構成・σ_ω 6倍レンジ）",
 b_x1="線幅 σ_ω（台上時計レートの空間std）", b_y1="寿命 τ_coh（C=1/2到達）",
 b_g="傾き−1（τ∝1/σ_ω）", b_y2="積 τ_coh·σ_ω", b_x2="構成",
 b_mean="平均 {m:.1f}（CV {cv:.0%}）",
 c_title="図V2c  分裂の読出し: 累積相対位相のπ超過＝2粒子化・t_split=π/Δω",
 c_x1="衝突数 j", c_y1="累積相対時計位相 |ΦA−ΦB|",
 c_pi="π（1粒子↔2粒子の自然基準）",
 c_x2="予言 π/Δω", c_y2="実測 t_split", c_diag="y=x",
 d_title="図V2d  定常粒子の処方: 全倍音・中心位相揃えでリンギング→0漸近",
 d_x="倍音幅 σ_k（積み上げ倍音数）", d_y="リンギング（瞬時場τ_tの時間std最大）",
 d_l1="全倍音・中心位相揃え（処方）", d_l2="偶数倍音のみ σ_k=16", d_l3="全倍音・位相ランダム σ_k=16",
 sfx="_ja"),
"en": dict(
 a_title="Fig. V2a  Asymmetry of the parity coupling: equal-amplitude f=0.500 is not a fixed point",
 a_x="collision j", a_y="fermionic fraction f",
 a_l1="f=0.5 (fixed point if symmetric)", a_l2="measured (evolved from f0=0.500000)",
 a_note="separation leaks of the half-translation projector P±=(1±T)/2: even {le:.0e} / odd {lo:.0e} (exact)",
 b_title="Fig. V2b  Linewidth-lifetime relation: τ_coh·σ_ω ≈ const (8 configs, 6× range of σ_ω)",
 b_x1="linewidth σ_ω (spatial std of clock rate on support)", b_y1="lifetime τ_coh (C=1/2 crossing)",
 b_g="slope −1 (τ∝1/σ_ω)", b_y2="product τ_coh·σ_ω", b_x2="configuration",
 b_mean="mean {m:.1f} (CV {cv:.0%})",
 c_title="Fig. V2c  Splitting readout: accumulated relative phase exceeding π = two particles; t_split=π/Δω",
 c_x1="collision j", c_y1="accumulated relative clock phase |ΦA−ΦB|",
 c_pi="π (natural 1-vs-2 particle criterion)",
 c_x2="prediction π/Δω", c_y2="measured t_split", c_diag="y=x",
 d_title="Fig. V2d  Stationary particles: ringing → 0 with full center-phased harmonics",
 d_x="harmonic width σ_k (number of stacked harmonics)", d_y="ringing (max time-std of τ_t)",
 d_l1="full harmonics, center-phased (prescription)", d_l2="even harmonics only, σ_k=16",
 d_l3="full harmonics, random phases, σ_k=16",
 sfx="_en"),
}

dA = json.loads((HERE / "pre_v2_parity_asymmetry_result_v1.json").read_text())
dB = json.loads((HERE / "pre_v2_linewidth_lifetime_result_v1.json").read_text())
dC = json.loads((HERE / "pre_v2_splitting_readout_result_v1.json").read_text())
dD = json.loads((HERE / "pre_v2_stationarity_result_v1.json").read_text())

for lang, t in T.items():
    plt.rcParams["font.family"] = JA_FONT if lang == "ja" else "DejaVu Sans"

    # V2a
    fs = np.array(dA["A"]["fs_decimated"]); js = np.arange(len(fs)) * dA["A"]["decimation"]
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.axhline(0.5, color="k", lw=1.0, ls="--", label=t["a_l1"])
    ax.plot(js, fs, lw=1.2, color="tab:red", label=t["a_l2"])
    ax.axhline(dA["A"]["f_star"], color="tab:blue", lw=1.0, ls=":",
               label=f"f* = {dA['A']['f_star']:.4f} ± {dA['A']['f_std']:.4f}")
    ax.set_xlabel(t["a_x"]); ax.set_ylabel(t["a_y"])
    ax.set_title(t["a_title"], fontsize=11)
    ax.text(0.02, 0.04, t["a_note"].format(le=dA["B"]["leak_even"], lo=dA["B"]["leak_odd"]),
            transform=ax.transAxes, fontsize=8.5,
            bbox=dict(boxstyle="round", fc="#f2f2f2", ec="#999999"))
    ax.legend(fontsize=9); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(HERE / f"fig_v2a_parity_asymmetry{t['sfx']}_v1.png", dpi=150); plt.close(fig)

    # V2b
    rows = [r for r in dB["rows"] if r["product"] is not None]
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.0))
    so = np.array([r["sigma_omega"] for r in rows]); tc = np.array([r["tau_coh"] for r in rows])
    axes[0].loglog(so, tc, "o", ms=8, color="tab:red")
    gx = np.array([so.min() * 0.8, so.max() * 1.2])
    axes[0].loglog(gx, dB["product_mean"] / gx, "k--", lw=1, label=t["b_g"])
    axes[0].set_xlabel(t["b_x1"]); axes[0].set_ylabel(t["b_y1"])
    axes[0].legend(fontsize=9); axes[0].grid(alpha=0.3, which="both")
    prods = np.array([r["product"] for r in rows])
    axes[1].bar(range(len(prods)), prods, color="tab:blue", alpha=0.8)
    axes[1].axhline(dB["product_mean"], color="tab:red", lw=1.5,
                    label=t["b_mean"].format(m=dB["product_mean"], cv=dB["product_cv"]))
    axes[1].set_xticks(range(len(rows)))
    axes[1].set_xticklabels([f"f{r['f_src']}\nA{r['amp']}\nσ{r['sig_k']:.0f}" for r in rows], fontsize=7)
    axes[1].set_ylabel(t["b_y2"]); axes[1].set_xlabel(t["b_x2"])
    axes[1].legend(fontsize=9); axes[1].grid(alpha=0.3, axis="y")
    fig.suptitle(t["b_title"])
    fig.tight_layout(); fig.savefig(HERE / f"fig_v2b_linewidth_lifetime{t['sfx']}_v1.png", dpi=150); plt.close(fig)

    # V2c
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.0))
    cmap = plt.get_cmap("viridis")
    for i, r in enumerate(dC["rows"]):
        ser = np.array(r["dphi_decimated"]); js = np.arange(len(ser)) * r["decimation"]
        axes[0].plot(js, ser, lw=1.3, color=cmap(i / max(len(dC["rows"]) - 1, 1)),
                     label=f"fA={r['fA']}, fB={r['fB']}")
    axes[0].axhline(np.pi, color="k", lw=1.2, ls="--", label=t["c_pi"])
    axes[0].set_xlabel(t["c_x1"]); axes[0].set_ylabel(t["c_y1"])
    axes[0].set_xlim(0, 300); axes[0].legend(fontsize=7); axes[0].grid(alpha=0.3)
    pred = [np.pi / r["delta_omega"] for r in dC["rows"] if r["t_split"]]
    meas = [r["t_split"] for r in dC["rows"] if r["t_split"]]
    axes[1].plot(pred, meas, "o", ms=9, color="tab:red")
    gmax = max(max(pred), max(meas)) * 1.2
    axes[1].plot([0, gmax], [0, gmax], "k--", lw=1, label=t["c_diag"])
    axes[1].set_xlabel(t["c_x2"]); axes[1].set_ylabel(t["c_y2"])
    axes[1].legend(fontsize=9); axes[1].grid(alpha=0.3)
    fig.suptitle(t["c_title"])
    fig.tight_layout(); fig.savefig(HERE / f"fig_v2c_splitting{t['sfx']}_v1.png", dpi=150); plt.close(fig)

    # V2d
    fig, ax = plt.subplots(figsize=(7, 4.2))
    sks, rings = [], []
    for name, c in dD["cases"].items():
        if name.startswith("full_phased"):
            sks.append(float(name.split("sk")[1])); rings.append(c["ringing"])
    order = np.argsort(sks)
    ax.semilogy(np.array(sks)[order], np.array(rings)[order], "o-", ms=8,
                color="tab:red", label=t["d_l1"])
    ax.axhline(dD["cases"]["even_only_sk16"]["ringing"], color="tab:gray", lw=1.4, ls="--",
               label=t["d_l2"])
    ax.axhline(dD["cases"]["full_random_sk16"]["ringing"], color="tab:blue", lw=1.2, ls=":",
               label=t["d_l3"])
    ax.set_xscale("log")
    ax.set_xlabel(t["d_x"]); ax.set_ylabel(t["d_y"])
    ax.set_title(t["d_title"], fontsize=11)
    ax.legend(fontsize=9); ax.grid(alpha=0.3, which="both")
    fig.tight_layout(); fig.savefig(HERE / f"fig_v2d_stationarity{t['sfx']}_v1.png", dpi=150); plt.close(fig)

print("figs V2a-V2d (ja/en) done")
