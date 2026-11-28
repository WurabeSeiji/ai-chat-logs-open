#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
第 2 速報用の図 2 種（日本語版・英語版、計 4 枚）を生成する専用スクリプト（計算なし）。
入力:
  - results/tri_loop_n64_case1_layer4.npz
  - results/tri_loop_n64_case2_layer4.npz
  - results/tri_loop_n64_probe_timeseries.npz
出力:
  - results/fig_r2_four_stage_cascade_ja.png / _en.png
  - results/fig_r2_probe_branching_ja.png / _en.png
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def setup_font(lang):
    if lang == "ja":
        plt.rcParams["font.sans-serif"] = ["Hiragino Sans", "Hiragino Maru Gothic Pro", "YuGothic", "DejaVu Sans"]
    else:
        plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Helvetica", "Arial"]
    plt.rcParams["axes.unicode_minus"] = False


def plot_fig1(lang):
    setup_font(lang)
    c1 = np.load(os.path.join("results", "tri_loop_n64_case1_layer4.npz"))
    c2 = np.load(os.path.join("results", "tri_loop_n64_case2_layer4.npz"))

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.2), dpi=200)
    ja = (lang == "ja")

    # (a) Case 1: ||P|| & SVD singular values (τ=0..1200)
    ax = axes[0, 0]
    t1 = c1["tau"]
    ax.plot(t1, c1["norm_P"], color="#1f77b4", lw=2.2, label=r"$\|P(\tau)\|$")
    ax.plot(t1, c1["sv4"][:, 0], color="#d62728", lw=1.6, ls="--", label=r"$\sigma_1$ (rank 1$\to$2$\to$4)")
    ax.plot(t1, c1["sv4"][:, 1], color="#ff7f0e", lw=1.4, ls=":", label=r"$\sigma_2$")
    ax.plot(t1, c1["sv4"][:, 2], color="#2ca02c", lw=1.6, ls="-.", label=r"$\sigma_3, \sigma_4$")
    ax.axhline(1.0 / np.sqrt(2), color="gray", lw=0.8, ls=":")
    ax.axhline(0.5, color="gray", lw=0.8, ls=":")
    ax.set_xlim(0, 1200)
    ax.set_ylim(-0.03, 1.08)
    ax.set_xlabel(r"更新ステップ $\tau$" if ja else r"Update step $\tau$", fontsize=11)
    ax.set_ylabel(r"$\|P\|$ および SVD 特異値 $\sigma_k$" if ja else r"$\|P\|$ and SVD singular values $\sigma_k$", fontsize=11)
    ax.set_title(
        "(a) Case 1（外側 $O_8$ 実数種）：4 段階・倍々相転移と静止固定点"
        if ja else
        "(a) Case 1 (Outer $O_8$ Seed): 4-Stage Doubling Cascade",
        fontsize=11.0, fontweight="bold"
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="center right", fontsize=9, framealpha=0.9)
    # 注釈
    ax.annotate("L2→L3\n(×2)", xy=(52, 0.48), xytext=(140, 0.56),
                arrowprops=dict(arrowstyle="->", color="#444"), fontsize=8.5)
    ax.annotate("L3→L4 (τ=928)\n||P||≈1.0, σ=0.5×4\nd1=0" if not ja else "第3→第4相転移\n(τ=928, ||P||≈1.0)\nσ=0.5×4, d1=0",
                xy=(928, 0.995), xytext=(550, 0.80),
                arrowprops=dict(arrowstyle="->", color="#1f77b4"), fontsize=8.5,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#1f77b4", alpha=0.9))

    # (b) Case 2: ||P|| & SVD singular values (τ=0..7200)
    ax = axes[0, 1]
    t2 = c2["tau"]
    ax.plot(t2, c2["norm_P"], color="#1f77b4", lw=2.2, label=r"$\|P(\tau)\|$")
    ax.plot(t2, c2["sv4"][:, 0], color="#d62728", lw=1.6, ls="--", label=r"$\sigma_1, \sigma_2$")
    ax.plot(t2, c2["sv4"][:, 2], color="#2ca02c", lw=1.6, ls="-.", label=r"$\sigma_3, \sigma_4$")
    ax.axhline(1.0 / np.sqrt(2), color="gray", lw=0.8, ls=":")
    ax.axhline(0.5, color="gray", lw=0.8, ls=":")
    ax.set_xlim(0, 7200)
    ax.set_ylim(-0.03, 1.02)
    ax.set_xlabel(r"更新ステップ $\tau$" if ja else r"Update step $\tau$", fontsize=11)
    ax.set_ylabel(r"$\|P\|$ および SVD 特異値 $\sigma_k$" if ja else r"$\|P\|$ and SVD singular values $\sigma_k$", fontsize=11)
    ax.set_title(
        "(b) Case 2（大域丸め誤差種）：約 5,000 歩の第 3 層からの第 4 相転移"
        if ja else
        "(b) Case 2 (Global Seed): ~5,000-Step L3 → L4 Transition",
        fontsize=11.0, fontweight="bold"
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="center left", fontsize=9, framealpha=0.9)
    ax.annotate("第3→第4相転移\n(τ=6,541)\n||P||: 0.431→0.865\nσ=0.5×4, d1=0" if ja else
                "L3→L4 (τ=6,541)\n||P||: 0.431→0.865\nσ=0.5×4, d1=0",
                xy=(6541, 0.865), xytext=(3600, 0.78),
                arrowprops=dict(arrowstyle="->", color="#1f77b4"), fontsize=8.5,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#1f77b4", alpha=0.9))

    # (c) Case 2: 非ゼロ自由度 nz(τ) の再加熱的バースト（初期拡大 + 全体）
    ax = axes[1, 0]
    ax.plot(t2, c2["nz"], color="#9467bd", lw=2.0, label="非ゼロ成分数 nz (>1e-5)" if ja else "Nonzero count nz (>1e-5)")
    ax.set_xlim(0, 7200)
    ax.set_ylim(0, 3400)
    ax.set_xlabel(r"更新ステップ $\tau$" if ja else r"Update step $\tau$", fontsize=11)
    ax.set_ylabel("非ゼロ成分数 nz" if ja else "Nonzero components nz", fontsize=11)
    ax.set_title(
        "(c) Case 2 の非ゼロ自由度推移：過冷却→再加熱的バースト→定着"
        if ja else
        "(c) Case 2 Active nz: Supercooling → Burst → Freeze",
        fontsize=11.0, fontweight="bold"
    )
    ax.grid(True, alpha=0.3)
    ax.annotate("再加熱的バースト\n(512→3,072→1,560)" if ja else "Reheating-like burst\n(512→3,072→1,560)",
                xy=(121, 3072), xytext=(700, 2800),
                arrowprops=dict(arrowstyle="->", color="#9467bd"), fontsize=8.5)
    ax.annotate("第4相転移後の定着\n(nz=2,534, d1=0)" if ja else "Layer 4 freeze\n(nz=2,534, d1=0)",
                xy=(6600, 2534), xytext=(4100, 2650),
                arrowprops=dict(arrowstyle="->", color="#9467bd"), fontsize=8.5)
    # インセット: τ=0..160 のバースト拡大
    ax_in = ax.inset_axes([0.22, 0.14, 0.42, 0.42])
    m_in = t2 <= 160
    ax_in.plot(t2[m_in], c2["nz"][m_in], color="#9467bd", lw=1.8)
    ax_in.set_title(r"初期拡大 ($\tau \leq 160$)" if ja else r"Early zoom ($\tau \leq 160$)", fontsize=8.5)
    ax_in.tick_params(labelsize=7.5)
    ax_in.grid(True, alpha=0.3)

    # (d) Case 2: 外側 8 軸エネルギー分布（τ=0, 1500, 7000）
    ax = axes[1, 1]
    x = np.arange(8)
    w = 0.26
    ax.bar(x - w, c2["eA_0"] * 100, width=w, color="#aec7e8", edgecolor="#444", lw=0.6,
           label=r"$\tau=0$ (第1層 $O_8\otimes M_8$: 均等 12.5%)" if ja else r"$\tau=0$ (Layer 1: uniform 12.5%)")
    ax.bar(x, c2["eA_1500"] * 100, width=w, color="#ffbb78", edgecolor="#444", lw=0.6,
           label=r"$\tau=1,500$ (第3層: 3軸 79.3%)" if ja else r"$\tau=1,500$ (Layer 3: 3 axes 79.3%)")
    ax.bar(x + w, c2["eA_7000"] * 100, width=w, color="#d62728", edgecolor="#444", lw=0.6,
           label=r"$\tau=7,000$ (第4層: A=4,2,7 が 84.8%, A=0 は 0.04%)" if ja else r"$\tau=7,000$ (Layer 4: A=4,2,7 = 84.8%, A=0 = 0.04%)")
    ax.set_xticks(x)
    ax.set_xticklabels([f"A={i}" + (" (Re)" if i == 0 else "") for i in range(8)], fontsize=9.5)
    ax.set_ylabel("外側 8 軸のエネルギー占有率 (%)" if ja else "Energy share across outer 8 axes (%)", fontsize=11)
    ax.set_title(
        "(d) Case 2 外側 8 軸の自発的対称性の破れ（実数軸消失・虚数 3 軸 84.8%）"
        if ja else
        "(d) Case 2 Outer 8 Axes (A=0→0.04%; 3 Imag Axes = 84.8%)",
        fontsize=11.0, fontweight="bold"
    )
    ax.set_ylim(0, 44)
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(loc="upper left", fontsize=8.5, framealpha=0.92)

    fig.suptitle(
        r"$N=64\;(O_8 \otimes O_8)$ 過冷却 512 マス状態からの 4 段階・倍々エンタングルメント相転移（$\tau = 0$〜$10,000$）"
        if ja else
        r"Four-Stage Doubling Entanglement Phase Transitions from Supercooled 512-Cell State in $N=64\;(O_8 \otimes O_8)$",
        fontsize=13.5, fontweight="bold", y=0.985
    )
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join("results", f"fig_r2_four_stage_cascade_{lang}.png")
    plt.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Saved {out}")


def plot_fig2(lang):
    setup_font(lang)
    d = np.load(os.path.join("results", "tri_loop_n64_probe_timeseries.npz"))
    tau = d["tau"]
    ja = (lang == "ja")

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.0), dpi=200)

    # (a) プローブ本体エネルギー E_P(τ) の 4 分岐
    ax = axes[0, 0]
    ax.plot(tau, d["P3a_eP"], color="#1f77b4", lw=2.2, label="P3a: 大域複素回転波（完全共存 0.0099 不変）" if ja else "P3a: Global complex wave (exact coexistence)")
    ax.plot(tau, d["P3b_eP"], color="#d62728", lw=2.0, ls="--", label="P3b: ヌル・カイラル波（τ=10 で吸収・消滅）" if ja else "P3b: Null/chiral wave (absorbed by τ=10)")
    ax.plot(tau, d["P3c_eP"], color="#2ca02c", lw=2.2, label="P3c: スピン2(σ_x)波（10倍超へエネルギー吸収）" if ja else "P3c: Spin-2 (σ_x) wave (10x energy growth)")
    ax.plot(tau, d["P5_eP"], color="#ff7f0e", lw=2.0, ls="-.", label="P5: 局在四元数凝縮体（0.0025 で安定束縛）" if ja else "P5: Local H4 condensate (stable bound)")
    ax.plot(tau, d["P4_eP"], color="#9467bd", lw=1.8, ls=":", label="P4: 局在単一サイト凝縮体（0.00014 で安定束縛）" if ja else "P4: Local 1-site condensate (stable bound)")
    ax.set_yscale("log")
    ax.set_xlim(0, 50)
    ax.set_ylim(1e-5, 0.25)
    ax.set_xlabel(r"更新ステップ $\tau$" if ja else r"Update step $\tau$", fontsize=11)
    ax.set_ylabel(r"プローブ本体エネルギー $E_P(\tau)$ (対数)" if ja else r"Probe core energy $E_P(\tau)$ (log)", fontsize=11)
    ax.set_title(
        "(a) 512 マス時間結晶背景 $H$ と近傍プローブの 4 分岐"
        if ja else
        "(a) 4-Way Branching of Probes Coupled to $H$",
        fontsize=11.0, fontweight="bold"
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower left", fontsize=8.3, framealpha=0.92)

    # (b) 全体の ||P||(τ)（P3b による真空の 2 倍相転移）
    ax = axes[0, 1]
    ax.plot(tau, d["P3b_nrm"], color="#d62728", lw=2.4, label="P3b 共存（||P||: 0.0469 → 0.0938 厳密2倍相転移！）" if ja else "With P3b (||P||: 0.0469 → 0.0938 exact 2× transition!)")
    ax.plot(tau, d["P3a_nrm"], color="#1f77b4", lw=1.8, ls="--", label="P3a 共存（||P|| = 0.0469 不変）" if ja else "With P3a (||P|| = 0.0469 invariant)")
    ax.plot(tau, d["P3c_nrm"], color="#2ca02c", lw=1.6, ls=":", label="P3c 共存（||P|| ≈ 0.0468 維持）" if ja else "With P3c (||P|| ≈ 0.0468)")
    ax.plot(tau, d["P5_nrm"], color="#ff7f0e", lw=1.6, ls="-.", label="P5 共存（||P|| = 0.0469 維持）" if ja else "With P5 (||P|| = 0.0469)")
    ax.set_xlim(0, 50)
    ax.set_ylim(0.04, 0.102)
    ax.set_xlabel(r"更新ステップ $\tau$" if ja else r"Update step $\tau$", fontsize=11)
    ax.set_ylabel(r"系全体の出力ノルム $\|P(\tau)\|$" if ja else r"Total output norm $\|P(\tau)\|$", fontsize=11)
    ax.set_title(
        "(b) ヌル・カイラル波 (P3b) が引き起こす背景真空の 2 倍相転移"
        if ja else
        "(b) Exact 2× Norm Transition by Null/Chiral Probe (P3b)",
        fontsize=11.0, fontweight="bold"
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="center right", fontsize=8.5, framealpha=0.92)

    # (c) 局在凝縮体 P4・P5 の「本体 E_P」vs「周囲の衣 E_oth」
    ax = axes[1, 0]
    ax.plot(tau, d["P5_eP"] * 1000, color="#ff7f0e", lw=2.2, label=r"P5 本体 $E_P$ ($\to 2.50\times 10^{-3}$)" if ja else r"P5 core $E_P$ ($\to 2.50\times 10^{-3}$)")
    ax.plot(tau, d["P5_eO"] * 1000, color="#ff7f0e", lw=2.2, ls="--", label=r"P5 周囲の衣 $E_{\mathrm{oth}}$ ($\to 2.50\times 10^{-3}$, 厳密 1:1)" if ja else r"P5 dressed cloud $E_{\mathrm{oth}}$ ($\to 2.50\times 10^{-3}$, 1:1)")
    ax.plot(tau, d["P4_eP"] * 1000, color="#9467bd", lw=2.0, label=r"P4 本体 $E_P$ ($\to 0.14\times 10^{-3}$)" if ja else r"P4 core $E_P$ ($\to 0.14\times 10^{-3}$)")
    ax.plot(tau, d["P4_eO"] * 1000, color="#9467bd", lw=2.0, ls="--", label=r"P4 周囲の衣 $E_{\mathrm{oth}}$ ($\to 0.83\times 10^{-3}$, 本体の約6倍)" if ja else r"P4 dressed cloud $E_{\mathrm{oth}}$ ($\to 0.83\times 10^{-3}$, ~6× core)")
    ax.set_xlim(0, 50)
    ax.set_ylim(0, 5.2)
    ax.set_xlabel(r"更新ステップ $\tau$" if ja else r"Update step $\tau$", fontsize=11)
    ax.set_ylabel(r"エネルギー ($\times 10^{-3}$)" if ja else r"Energy ($\times 10^{-3}$)", fontsize=11)
    ax.set_title(
        "(c) 局在凝縮体が $O_8$ 格子にまとう「衣（ドレスド場）」の形成"
        if ja else
        "(c) Formation of Dressed Cloud ($E_{\\mathrm{oth}}$) Around Localized Condensates",
        fontsize=11.5, fontweight="bold"
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=8.5, framealpha=0.92)

    # (d) ステップ間距離 d1_P(τ): P4 の減速と P3c の周期 2 ロック
    ax = axes[1, 1]
    ax.plot(tau[1:], d["P4_d1P"][1:], color="#9467bd", lw=2.2, label="P4 局在凝縮体（τ=2〜10 で d1≈0.12 へ大幅減速）" if ja else "P4 local condensate (strong slowdown d1≈0.12 at τ=2–10)")
    ax.plot(tau[1:], d["P3c_d1P"][1:], color="#2ca02c", lw=2.0, ls="--", label="P3c スピン2波（d1=0.661 厳密周期2ロック）" if ja else "P3c spin-2 wave (locked at d1=0.661, period-2)")
    ax.plot(tau[1:], d["P3a_d1P"][1:], color="#1f77b4", lw=1.5, alpha=0.75, label="P3a 大域回転波（自由回転継続）" if ja else "P3a global wave (free rotation)")
    ax.axhline(0.904, color="#9467bd", lw=0.9, ls=":", alpha=0.7)
    ax.text(30, 0.92, "P4 単体走行時 (0.904)" if ja else "P4 solo (0.904)", color="#9467bd", fontsize=8.5)
    ax.set_xlim(1, 50)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel(r"更新ステップ $\tau$" if ja else r"Update step $\tau$", fontsize=11)
    ax.set_ylabel(r"プローブ部の1歩変化量 $d_1(\tau)$" if ja else r"Probe 1-step projective shift $d_1(\tau)$", fontsize=11)
    ax.set_title(
        "(d) 背景結合に伴うダイナミクスの変化（P4 の減速・P3c の周期ロック）"
        if ja else
        "(d) Dynamical Response: Slowdown of P4 & Period-2 Locking of P3c",
        fontsize=11.5, fontweight="bold"
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right", fontsize=8.3, framealpha=0.92)

    fig.suptitle(
        r"512 マス時間結晶背景 $H\;(= O_8 \otimes M_8)$ による近傍プローブのスピン・局在性選別と「衣」の形成"
        if ja else
        r"Spin & Localization Selection and Dressed-Cloud Formation by 512-Cell Time-Crystal Background $H\;(= O_8 \otimes M_8)$",
        fontsize=13.5, fontweight="bold", y=0.985
    )
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join("results", f"fig_r2_probe_branching_{lang}.png")
    plt.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Saved {out}")


if __name__ == "__main__":
    for lg in ("ja", "en"):
        plot_fig1(lg)
        plot_fig2(lg)
