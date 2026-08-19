#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第3編（実験論文 DL0〜DL3）用の図 H1〜H7 を生成する（PNG/SVG）

データ源は保存済み走行系列（再計算しない）:
  dl0_series_v1.npz      … run_dl0_vacuum_v1.py（真空・v1・N=16・T=20000）
  dl23_series_v1.npz     … run_dl23_matter_v1.py（物質・δ=0.1・同上）
  dl3_perstep_series_v1.npz … probe_dl3_transport_perstep_v1.py（毎步輸送比）

H1 真空正単体: 固有値帯（×2M/c）の縮退と理論値 1
H2 真空の緩和と枠の非持続: top3 占有 0.343→0.200、重なり（ラグ120步）
H3 レジスタのジャンプと保存: Σx² の Re/Im（点火貫通・式値と一致）
H4 頂点パターンの出現: CV(A_v) 真空水準→0.54
H5 枠の誕生: g_min の開通（対数）と虚方向本数
H6 輸送比: 10步間隔（域外）vs 毎步（保証域の縁）——0.5 線
H7 尾部質量と楕円体偏差: tail・s_max/s_min の時系列

使い方: python3 make_paper3_figures_v1.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
plt.rcParams["font.family"] = ["Hiragino Sans", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False

D0 = np.load(HERE / "dl0_series_v1.npz")
D23 = np.load(HERE / "dl23_series_v1.npz")
DP = np.load(HERE / "dl3_perstep_series_v1.npz")
N, M = 16, 120


def _save(fig, stem):
    for ext in ("png", "svg"):
        fig.savefig(HERE / f"{stem}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  {stem}.png/.svg")


# ---------------------------------------------------------------- H1 真空正単体
tau0 = D0["tau"]
sc = 2 * M / D0["c"]
fig, ax = plt.subplots(figsize=(7, 4.2))
ax.fill_between(tau0, D0["lam_min"] * sc, D0["lam_max"] * sc,
                alpha=0.35, color="tab:blue", label="固有値帯（最小〜最大）")
ax.plot(tau0, D0["lam_med"] * sc, color="tab:blue", lw=1.2, label="中央値")
ax.axhline(1.0, color="k", ls="--", lw=1, label=r"理論値 $\lambda\cdot 2M/c=1$（正単体）")
ax.set_xlabel(r"$\tau$（步）")
ax.set_ylabel(r"$\lambda_i \cdot 2M/c$")
ax.set_title("H1 真空スペクトルの縮退（DL0・判定 J1/J2）")
ax.legend(loc="upper right", fontsize=9)
_save(fig, "fig_h1_vacuum_simplex")

# ------------------------------------------- H2 緩和と枠の非持続（2パネル）
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.0))
a1.plot(tau0, D0["top3"], color="tab:blue", lw=1.0)
a1.axhline(0.200, color="k", ls="--", lw=1, label=r"等方値 $3/(N-1)=0.200$")
a1.set_xlabel(r"$\tau$（步）")
a1.set_ylabel("上位3占有率")
a1.set_title("緩和 0.343→0.200（判定 J3/J7）")
a1.legend(fontsize=9)
ov = D0["overlap"]
ok = ~np.isnan(ov)
a2.plot(tau0[ok], ov[ok], color="tab:green", lw=0.8)
a2.axhline(0.200, color="k", ls="--", lw=1, label="乱数基準 0.200")
a2.axhline(0.5, color="tab:red", ls=":", lw=1.2, label="持続閾値 0.5")
a2.set_xlabel(r"$\tau$（步）")
a2.set_ylabel("上位3部分空間の重なり（ラグ120步）")
a2.set_title("枠の非持続（判定 J8/K5）")
a2.set_ylim(0, 1.0)
a2.legend(fontsize=9)
fig.suptitle("H2 真空に枠は立たない（DL0/DL1）", y=1.02)
_save(fig, "fig_h2_vacuum_relax_noframe")

# ------------------------------------------- H3 レジスタのジャンプと保存
tau2 = D23["tau"]
fig, ax = plt.subplots(figsize=(7, 4.2))
ax.plot(tau2, D23["bil_re"] * 1e3, color="tab:blue", lw=1.2,
        label=r"$\mathrm{Re}\,\Sigma x^2\ (\times10^{3})$")
ax.plot(tau2, D23["bil_im"] * 1e3, color="tab:orange", lw=1.2,
        label=r"$\mathrm{Im}\,\Sigma x^2\ (\times10^{3})$")
ax.axhline(7.168366670, color="tab:blue", ls="--", lw=0.8)
ax.axhline(2.991117345, color="tab:orange", ls="--", lw=0.8)
ax.set_xlabel(r"$\tau$（步）")
ax.set_ylabel(r"双線形レジスタ $\Sigma x_e^2\ (\times10^{3})$")
ax.set_title("H3 ジャンプ式の値が点火を貫通して保存（DL2・判定 L1/L2）\n"
             "破線＝初期データからの式値（実測との差 7×10⁻¹⁸・ドリフト 10⁻⁹）")
ax.legend(fontsize=9)
_save(fig, "fig_h3_register_jump_conservation")

# ------------------------------------------- H4 頂点パターン
fig, ax = plt.subplots(figsize=(7, 4.2))
ax.plot(tau2, D23["cvA"], color="tab:purple", lw=1.0)
ax.axhline(0.0539, color="k", ls="--", lw=1, label="真空水準 0.054（対照）")
ax.set_xlabel(r"$\tau$（步）")
ax.set_ylabel(r"$\mathrm{CV}(A_v)$")
ax.set_title("H4 頂点モーメントの空間パターン出現（DL2・判定 L6）")
ax.legend(fontsize=9)
_save(fig, "fig_h4_vertex_pattern")

# ------------------------------------------- H5 枠の誕生
fig, ax = plt.subplots(figsize=(7, 4.2))
ax.semilogy(tau2, np.maximum(D23["gmin"], 1e-9), color="tab:blue", lw=0.9,
            label=r"$g_{\min}$")
ax.set_xlabel(r"$\tau$（步）")
ax.set_ylabel(r"$g_{\min}$（対数）", color="tab:blue")
ax2 = ax.twinx()
ax2.plot(tau2, D23["n_neg"], color="tab:red", lw=0.7, alpha=0.6, label="虚方向本数")
ax2.set_ylabel("虚方向本数", color="tab:red")
ax.set_title("H5 枠の誕生: ギャップ開通と虚方向（DL3・判定 M1/M8）")
_save(fig, "fig_h5_frame_birth")

# ------------------------------------------- H6 輸送比
rho10 = D23["dB2"][1:] / np.maximum(D23["gmin"][:-1], 1e-300)
rho1 = DP["rho"]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.0))
bins = np.logspace(-2, 2, 60)
a1.hist(rho10, bins=bins, color="tab:gray", alpha=0.8)
a1.axvline(0.5, color="tab:red", ls=":", lw=1.5, label="保証境界 1/2")
a1.set_xscale("log")
a1.set_xlabel(r"$\rho=\|\Delta B\|_2/g_{\min}$")
a1.set_title("10步間隔: 中央値 3.74・全域外（定理の保証なし）")
a1.legend(fontsize=9)
a2.hist(rho1, bins=bins, color="tab:blue", alpha=0.8)
a2.axvline(0.5, color="tab:red", ls=":", lw=1.5, label="保証境界 1/2")
a2.axvline(float(np.median(rho1)), color="k", ls="--", lw=1,
           label=f"中央値 {np.median(rho1):.3f}")
a2.set_xscale("log")
a2.set_xlabel(r"$\rho=\|\Delta B\|_2/g_{\min}$")
a2.set_title("毎步: 中央値 0.477・域内 52%（保証域の縁）")
a2.legend(fontsize=9)
fig.suptitle("H6 枠輸送の輸送比（DL3・判定 M3。決着は直接検定＝図H9）", y=1.02)
_save(fig, "fig_h6_transport_ratio")

# ------------------------------------------- H7 尾部質量と楕円体偏差
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.0))
a1.plot(tau2, D23["tail"], color="tab:blue", lw=0.9)
a1.set_xlabel(r"$\tau$（步）")
a1.set_ylabel("尾部質量")
a1.set_title("尾部質量 0.62→0.49（判定 M6）")
a2.plot(tau2, D23["smax_smin"], color="tab:orange", lw=0.9)
a2.axhline(11.0, color="k", ls="--", lw=1, label="[V4] 実測水準 ≈11")
a2.set_xlabel(r"$\tau$（步）")
a2.set_ylabel(r"$s_{\max}/s_{\min}$")
a2.set_title("楕円体偏差（判定 M9・[V4] と同域）")
a2.legend(fontsize=9)
fig.suptitle("H7 rank-3 支配の定量と過渡の偏差（DL3）", y=1.02)
_save(fig, "fig_h7_tail_ellipsoid")

# ------------------------------------------- H8 点火の順序構造（毎步・§7 新事実2）
DI = np.load(HERE / "dl3_ignition_series_v1.npz")
fig, ax = plt.subplots(figsize=(8.5, 4.4))
ax.semilogy(DI["tau"], np.maximum(DI["gmin"], 1e-9), color="tab:blue", lw=0.7,
            label=r"$g_{\min}$（毎步）")
ax2 = ax.twinx()
ax2.plot(DI["tau"], DI["n_neg"], color="tab:red", lw=0.8, alpha=0.7,
         label=r"$\nu_-$（毎步）")
ax.axvline(9090, color="tab:red", ls="--", lw=1.2)
ax.axvline(10118, color="tab:blue", ls="--", lw=1.2)
ax.annotate(r"$\tau_-=9090$", (9090, 3e-3), color="tab:red", fontsize=10,
            xytext=(8350, 5e-3))
ax.annotate(r"$\tau_g=10118$", (10118, 1e-5), color="tab:blue", fontsize=10,
            xytext=(10180, 8e-6))
ax.set_xlabel(r"$\tau$（步）")
ax.set_ylabel(r"$g_{\min}$（対数）", color="tab:blue")
ax2.set_ylabel(r"虚方向本数 $\nu_-$", color="tab:red")
ax.set_title("H8 点火の順序構造: 虚方向が枠に約1000步先行する（§7 新事実2）")
_save(fig, "fig_h8_ignition_order")

# ------------------------------------------- H9 輸送の直接検定（§7 新事実1）
DT = np.load(HERE / "dl3_transport_compose_v1.npz")
th = DT["theta_disc_deg"]
de = DT["det_disc"]
fig, ax = plt.subplots(figsize=(8.5, 4.2))
w = np.arange(1, len(de) + 1)
ok = de > 0
ax.bar(w[ok], th[ok], color="tab:blue", label="一致窓（差角）")
ax.bar(w[~ok], [32] * int((~ok).sum()), color="tab:red",
       label="枝不一致（反射）イベント")
ax.axhline(30, color="k", ls=":", lw=1.2, label="事前登録しきい値 30°")
ax.set_xlabel("窓番号（各10步・物質相後期）")
ax.set_ylabel("合成輸送と直接輸送の差角（度）")
ax.set_title("H9 輸送の直接検定: 38/40窓で一致（中央値0.9°）・枝不一致2窓（§7 新事実1）")
ax.legend(fontsize=9)
_save(fig, "fig_h9_transport_direct")

print("完了")
