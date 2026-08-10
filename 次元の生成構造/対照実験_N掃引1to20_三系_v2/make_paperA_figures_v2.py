#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""論文A の図を作る v2（全 15 枚）

v1（8枚）からの変更: 主張1つ1枚では実験の記録が本文に現れないため、
系の構成・実験行列の全体・通し走行の実体・帯の時間発展・全40条件の掃引・
N掃引の誕生マトリクス・再現対照を追加した。

  figA01 系の構成（128セル格子・ポンプ・シード・相棒セル）      §3
  figA02 実験行列の全体（5シード型 × 8強度）                  §3.9
  figA03 代表的な通し走行（4パネル）                          §3.6
  figA04 帯の時間発展と和則の実演                             §3.4
  figA05 三つの誕生                                          §4
  figA06 支持構造（128セル地図）                              §5
  figA07 偶奇選択則                                          §6
  figA08 強度の窓と山                                        §7
  figA09 空間の誕生時刻の全40条件                             §7,§8
  figA10 助走の対数則                                        §8.1
  figA11 位相相殺                                            §8.3
  figA12 分解能とシードによる縮退                             §9
  figA13 N掃引の誕生マトリクス                                §9
  figA14 約数類定理の位数追随                                 §10
  figA15 再現対照（40条件のτ照合・ne=16のビット一致）           §14

入力は既存の測定結果のみ（read-only）。新規走行はしない。

使い方: python3 make_paperA_figures_v2.py
"""
from __future__ import annotations
import json
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

plt.rcParams["font.family"] = ["Hiragino Sans", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False

HERE = Path(__file__).resolve().parent
CLAIMS = json.loads((HERE / "result_paperA_claims_v1.json").read_text())
DIVC = json.loads((HERE / "result_divisor_class_register_order_v1.json").read_text())
T_LONG = 42000
DELTAS = [1e-15, 1e-08, 1e-04, 1e-03, 1e-02,
          0.03162277660168379, 0.04357, 0.1]
MODES = ["neutral", "electron", "fermion_family", "boson_family", "mixed"]
LABEL = {"neutral": "中性型 (1,0)", "electron": "電子型 (1,3)",
         "fermion_family": "F5 フェルミオン5セル", "boson_family": "B3 ボゾン3セル",
         "mixed": "混合 8セル", "vacuum": "シードなし"}
COLOR = {"neutral": "tab:blue", "electron": "tab:orange",
         "fermion_family": "tab:green", "boson_family": "tab:red",
         "mixed": "tab:purple", "vacuum": "gray"}
SEED_CELLS = {"neutral": [(1, 0)], "electron": [(1, 3)],
              "fermion_family": [(1, e) for e in (0, 1, 3, 5, 6)],
              "boson_family": [(6, e) for e in (0, 3, 5)],
              "mixed": [(1, e) for e in (0, 1, 3, 5, 6)] + [(6, e) for e in (0, 3, 5)]}
SPECIES = {0: "ν型/γ・Z型", -1: "d型", -3: "e型/W⁻", 2: "u型", 3: "陽電子型/W⁺"}


def signed(e):
    return ((e + 4) % 8) - 4


def npz(mode, delta=None, suffix=""):
    tag = f"_rep-{suffix}" if suffix else ""
    if mode == "vacuum":
        return HERE / f"nsweep_vacuum_T{T_LONG}_N12_v2.npz"
    if not suffix and delta is not None and abs(delta - 0.01) < 1e-18 \
            and mode in ("mixed", "neutral"):
        return HERE / f"nsweep_{mode}_T{T_LONG}_N12_v2.npz"
    return HERE / f"nsweep_{mode}_T{T_LONG}_d{delta:g}{tag}_N12_v2.npz"


def ledger_axis(z, n):
    if "ledger_t" in z.files:
        return z["ledger_t"]
    a = np.arange(n, dtype=float) * 50.0
    a[0] = 1.0
    return a


def late_cell_mean(z):
    led = z["rec_m_ledger"]
    t = ledger_axis(z, led.shape[0])
    return led[(t >= T_LONG // 2) & (t < T_LONG)].mean(axis=0)


def save(fig, name, title, y=0.98, tight=(0, 0, 1, 0.94)):
    fig.suptitle(title, fontsize=12, y=y)
    fig.tight_layout(rect=tight)
    fig.savefig(HERE / name, dpi=150)
    plt.close(fig)
    print(f"  → {name}")


# ===================== 図A01 系の構成 =====================
def figA01():
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
    for a, mode in zip(axes, ("neutral", "electron", "boson_family")):
        grid = np.zeros((16, 8))
        for k in range(16):
            grid[k, :] = 0.12 if k % 2 == 1 else 0.0
        a.imshow(grid.T, origin="lower", aspect="auto", cmap="Blues",
                 vmin=0, vmax=1)
        a.add_patch(Rectangle((2 - .5, -.5), 1, 1, fill=True, fc="black",
                              alpha=0.85))
        a.text(2, 0, "P", color="w", ha="center", va="center",
               fontsize=11, fontweight="bold")
        for (k, e) in SEED_CELLS[mode]:
            a.add_patch(Rectangle((k - .5, e - .5), 1, 1, fill=True,
                                  fc="crimson", alpha=0.9))
            ks = (2 * 2 - k) % 16
            ms = (-signed(e)) % 8
            a.add_patch(Rectangle((ks - .5, ms - .5), 1, 1, fill=False,
                                  ec="darkgreen", lw=2.4))
            a.annotate("", xy=(ks, ms), xytext=(k, e),
                       arrowprops=dict(arrowstyle="->", color="darkgreen",
                                       lw=1.4, connectionstyle="arc3,rad=0.28"))
            if ks == 3:      # 粒子種のラベルが付いているのは k=3 帯だけ
                a.text(ks, ms + 0.62, SPECIES.get(signed(ms), ""),
                       ha="center", fontsize=8, color="darkgreen")
        a.set_xticks(range(0, 16, 2))
        a.set_yticks(range(8))
        a.set_yticklabels([f"{e}\n({signed(e):+d})" for e in range(8)], fontsize=7)
        a.set_xlabel("巻き数 k　（濃い列 = 奇数 = フェルミオン的）")
        a.set_title(f"{LABEL[mode]}　→　相棒帯 k*=(4−k) mod 16", fontsize=10)
        a.grid(which="minor", color="w", lw=0.5)
        a.set_xticks(np.arange(-.5, 16, 1), minor=True)
        a.set_yticks(np.arange(-.5, 8, 1), minor=True)
    axes[0].set_ylabel("位相のずれ η（符号付き m）")
    axes[0].text(0.02, 1.10, "■ P = ポンプ (2,0)　■ = シード　□ = 狙う相棒セル",
                 transform=axes[0].transAxes, fontsize=9,
                 bbox=dict(fc="lightyellow", ec="gray"))
    save(fig, "figA01_lattice_and_seeds_v1.png",
         "図A01　系の構成：各関係が持つ 16×8 = 128 セルと、ポンプ・シード・相棒セルの関係"
         "（和則 k* = 2·k_pump − k, m* = −m）", tight=(0, 0, 1, 0.90))


# ===================== 図A02 実験行列 =====================
def figA02():
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.2))
    M_ts = np.full((len(MODES), len(DELTAS)), np.nan)
    M_tt = np.full_like(M_ts, np.nan)
    M_sp = np.full_like(M_ts, np.nan)
    for i, m in enumerate(MODES):
        for j, d in enumerate(DELTAS):
            r = CLAIMS["runs"].get(f"{m}|{d:g}")
            if r is None:
                continue
            M_ts[i, j] = r["tau_space"] if r["tau_space"] else np.nan
            M_tt[i, j] = r["tau_time"] if r["tau_time"] else np.nan
            M_sp[i, j] = r["support_cells"]
    specs = [(M_ts, "空間が生まれる更新回数", "viridis_r", True),
             (M_tt, "時計が生まれる更新回数", "plasma_r", True),
             (M_sp, "支持セル数（128 中）", "cividis", False)]
    for a, (Mx, ttl, cm, logn) in zip(axes, specs):
        norm = matplotlib.colors.LogNorm() if logn else None
        im = a.imshow(Mx, aspect="auto", cmap=cm, norm=norm)
        for i in range(len(MODES)):
            for j in range(len(DELTAS)):
                v = Mx[i, j]
                a.text(j, i, "—" if np.isnan(v) else f"{int(v)}",
                       ha="center", va="center", fontsize=8,
                       color="crimson" if np.isnan(v) else "w")
        a.set_xticks(range(len(DELTAS)))
        a.set_xticklabels([f"{d:g}" for d in DELTAS], rotation=45, fontsize=8)
        a.set_yticks(range(len(MODES)))
        a.set_yticklabels([LABEL[m] for m in MODES], fontsize=9)
        a.set_xlabel("シード強度 δ")
        a.set_title(ttl, fontsize=10)
        fig.colorbar(im, ax=a, fraction=0.04)
    axes[1].text(0.5, -0.42, "赤い「—」＝ 生まれない",
                 transform=axes[1].transAxes, ha="center", fontsize=9,
                 color="crimson")
    save(fig, "figA02_experiment_matrix_v1.png",
         "図A02　実験行列の全体：5 シード型 × 8 強度 = 40 条件（N=12・42000 回更新・同一力学）",
         tight=(0, 0.04, 1, 0.93))


# ===================== 図A03 代表的な通し走行 =====================
def figA03():
    z = np.load(npz("mixed", 0.01))
    t = np.arange(1, T_LONG + 1)
    fig, ax = plt.subplots(4, 1, figsize=(11, 9), sharex=True)
    f2 = np.nan_to_num(z["m_f2"], nan=0.0)
    v2 = np.nan_to_num(z["v_f2"], nan=0.0)
    ax[0].semilogx(t, f2, lw=1.0, color="tab:blue", label="シードあり（混合 δ=0.01）")
    ax[0].semilogx(t, v2, lw=1.0, color="gray", ls="--", label="同走行の内部シードなし")
    ax[0].axhline(0.05, color="k", ls=":", lw=0.8)
    ts = int(np.flatnonzero(f2 > 0.05)[0]) + 1
    tv = int(np.flatnonzero(v2 > 0.05)[0]) + 1
    ax[0].axvline(ts, color="tab:blue", ls="--", lw=1.1)
    ax[0].axvline(tv, color="gray", ls="--", lw=1.1)
    ax[0].set_ylabel("空間 f₂")
    ax[0].legend(fontsize=8, loc="lower right")
    ax[0].set_title(f"空間が生まれる：{ts} 回目（シードなしなら {tv} 回目）", fontsize=10)
    ax[1].loglog(t, np.maximum(z["rec_m_odd_power"], 1e-45), lw=1.0,
                 color="tab:red", label="奇数 k（フェルミオン的）")
    ax[1].loglog(t, np.maximum(z["rec_m_even_power"], 1e-45), lw=1.0,
                 color="tab:cyan", label="偶数 k（ボゾン的）")
    ax[1].set_ylabel("帯パワー")
    ax[1].legend(fontsize=8, loc="lower right")
    ax[1].set_title("物質が生まれる：奇数 k のパワーが立つ", fontsize=10)
    acq = z["m_acq"]
    ax[2].semilogx(t, acq.astype(float), lw=1.0, color="tab:green")
    tt = int(np.flatnonzero(acq)[0]) + 1
    ax[2].set_ylabel("時計の取得可否")
    ax[2].set_ylim(-0.1, 1.15)
    ax[2].set_title(f"時計が生まれる：{tt} 回目以降、固有時間が取得可能", fontsize=10)
    ax[3].semilogx(t, z["m_n_eff"], lw=1.0, color="tab:purple", label="n_eff")
    ax[3].semilogx(t, z["m_align"], lw=1.0, color="tab:brown", label="frame 整合度")
    ax[3].axhline(1.8, color="k", ls=":", lw=0.8)
    ax[3].set_ylabel("平面の縮退度／整合度")
    ax[3].set_xlabel("更新回数（対数軸）")
    ax[3].legend(fontsize=8)
    ax[3].set_title("平面が一意に選ばれる（n_eff → 1）", fontsize=10)
    for a in ax:
        a.grid(alpha=0.3, which="both")
        a.set_xlim(1, T_LONG)
    save(fig, "figA03_full_run_example_v1.png",
         "図A03　代表的な通し走行（混合シード・δ=0.01・N=12・42000 回）——"
         "三つの誕生が一本の走行の中で順に起こる")


# ===================== 図A04 帯の時間発展と和則 =====================
def figA04():
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.6))
    for a, (mode, d) in zip(axes[:2], (("neutral", 0.01), ("boson_family", 0.1))):
        z = np.load(npz(mode, d))
        B = z["rec_m_bands"]
        t = np.arange(1, B.shape[0] + 1)
        k_seed = SEED_CELLS[mode][0][0]
        k_star = (4 - k_seed) % 16
        for k in sorted({k_seed, 2, k_star, 4, 5, 6, 14}):
            style = "-" if k in (k_seed, 2, k_star) else ":"
            lw = 1.8 if k in (k_seed, 2, k_star) else 0.9
            lab = (f"k={k}（シード）" if k == k_seed else
                   f"k={k}（ポンプ）" if k == 2 else
                   f"k={k}（相棒 k*）" if k == k_star else f"k={k}")
            a.loglog(t, np.maximum(B[:, k], 1e-45), style, lw=lw, label=lab)
        nz = int(sum(1 for k in range(16) if float(np.max(B[:, k])) > 0))
        a.set_ylim(1e-45, 1e2)
        a.set_xlim(1, B.shape[0])
        a.set_xlabel("更新回数（対数軸）")
        a.set_ylabel("帯パワー")
        a.set_title(f"{LABEL[mode]}・δ={d:g}：シード k={k_seed} → 相棒 k*={k_star}\n"
                    f"非ゼロになった帯：16 本中 {nz} 本", fontsize=10)
        a.legend(fontsize=7, ncol=2, loc="lower left")
        a.grid(alpha=0.3, which="both")
    axes[1].text(0.04, 0.55,
                 "偶数 k だけのシードでは\n相棒帯 k=14 すら立たない\n"
                 "非ゼロはポンプとシードの 2 本だけ",
                 transform=axes[1].transAxes, fontsize=9,
                 bbox=dict(fc="mistyrose", ec="gray"))
    b = axes[2]
    for mode, d in (("neutral", 0.01), ("fermion_family", 0.1),
                    ("mixed", 0.1), ("boson_family", 0.1)):
        z = np.load(npz(mode, d))
        r = z["rec_m_r_mean"]
        t = np.arange(1, len(r) + 1)
        b.loglog(t, np.maximum(r, 1e-30), lw=1.3, color=COLOR[mode],
                 label=f"{LABEL[mode]}（δ={d:g}）")
    b.set_ylim(1e-30, 1e0)
    b.set_xlim(1, T_LONG)
    b.set_xlabel("更新回数（対数軸）")
    b.set_ylabel("非線形部に掛かる読出し R")
    b.set_title("B3 では R が全区間で厳密に 0\n"
                "＝ 非線形頂点が一度も点火しない", fontsize=10)
    b.legend(fontsize=7, loc="lower right")
    b.grid(alpha=0.3, which="both")
    save(fig, "figA04_band_evolution_sumrule_v1.png",
         "図A04　帯の時間発展と和則、および偶数 k だけのシードで非線形部が点火しないこと")


# ===================== 図A05 三つの誕生 =====================
def figA05():
    conds = [("vacuum", None, "シードなし\n（ポンプのみ）"),
             ("boson_family", 0.1, "偶数 k だけにシード\n（δ=0.1）"),
             ("fermion_family", 1e-08, "奇数 k にシード\n（δ=10⁻⁸）"),
             ("fermion_family", 1e-04, "奇数 k にシード\n（δ=10⁻⁴）")]
    fig, axes = plt.subplots(2, 4, figsize=(15, 6.4), sharex=True)
    for j, (mode, d, name) in enumerate(conds):
        z = np.load(npz(mode, d))
        t = np.arange(1, len(z["m_f2"]) + 1)
        f2 = np.nan_to_num(z["m_f2"], nan=0.0)
        ts = np.flatnonzero(f2 > 0.05)
        a = axes[0, j]
        a.semilogx(t, f2, lw=1.0, color="tab:blue")
        a.axhline(0.05, color="k", ls=":", lw=0.8)
        if len(ts):
            a.axvline(ts[0] + 1, color="tab:blue", lw=1.2, ls="--")
            a.text(ts[0] + 1, 0.40, f" 空間\n {ts[0]+1} 回目", fontsize=8,
                   color="tab:blue")
        a.set_title(name, fontsize=10)
        a.set_ylim(-0.02, 1.05)
        a.set_xlim(1, len(t))
        a.grid(alpha=0.3, which="both")
        b = axes[1, j]
        b.loglog(t, np.maximum(z["rec_m_odd_power"], 1e-40), lw=1.0,
                 color="tab:red", label="奇数 k のパワー")
        idx = np.flatnonzero(z["m_acq"])
        if len(idx):
            b.axvline(idx[0] + 1, color="tab:green", lw=1.6,
                      label=f"時計（{idx[0]+1} 回目）")
        else:
            b.text(0.5, 0.62, "時計は生まれない", transform=b.transAxes,
                   ha="center", fontsize=10, color="tab:green",
                   bbox=dict(fc="honeydew", ec="tab:green"))
        if float(np.max(z["rec_m_odd_power"])) == 0.0:
            b.text(0.5, 0.34, "奇数 k は全区間で厳密に 0", transform=b.transAxes,
                   ha="center", fontsize=9,
                   bbox=dict(fc="mistyrose", ec="gray"))
        b.set_ylim(1e-43, 1e2)
        b.set_xlim(1, len(t))
        b.grid(alpha=0.3, which="both")
        b.set_xlabel("更新回数（対数軸）")
        b.legend(fontsize=7, loc="upper left")
        sp = "○" if len(ts) else "×"
        mt = "○" if float(np.max(z["rec_m_odd_power"])) > 1e-30 else "×"
        ck = "○" if z["m_acq"].any() else "×"
        a.text(0.03, 0.86, f"空間 {sp}　物質 {mt}　時計 {ck}", transform=a.transAxes,
               fontsize=10, bbox=dict(fc="lightyellow", ec="gray", alpha=0.9))
    axes[0, 0].set_ylabel("空間 f₂（面外分率）")
    axes[1, 0].set_ylabel("物質 = 奇数 k のパワー")
    save(fig, "figA05_three_births_v1.png",
         "図A05　空間・物質・時計は独立に起きる——4 通りの状況がすべて実現する")


# ===================== 図A06 支持構造 =====================
def figA06():
    fig, axes = plt.subplots(1, 5, figsize=(16, 3.8))
    for a, mode in zip(axes, MODES):
        z = np.load(npz(mode, 0.1))
        cm = late_cell_mean(z)
        im = a.imshow(np.log10(np.maximum(cm, 1e-40)).T, aspect="auto",
                      origin="lower", cmap="viridis", vmin=-40, vmax=0)
        sup = int(np.count_nonzero(cm > 0.0))
        s1, s2 = cm.sum(), (cm ** 2).sum()
        a.set_title(f"{LABEL[mode]}\n支持 {sup}/128・n_eff {s1*s1/s2:.3f}", fontsize=9)
        a.set_xlabel("巻き数 k")
        a.set_xticks(range(0, 16, 4))
        a.set_yticks(range(8))
        for (k, e) in SEED_CELLS[mode]:
            a.add_patch(Rectangle((k - .5, e - .5), 1, 1, fill=False,
                                  ec="crimson", lw=1.6))
    axes[0].set_ylabel("位相のずれ η")
    fig.colorbar(im, ax=axes, fraction=0.014, label="log₁₀ パワー（後半時間平均）")
    fig.suptitle("図A06　同じ強度・同じ誕生時刻でも、シードの置き場所で占有セルが変わる"
                 "（δ=0.1・42000 回更新後・赤枠＝シード位置）", fontsize=12, y=1.16)
    fig.savefig(HERE / "figA06_support_structure_v1.png", dpi=150,
                bbox_inches="tight")
    plt.close(fig)
    print("  → figA06_support_structure_v1.png")


# ===================== 図A07 偶奇選択則 =====================
def figA07():
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4))
    a = axes[0]
    for mode in ("boson_family", "fermion_family"):
        z = np.load(npz(mode, 0.1))
        t = np.arange(1, len(z["rec_m_odd_power"]) + 1)
        a.loglog(t, np.maximum(z["rec_m_odd_power"], 1e-45), lw=1.1,
                 color=COLOR[mode], label=LABEL[mode])
    a.set_ylim(1e-45, 1e2)
    a.set_xlabel("更新回数（対数軸）")
    a.set_ylabel("奇数 k のパワー")
    a.set_title("偶数 k だけのシードでは、奇数 k は厳密に 0 のまま", fontsize=10)
    a.legend(fontsize=8, loc="lower right")
    a.grid(alpha=0.3, which="both")
    a.text(0.05, 0.10, "B3 は全 42000 回で厳密に 0\n（下端の線は 10⁻⁴⁵ の描画床）",
           transform=a.transAxes, fontsize=8,
           bbox=dict(fc="mistyrose", ec="gray"))
    b = axes[1]
    per = CLAIMS["checks"]["claim3_odd_exact_zero"]["per_delta"]
    xs = [float(k) for k in per]
    b.semilogx(xs, [max(per[k], 1e-45) for k in per], "o", ms=10,
               color="tab:red", mec="k")
    for mode in ("fermion_family",):
        ys = [CLAIMS["runs"][f"{mode}|{d:g}"]["odd_power_max"] for d in DELTAS]
        b.semilogx(DELTAS, ys, "s", ms=8, color=COLOR[mode], mec="k",
                   label=LABEL[mode])
    b.set_yscale("log")
    b.set_ylim(1e-45, 1e2)
    b.set_xlabel("シード強度 δ")
    b.set_ylabel("奇数 k のパワーの全区間最大")
    b.set_title("8 強度すべてで厳密に 0（丸め誤差ではなく選択則）", fontsize=10)
    b.legend(handles=[plt.Line2D([], [], marker="o", ls="", color="tab:red",
                                 mec="k", label="B3（偶数 k のみ）"),
                      plt.Line2D([], [], marker="s", ls="", color="tab:green",
                                 mec="k", label="F5（奇数 k のみ）")], fontsize=8)
    b.grid(alpha=0.3, which="both")
    save(fig, "figA07_parity_selection_rule_v1.png",
         "図A07　偶奇選択則：出力の偶奇は入力 3 つの偶奇の和で決まるため、"
         "偶数側は奇数側を作れない")


# ===================== 図A08 強度の窓と山 =====================
def figA08():
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4))
    a = axes[0]
    for mode in MODES:
        xs, ys = [], []
        for d in DELTAS:
            r = CLAIMS["runs"].get(f"{mode}|{d:g}")
            if r is None:
                continue
            xs.append(d)
            ys.append(r["tau_time"] if r["tau_time"] is not None else np.nan)
        a.semilogx(xs, ys, "o-", ms=6, color=COLOR[mode], label=LABEL[mode],
                   alpha=0.85)
    a.axvspan(1e-16, 3e-8, color="lightgray", alpha=0.6)
    a.text(3e-13, 12, "時計が生まれない領域", fontsize=9, ha="center")
    a.set_xlabel("シード強度 δ")
    a.set_ylabel("時計が生まれる更新回数")
    a.set_title("下限：強度 10⁻⁸ 以下では時計が一度も生まれない", fontsize=10)
    a.legend(fontsize=7)
    a.grid(alpha=0.3, which="both")
    b = axes[1]
    xs = [d for d in DELTAS if f"mixed|{d:g}" in CLAIMS["runs"]]
    ys = [CLAIMS["runs"][f"mixed|{d:g}"]["r_nopump_peak"] for d in xs]
    b.semilogx(xs, ys, "o-", ms=8, color="tab:purple")
    imax = int(np.argmax(ys))
    b.plot(xs[imax], ys[imax], "*", ms=20, color="crimson", zorder=5)
    b.annotate(f"山: δ={xs[imax]:g}\nピーク {ys[imax]:.4f}",
               (xs[imax], ys[imax]), textcoords="offset points",
               xytext=(-95, -36), fontsize=9,
               arrowprops=dict(arrowstyle="->", color="crimson"))
    b.axhline(0.7, color="k", ls=":", lw=0.9)
    b.text(2e-14, 0.706, "0.7", fontsize=9)
    b.set_xlabel("シード強度 δ")
    b.set_ylabel("混合率のピーク高")
    b.set_title("強い側：単調でなく、中間の強度に山がある", fontsize=10)
    b.grid(alpha=0.3, which="both")
    save(fig, "figA08_amplitude_window_v1.png",
         "図A08　シード強度には有効な幅がある——下限（時計が立たない）と山（単調でない）")


# ===================== 図A09 空間の誕生時刻の全40条件 =====================
def figA09():
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
    a = axes[0]
    for mode in MODES:
        xs = [d for d in DELTAS if f"{mode}|{d:g}" in CLAIMS["runs"]]
        ys = [CLAIMS["runs"][f"{mode}|{d:g}"]["tau_space"] for d in xs]
        a.loglog(xs, ys, "o-", ms=6, color=COLOR[mode], label=LABEL[mode],
                 alpha=0.9)
    vac = CLAIMS["runs"]["mixed|0.1"]["tau_space_vacuum"]
    a.axhline(vac, color="gray", ls="--", lw=1.2)
    a.text(1.5e-15, vac * 1.08, f"シードなし {vac} 回", fontsize=9, color="gray")
    a.axvspan(1e-16, 3e-15, color="lightgray", alpha=0.5)
    a.text(1.3e-15, 300, "シードなしと\n区別がつかない", fontsize=8, ha="center")
    a.set_xlabel("シード強度 δ")
    a.set_ylabel("空間が生まれる更新回数")
    a.set_title("全 40 条件：強度を上げるほど早く立ち上がる", fontsize=10)
    a.legend(fontsize=7)
    a.grid(alpha=0.3, which="both")
    b = axes[1]
    for mode in MODES:
        xs = [d for d in DELTAS if f"{mode}|{d:g}" in CLAIMS["runs"]]
        ys = [CLAIMS["runs"][f"{mode}|{d:g}"]["tau_space"] for d in xs]
        n = len(SEED_CELLS[mode])
        b.semilogx([n * d for d in xs], ys, "o", ms=8, color=COLOR[mode],
                   label=f"{LABEL[mode]}（{n} セル）", alpha=0.9)
    b.set_xlabel("複素振幅の和  A = (セル数) × δ")
    b.set_ylabel("空間が生まれる更新回数")
    b.set_title("同じ量を A で並べ直すと 5 シード型が 1 本に重なる", fontsize=10)
    b.legend(fontsize=7)
    b.grid(alpha=0.3, which="both")
    save(fig, "figA09_tau_space_all_conditions_v1.png",
         "図A09　空間の誕生時刻（全 40 条件）——左：強度別／右：複素振幅の和で整理")


# ===================== 図A10 対数則 =====================
def figA10():
    c5 = CLAIMS["checks"]["claim5_fit"]
    pts, fa, fp = c5["points"], c5["fit_A_coh"], c5["fit_P_seed"]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), sharey=True)
    for a, key, fit, name in ((axes[0], "A_coh", fa, "複素振幅の和  A = (セル数)×δ"),
                              (axes[1], "P_seed", fp, "パワー合計  P = (セル数)×δ²")):
        for mode in MODES:
            xs = [math.log(p[key]) for p in pts if p["mode"] == mode]
            ys = [p["tau_space"] for p in pts if p["mode"] == mode]
            if xs:
                a.plot(xs, ys, "o", ms=9, color=COLOR[mode], mec="k", mew=0.5,
                       label=LABEL[mode])
        allx = np.array([math.log(p[key]) for p in pts])
        gx = np.linspace(allx.min(), allx.max(), 100)
        a.plot(gx, fit["a"] + fit["b"] * gx, "-", color="k", lw=1.4)
        a.set_xlabel(f"ln {name}")
        a.set_title(f"R² = {fit['r2']:.6f}　誤差 {fit['rmse']:.2f} 回", fontsize=11)
        a.grid(alpha=0.3)
    axes[0].set_ylabel("空間が生まれる更新回数")
    axes[0].legend(fontsize=8)
    axes[0].text(0.03, 0.06,
                 f"τ = {fa['a']:.3f} {fa['b']:+.3f}·ln A\n"
                 f"最大残差 {fa['max_abs_resid']:.2f} 回（{fa['n']} 点）",
                 transform=axes[0].transAxes, fontsize=9,
                 bbox=dict(fc="honeydew", ec="gray"))
    axes[1].text(0.03, 0.06, "同じ 15 点でも\n誤差が 7.7 倍に悪化する",
                 transform=axes[1].transAxes, fontsize=9,
                 bbox=dict(fc="mistyrose", ec="gray"))
    axes[0].text(0.55, 0.90, "中性型は電子型の下に完全に重なる",
                 transform=axes[0].transAxes, fontsize=8, color="gray")
    save(fig, "figA10_runup_log_law_v1.png",
         "図A10　助走期間を決めるのはパワーの合計ではなく、複素振幅の和である（弱域 15 点）")


# ===================== 図A11 位相相殺 =====================
def figA11():
    c6 = CLAIMS["checks"]["claim6_phase_balanced"]["detail"]
    ks = sorted(c6, key=float)
    x = np.arange(len(ks))
    w = 0.26
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.8),
                             gridspec_kw={"width_ratios": [1.35, 1]})
    ax = axes[0]
    series = [("control", "通常の混合シード（8 セル）", "tab:purple"),
              ("phase_balanced", "位相を打ち消した混合シード（同パワー・同住所）", "tab:cyan"),
              ("F5", "F5（偶数成分なし・5 セル）", "tab:green")]
    for i, (key, lab, col) in enumerate(series):
        vals = [c6[k][key]["got"] for k in ks]
        bars = ax.bar(x + (i - 1) * w, vals, w, label=lab, color=col,
                      edgecolor="k", lw=0.5)
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v + 1.5, str(v),
                    ha="center", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels([f"δ = {float(k):g}" for k in ks])
    ax.set_ylabel("空間が生まれる更新回数")
    ax.set_ylim(0, 138)
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.3, axis="y")
    ax.text(0.02, 0.72, "位相を打ち消すと、パワーも置き場所も変えていないのに\n"
                        "F5 と 1 回単位で完全に一致する（4 強度すべて）\n"
                        "測定前に事前登録した予言で 4/4 的中",
            transform=ax.transAxes, fontsize=9,
            bbox=dict(fc="honeydew", ec="gray"))
    b = axes[1]
    z0 = np.load(npz("mixed", 0.04357, "cohmixv1-ctl1"))
    z1 = np.load(npz("mixed", 0.04357, "pbmixv1-r1"))
    z2 = np.load(npz("fermion_family", 0.04357))
    t = np.arange(1, T_LONG + 1)
    for z, lab, col in ((z0, "通常の混合シード", "tab:purple"),
                        (z1, "位相を打ち消した混合シード", "tab:cyan"),
                        (z2, "F5", "tab:green")):
        b.semilogx(t, np.nan_to_num(z["m_f2"], nan=0.0), lw=1.2, color=col,
                   label=lab, alpha=0.85)
    b.axhline(0.05, color="k", ls=":", lw=0.8)
    b.set_xlim(1, 400)
    b.set_ylim(-0.02, 1.02)
    b.set_xlabel("更新回数（対数軸）")
    b.set_ylabel("空間 f₂")
    b.set_title("δ=0.04357 の立ち上がり（相殺後は F5 に重なる）", fontsize=10)
    b.legend(fontsize=8, loc="lower right")
    b.grid(alpha=0.3, which="both")
    save(fig, "figA11_phase_cancellation_v1.png",
         "図A11　効くのは置いたパワーの量ではなく、複素振幅の和である")


# ===================== 図A12 分解能とシード =====================
def figA12():
    c7 = CLAIMS["checks"]["claim7_resolution"]
    order = [("vacuum", 0), ("neutral", 1), ("mixed", 8)]
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.4))
    failed = c7["mixed"]["failed_N"] if "mixed" in c7 else []
    deg_counts = []
    for mode, ncell in order:
        rows = [r for r in c7[mode]["rows"] if r["built"]]
        Ns = [r["N"] for r in rows]
        ne = [r["n_eff_med"] for r in rows]
        al = [r["align_med"] for r in rows]
        lab = f"{LABEL[mode]}（{ncell} セル）"
        axes[0].plot(Ns, ne, "o-", ms=6, color=COLOR[mode], label=lab, alpha=0.9)
        axes[1].plot(Ns, al, "o-", ms=6, color=COLOR[mode], label=lab, alpha=0.9)
        deg_counts.append((mode, ncell, [n for n, v in zip(Ns, ne) if v > 1.8],
                           len(Ns)))
    for a, ylab, ttl in (
            (axes[0], "n_eff（平面の縮退度）", "縮退する N はシードのセル数が多いほど減る"),
            (axes[1], "frame 整合度", "整合度は n_eff と逆に動く（同じ現象の別読み）")):
        for nf in failed:
            a.axvspan(nf - 0.4, nf + 0.4, color="lightcoral", alpha=0.35)
        a.axvspan(2.5, 4.5, color="khaki", alpha=0.35)
        a.set_xlabel("分解能 N（ノード数）")
        a.set_ylabel(ylab)
        a.set_title(ttl, fontsize=10)
        a.set_xticks(range(1, 21, 2))
        a.legend(fontsize=8, loc="upper right")
        a.grid(alpha=0.3)
    axes[0].axhline(1.8, color="k", ls=":", lw=0.9)
    axes[0].text(0.02, 0.02, "黄 = N=3,4（全条件で縮退）\n"
                             "赤 = 構築に失敗（N=1,2,8,10）\n点線 = 判定 n_eff = 1.8",
                 transform=axes[0].transAxes, fontsize=8, va="bottom",
                 bbox=dict(fc="lightyellow", ec="gray"))
    b = axes[2]
    ys = [len(d[2]) for d in deg_counts]
    b.bar(range(len(ys)), ys, color=[COLOR[d[0]] for d in deg_counts],
          edgecolor="k", lw=0.6, width=0.55)
    b.set_xticks(range(len(ys)))
    b.set_xticklabels([f"{LABEL[d[0]]}\nシード {d[1]} セル" for d in deg_counts],
                      fontsize=9)
    for i, d in enumerate(deg_counts):
        b.text(i, len(d[2]) + 0.15,
               f"{len(d[2])} / {d[3]}\nN = {','.join(map(str, d[2]))}",
               ha="center", fontsize=9)
    b.set_ylim(0, max(ys) + 2.2)
    b.set_ylabel("平面が縮退している N の個数")
    b.set_title("シードは平面の縮退を解く", fontsize=10)
    b.grid(alpha=0.3, axis="y")
    save(fig, "figA12_resolution_regimes_v1.png",
         "図A12　分解能とシードの両方が平面の縮退を決める——"
         "N=3,4 は解けず、それ以外はシードのセル数が多いほど解ける")


# ===================== 図A13 N掃引の誕生マトリクス =====================
def figA13():
    c7 = CLAIMS["checks"]["claim7_resolution"]
    modes = ["vacuum", "neutral", "electron", "mixed"]
    Ns = list(range(1, 21))
    fig, axes = plt.subplots(1, 2, figsize=(15, 4.6),
                             gridspec_kw={"width_ratios": [1.5, 1]})
    a = axes[0]
    rowlab, Mx = [], []
    for m in modes:
        rows = {r["N"]: r for r in c7[m]["rows"]}
        for ev, key in (("空間", "tau_space"), ("物質", "matter_born"),
                        ("時計", "time_born")):
            line = []
            for N in Ns:
                r = rows.get(N)
                if r is None or not r["built"]:
                    line.append(np.nan)
                elif key == "tau_space":
                    line.append(1.0 if r["tau_space"] else 0.0)
                else:
                    line.append(1.0 if r[key] else 0.0)
            Mx.append(line)
            rowlab.append(f"{LABEL[m]}／{ev}")
    Mx = np.array(Mx)
    cmap = matplotlib.colors.ListedColormap(["#d9d9d9", "#2b8cbe"])
    a.imshow(np.nan_to_num(Mx, nan=0.5), aspect="auto", cmap=cmap, vmin=0, vmax=1)
    for i in range(Mx.shape[0]):
        for j in range(Mx.shape[1]):
            v = Mx[i, j]
            a.text(j, i, "—" if np.isnan(v) else ("○" if v > 0.5 else "×"),
                   ha="center", va="center", fontsize=9,
                   color="crimson" if np.isnan(v) else "w")
    a.set_xticks(range(len(Ns)))
    a.set_xticklabels(Ns, fontsize=8)
    a.set_yticks(range(len(rowlab)))
    a.set_yticklabels(rowlab, fontsize=8)
    a.set_xlabel("分解能 N")
    a.set_title("赤い「—」＝ 系が構成できない（N=1,2,8,10）", fontsize=10)
    for y in (2.5, 5.5, 8.5):
        a.axhline(y, color="k", lw=0.8)
    b = axes[1]
    for m in modes:
        rows = [r for r in c7[m]["rows"] if r["built"] and r["tau_space"]]
        b.semilogy([r["N"] for r in rows], [r["tau_space"] for r in rows],
                   "o-", ms=6, color=COLOR[m], label=LABEL[m], alpha=0.9)
    b.set_xlabel("分解能 N")
    b.set_ylabel("空間が生まれる更新回数")
    b.set_xticks(range(1, 21, 2))
    b.set_title("N≥12 で助走期間は大きく変わらない", fontsize=10)
    b.legend(fontsize=8)
    b.grid(alpha=0.3, which="both")
    save(fig, "figA13_nsweep_birth_matrix_v1.png",
         "図A13　N 掃引（N=1〜20 × 4 条件）：シードなしでは空間だけが生まれる")


# ===================== 図A14 約数類定理 =====================
def figA14():
    nes = [16, 8, 12, 32]
    cmap = {1: "tab:red", 2: "tab:blue", 3: "tab:purple",
            4: "tab:green", 8: "tab:orange"}
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.2), sharey=True)
    for a, ne in zip(axes, nes):
        rows = [r for r in DIVC["runs"][str(ne)]["rows"]
                if r["settle"] == 4000 and r["m"] >= 1]
        for r in rows:
            a.plot(r["m"], r["mass2"], "o", ms=11, color=cmap.get(r["gcd"], "gray"),
                   mec="k", mew=0.6, zorder=3)
        m1 = next((r for r in rows if r["m"] == 1), None)
        m3 = next((r for r in rows if r["m"] == 3), None)
        for r in (m1, m3):
            if r:
                a.annotate(f"m={r['m']}", (r["m"], r["mass2"]),
                           textcoords="offset points", xytext=(0, -24),
                           ha="center", fontsize=10, fontweight="bold")
        if m1 and m3:
            rel = abs(m1["mass2"] - m3["mass2"]) / abs(m3["mass2"])
            same = rel < 1e-5
            a.set_title(f"位数 ne = {ne}　"
                        f"gcd(1,{ne})={math.gcd(1,ne)} / gcd(3,{ne})={math.gcd(3,ne)}\n"
                        f"m=1 と m=3 は{'同類' if same else '別類'}（差 {rel:.1e}）",
                        fontsize=10, pad=8,
                        color="crimson" if not same else "black")
        a.set_xlabel("巻き数 m")
        a.margins(y=0.24)
        a.grid(alpha=0.3)
    axes[0].set_ylabel("質量²（分散補償・海中）")
    handles = [plt.Line2D([], [], marker="o", ls="", ms=10, mec="k", color=c,
                          label=f"gcd(m, ne) = {g}") for g, c in sorted(cmap.items())]
    axes[-1].legend(handles=handles, fontsize=8, loc="lower right")
    save(fig, "figA14_divisor_class_register_order_v1.png",
         "図A14　約数類定理はレジスタ位数に追随する——"
         "同じ m=1 と m=3 が位数 8・16・32 で 14 桁一致し、位数 12 でだけ 60% 違う")


# ===================== 図A15 再現対照 =====================
def figA15():
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
    a = axes[0]
    det = CLAIMS["checks"]["tau_vs_result_json"]["details"]
    xs = [d["tau_space"]["ref"] for d in det if d["tau_space"]["ref"]]
    ys = [d["tau_space"]["got"] for d in det if d["tau_space"]["ref"]]
    a.loglog(xs, ys, "o", ms=9, color="tab:blue", mec="k", mew=0.5,
             label=f"τ_space（{len(xs)} 条件）")
    xs2 = [d["tau_time"]["ref"] for d in det if d["tau_time"]["ref"]]
    ys2 = [d["tau_time"]["got"] for d in det if d["tau_time"]["ref"]]
    a.loglog(xs2, ys2, "s", ms=8, color="tab:green", mec="k", mew=0.5,
             label=f"τ_time（{len(xs2)} 条件）")
    lim = [1, 3000]
    a.plot(lim, lim, "-", color="k", lw=1.0, zorder=0)
    a.set_xlim(*lim)
    a.set_ylim(*lim)
    a.set_xlabel("既存の結果 JSON の値")
    a.set_ylabel("NPZ から独立に再算出した値")
    a.set_title(f"40 条件・不一致 0（全点が対角線上）", fontsize=10)
    a.legend(fontsize=8)
    a.grid(alpha=0.3, which="both")
    b = axes[1]
    ref = json.loads((HERE.parent / "波の周期表検討" /
                      "pre_v10b_longcoupling_result_v1.json").read_text())
    idx = {(r["settle"], r["m"]): r for r in ref["rows"]}
    gx, gy = [], []
    for r in DIVC["runs"]["16"]["rows"]:
        q = idx.get((r["settle"], r["m"]))
        if q is None:
            continue
        for f in ("mass2", "S", "retention"):
            if q.get(f) is not None and r.get(f) is not None:
                gx.append(q[f])
                gy.append(r[f])
    b.plot(gx, gy, "o", ms=9, color="tab:red", mec="k", mew=0.5)
    lim2 = [min(gx) * 0.9, max(gx) * 1.1]
    b.plot(lim2, lim2, "-", color="k", lw=1.0, zorder=0)
    b.set_xlim(*lim2)
    b.set_ylim(*lim2)
    b.set_xlabel("前論文の保存済み結果")
    b.set_ylabel("本稿での再走行")
    rep = DIVC["reproduction_ne16"]
    b.set_title(f"ne=16 の再現：{rep['n_compared']} 件・"
                f"最大相対差 {rep['worst_rel']:.1e}（ビット一致）", fontsize=10)
    b.grid(alpha=0.3)
    save(fig, "figA15_reproduction_controls_v1.png",
         "図A15　再現対照：主張の数値は正本データから独立に再算出して照合してある")


if __name__ == "__main__":
    print("=== 論文A の図（15 枚）===")
    for f in (figA01, figA02, figA03, figA04, figA05, figA06, figA07, figA08,
              figA09, figA10, figA11, figA12, figA13, figA14, figA15):
        f()
    print("完了：15 枚")
