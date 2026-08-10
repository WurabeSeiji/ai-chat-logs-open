#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""論文A の主張図を作る v1（主張1つにつき1枚・計8枚）

入力（read-only）:
  - result_paperA_claims_v1.json（aggregate_paperA_claims_v1.py の出力）
  - 正本NPZ（時系列が要る図のみ）
  - result_divisor_class_register_order_v1.json（主張8）

出力: figA1_… 〜 figA8_….png

新規走行はしない。すべて既存の測定結果から描く。

使い方: python3 make_paperA_figures_v1.py
"""
from __future__ import annotations
import json
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = ["Hiragino Sans", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False

HERE = Path(__file__).resolve().parent
CLAIMS = json.loads((HERE / "result_paperA_claims_v1.json").read_text())
T_LONG = 42000
DELTAS = [1e-15, 1e-08, 1e-04, 1e-03, 1e-02,
          0.03162277660168379, 0.04357, 0.1]
MODES = ["neutral", "electron", "fermion_family", "boson_family", "mixed"]
LABEL = {"neutral": "中性型 (1,0)", "electron": "電子型 (1,3)",
         "fermion_family": "フェルミオン5セル型 k=1", "boson_family": "ボゾン3セル型 k=6",
         "mixed": "混合シード 8セル", "vacuum": "シードなし"}
COLOR = {"neutral": "tab:blue", "electron": "tab:orange",
         "fermion_family": "tab:green", "boson_family": "tab:red",
         "mixed": "tab:purple", "vacuum": "gray"}


def npz(mode, delta=None, suffix=""):
    tag = f"_rep-{suffix}" if suffix else ""
    if mode == "vacuum":
        return HERE / f"nsweep_vacuum_T{T_LONG}_N12_v2.npz"
    if not suffix and delta is not None and abs(delta - 0.01) < 1e-18 \
            and mode in ("mixed", "neutral"):
        return HERE / f"nsweep_{mode}_T{T_LONG}_N12_v2.npz"
    return HERE / f"nsweep_{mode}_T{T_LONG}_d{delta:g}{tag}_N12_v2.npz"


def ledger_time_axis(z, n_snap):
    """帳簿スナップショットの時刻。母体は t==1 と t%50==0 で保存する。
    `ledger_t` が保存されていない正本があるため、同じ規約で再構成する。"""
    if "ledger_t" in z.files:
        return z["ledger_t"]
    ax = np.arange(n_snap, dtype=float) * 50.0
    ax[0] = 1.0
    return ax


def late_cell_mean(z):
    """後半窓 [T/2, T) の 128 セル別時間平均。"""
    led = z["rec_m_ledger"]
    lt = ledger_time_axis(z, led.shape[0])
    sel = (lt >= T_LONG // 2) & (lt < T_LONG)
    return led[sel].mean(axis=0)


def save(fig, name, title):
    fig.suptitle(title, fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p = HERE / name
    fig.savefig(p, dpi=150)
    plt.close(fig)
    print(f"  → {name}")


# =============== 図A1: 三つの誕生は分離する ===============
def figA1():
    conds = [("vacuum", None, "シードなし\n（背景振動のみ）"),
             ("boson_family", 0.1, "偶数 k だけにシード\n（強度 0.1）"),
             ("fermion_family", 1e-08, "奇数 k にシード\n（強度 10⁻⁸）"),
             ("fermion_family", 1e-04, "奇数 k にシード\n（強度 10⁻⁴）")]
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
            a.text(ts[0] + 1, 0.42, f" 空間が生まれる\n {ts[0]+1} 回目",
                   fontsize=8, color="tab:blue")
        a.set_title(name, fontsize=10)
        a.set_ylim(-0.02, 1.05)
        a.set_xlim(1, len(t))
        a.grid(alpha=0.3, which="both")
        if j == 0:
            a.set_ylabel("空間 f₂\n（面外分率）")
        b = axes[1, j]
        odd = np.maximum(z["rec_m_odd_power"], 1e-40)
        b.loglog(t, odd, lw=1.0, color="tab:red", label="奇数 k のパワー")
        acq = z["m_acq"]
        idx = np.flatnonzero(acq)
        if len(idx):
            b.axvline(idx[0] + 1, color="tab:green", lw=1.6,
                      label=f"時計が生まれる（{idx[0]+1} 回目）")
        else:
            b.text(0.5, 0.62, "時計は生まれない", transform=b.transAxes,
                   ha="center", fontsize=10, color="tab:green",
                   bbox=dict(fc="honeydew", ec="tab:green"))
        if float(np.max(z["rec_m_odd_power"])) == 0.0:
            b.text(0.5, 0.36, "奇数 k のパワーは全区間で厳密に 0\n（下端は描画床 10⁻⁴⁰）",
                   transform=b.transAxes, ha="center", fontsize=9,
                   bbox=dict(fc="mistyrose", ec="gray"))
        b.set_ylim(1e-43, 1e2)
        b.set_xlim(1, len(t))
        b.grid(alpha=0.3, which="both")
        b.set_xlabel("更新回数（対数軸）")
        b.legend(fontsize=7, loc="upper left")
        if j == 0:
            b.set_ylabel("物質 = 奇数 k のパワー")
        sp = "○" if np.nanmax(np.nan_to_num(z["m_f2"], nan=0)) > 0.05 else "×"
        mt = "○" if float(np.max(z["rec_m_odd_power"])) > 1e-30 else "×"
        ck = "○" if acq.any() else "×"
        a.text(0.03, 0.86, f"空間 {sp}　物質 {mt}　時計 {ck}",
               transform=a.transAxes, fontsize=10,
               bbox=dict(fc="lightyellow", ec="gray", alpha=0.9))
    save(fig, "figA1_three_births_v1.png",
         "図A1　「粒子が生まれる」の中身は空間・物質・時計の三つに分かれ、独立に起きる")


# =============== 図A2: 置き場所が生成後の占有を決める ===============
def figA2():
    fig, axes = plt.subplots(1, 5, figsize=(16, 3.6))
    for a, mode in zip(axes, MODES):
        z = np.load(npz(mode, 0.1))
        cm = late_cell_mean(z)
        im = a.imshow(np.log10(np.maximum(cm, 1e-40)).T, aspect="auto",
                      origin="lower", cmap="viridis", vmin=-40, vmax=0)
        sup = int(np.count_nonzero(cm > 0.0))
        s1, s2 = cm.sum(), (cm ** 2).sum()
        a.set_title(f"{LABEL[mode]}\n支持 {sup}/128・N_eff {s1*s1/s2:.3f}",
                    fontsize=9)
        a.set_xlabel("巻き数 k")
        if mode == MODES[0]:
            a.set_ylabel("位相のずれ η")
        a.set_xticks(range(0, 16, 4))
        a.set_yticks(range(8))
    fig.colorbar(im, ax=axes, fraction=0.015, label="log₁₀ パワー（後半時間平均）")
    fig.suptitle("図A2　同じ強度・同じ誕生時刻でも、シードの置き場所で占有セルが変わる"
                 "（強度 0.1・42000 回更新後）", fontsize=13, y=1.16)
    p = HERE / "figA2_support_structure_v1.png"
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  → figA2_support_structure_v1.png")


# =============== 図A3: 偶奇選択則 ===============
def figA3():
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    a = axes[0]
    for mode in ("boson_family", "fermion_family"):
        z = np.load(npz(mode, 0.1))
        t = np.arange(1, len(z["rec_m_odd_power"]) + 1)
        a.semilogy(t, np.maximum(z["rec_m_odd_power"], 1e-45), lw=0.9,
                   color=COLOR[mode], label=LABEL[mode])
    a.set_ylim(1e-45, 1e2)
    a.set_xlabel("更新回数")
    a.set_ylabel("奇数 k のパワー")
    a.set_title("偶数 k だけのシードでは、奇数 k は厳密に 0 のまま", fontsize=10)
    a.legend(fontsize=8)
    a.grid(alpha=0.3)
    a.text(0.30, 0.12, "ボゾン3セル型は全 42000 回で厳密に 0\n"
                       "（下限に貼り付いた線は 10⁻⁴⁵ の描画床）",
           transform=a.transAxes, fontsize=8,
           bbox=dict(fc="mistyrose", ec="gray"))
    b = axes[1]
    per = CLAIMS["checks"]["claim3_odd_exact_zero"]["per_delta"]
    xs = [float(k) for k in per]
    ys = [per[k] for k in per]
    b.semilogx(xs, [max(v, 1e-45) for v in ys], "o", ms=9, color="tab:red")
    b.set_yscale("log")
    b.set_ylim(1e-45, 1e2)
    b.set_xlabel("シード強度 δ")
    b.set_ylabel("奇数 k のパワーの全区間最大")
    b.set_title("8 強度すべてで厳密に 0（丸め誤差ではなく選択則）", fontsize=10)
    b.grid(alpha=0.3, which="both")
    save(fig, "figA3_parity_selection_rule_v1.png",
         "図A3　偶奇選択則：出力の偶奇は入力3つの偶奇の和で決まるため、"
         "偶数側は奇数側を作れない")


# =============== 図A4: 強度の窓と山 ===============
def figA4():
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    a = axes[0]
    for i, mode in enumerate(MODES):
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
    xs, ys = [], []
    for d in DELTAS:
        r = CLAIMS["runs"].get(f"mixed|{d:g}")
        if r:
            xs.append(d)
            ys.append(r["r_nopump_peak"])
    b.semilogx(xs, ys, "o-", ms=8, color="tab:purple")
    imax = int(np.argmax(ys))
    b.plot(xs[imax], ys[imax], "*", ms=20, color="crimson", zorder=5)
    b.annotate(f"山: δ={xs[imax]:g}\nピーク {ys[imax]:.4f}",
               (xs[imax], ys[imax]), textcoords="offset points",
               xytext=(-90, -34), fontsize=9,
               arrowprops=dict(arrowstyle="->", color="crimson"))
    b.axhline(0.7, color="k", ls=":", lw=0.9)
    b.text(2e-14, 0.706, "0.7", fontsize=9)
    b.set_xlabel("シード強度 δ")
    b.set_ylabel("混合率のピーク高")
    b.set_title("強い側：単調でなく、中間の強度に山がある", fontsize=10)
    b.grid(alpha=0.3, which="both")
    save(fig, "figA4_amplitude_window_v1.png",
         "図A4　シード強度には有効な幅がある——下限（時計が立たない）と山（単調でない）")


# =============== 図A5: 助走の対数則 ===============
def figA5():
    c5 = CLAIMS["checks"]["claim5_fit"]
    pts, fa, fp = c5["points"], c5["fit_A_coh"], c5["fit_P_seed"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
    for a, key, fit, name in ((axes[0], "A_coh", fa, "コヒーレント振幅和  A = (セル数)×δ"),
                              (axes[1], "P_seed", fp, "パワー合計  P = (セル数)×δ²")):
        for mode in MODES:
            xs = [math.log(p[key]) for p in pts if p["mode"] == mode]
            ys = [p["tau_space"] for p in pts if p["mode"] == mode]
            if xs:
                a.plot(xs, ys, "o", ms=8, color=COLOR[mode], label=LABEL[mode])
        allx = np.array([math.log(p[key]) for p in pts])
        gx = np.linspace(allx.min(), allx.max(), 100)
        a.plot(gx, fit["a"] + fit["b"] * gx, "-", color="k", lw=1.4)
        a.set_xlabel(f"ln {name}")
        a.set_title(f"R² = {fit['r2']:.6f}　誤差 {fit['rmse']:.2f} 回",
                    fontsize=11)
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
    save(fig, "figA5_runup_log_law_v1.png",
         "図A5　助走期間を決めるのはパワーの合計ではなく、複素振幅の和である")


# =============== 図A6: 位相相殺 ===============
def figA6():
    c6 = CLAIMS["checks"]["claim6_phase_balanced"]["detail"]
    ks = sorted(c6, key=float)
    x = np.arange(len(ks))
    w = 0.26
    fig, ax = plt.subplots(figsize=(10, 4.6))
    series = [("control", "通常の混合シード（8セル）", "tab:purple"),
              ("phase_balanced", "位相を打ち消した混合シード（8セル・同パワー）", "tab:cyan"),
              ("F5", "フェルミオン5セル型（偶数成分なし）", "tab:green")]
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
    ax.set_ylim(0, 135)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, axis="y")
    ax.text(0.5, 0.80,
            "位相を打ち消すと、パワーも置き場所も変えていないのに\n"
            "「偶数成分を持たないシード」と 1 回単位で完全に一致する\n"
            "（4 強度すべて。測定前に事前登録した予言で 4/4 的中）",
            transform=ax.transAxes, ha="center", fontsize=10,
            bbox=dict(fc="honeydew", ec="gray"))
    save(fig, "figA6_phase_cancellation_v1.png",
         "図A6　効くのは置いたパワーの量ではなく、複素振幅の和である")


# =============== 図A7: 分解能の三段階 ===============
def figA7():
    c7 = CLAIMS["checks"]["claim7_resolution"]
    order = [("vacuum", 0), ("neutral", 1), ("mixed", 8)]
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.4))
    failed = c7["mixed"]["failed_N"] if "mixed" in c7 else []
    deg_counts = []
    for mode, ncell in order:
        if mode not in c7:
            continue
        rows = [r for r in c7[mode]["rows"] if r["built"]]
        Ns = [r["N"] for r in rows]
        ne = [r["n_eff_med"] for r in rows]
        al = [r["align_med"] for r in rows]
        lab = f"{LABEL[mode]}（{ncell} セル）"
        axes[0].plot(Ns, ne, "o-", ms=6, color=COLOR[mode], label=lab, alpha=0.9)
        axes[1].plot(Ns, al, "o-", ms=6, color=COLOR[mode], label=lab, alpha=0.9)
        deg = [n for n, v in zip(Ns, ne) if v > 1.8]
        deg_counts.append((mode, ncell, len(deg), len(Ns), deg))
    for a, ylab, ttl in (
            (axes[0], "n_eff（平面の縮退度）",
             "縮退する N はシードのセル数が多いほど減る"),
            (axes[1], "frame 整合度",
             "整合度は n_eff と逆に動く（同じ現象の別読み）")):
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
    axes[0].text(0.02, 0.02,
                 "黄 = N=3,4（全条件で縮退）\n赤 = 初期状態の構築に失敗（N=1,2,8,10）\n"
                 "点線 = 縮退の判定 n_eff = 1.8",
                 transform=axes[0].transAxes, fontsize=8, va="bottom",
                 bbox=dict(fc="lightyellow", ec="gray"))
    b = axes[2]
    xs = [d[1] for d in deg_counts]
    ys = [d[2] for d in deg_counts]
    cols = [COLOR[d[0]] for d in deg_counts]
    b.bar(range(len(xs)), ys, color=cols, edgecolor="k", lw=0.6, width=0.55)
    b.set_xticks(range(len(xs)))
    b.set_xticklabels([f"{LABEL[d[0]]}\nシード {d[1]} セル" for d in deg_counts],
                      fontsize=9)
    for i, d in enumerate(deg_counts):
        b.text(i, d[2] + 0.15, f"{d[2]} / {d[3]}\nN = {','.join(map(str, d[4]))}",
               ha="center", fontsize=9)
    b.set_ylim(0, max(ys) + 2.2)
    b.set_ylabel("平面が縮退している N の個数")
    b.set_title("シードは平面の縮退を解く", fontsize=10)
    b.grid(alpha=0.3, axis="y")
    save(fig, "figA7_resolution_regimes_v1.png",
         "図A7　分解能とシードの両方が平面の縮退を決める——"
         "N=3,4 は解けず、それ以外はシードのセル数が多いほど解ける")


# =============== 図A8: 約数類定理の位数追随 ===============
def figA8():
    d = json.loads((HERE / "result_divisor_class_register_order_v1.json").read_text())
    nes = [16, 8, 12, 32]
    cmap = {1: "tab:red", 2: "tab:blue", 3: "tab:purple",
            4: "tab:green", 8: "tab:orange"}
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.0), sharey=True)
    for a, ne in zip(axes, nes):
        rows = [r for r in d["runs"][str(ne)]["rows"]
                if r["settle"] == 4000 and r["m"] >= 1]
        for r in rows:
            a.plot(r["m"], r["mass2"], "o", ms=11,
                   color=cmap.get(r["gcd"], "gray"),
                   mec="k", mew=0.6, zorder=3)
        m1 = next((r for r in rows if r["m"] == 1), None)
        m3 = next((r for r in rows if r["m"] == 3), None)
        for r in (m1, m3):
            if r:
                a.annotate(f"m={r['m']}", (r["m"], r["mass2"]),
                           textcoords="offset points", xytext=(0, -22),
                           ha="center", fontsize=10, fontweight="bold")
        if m1 and m3:
            rel = abs(m1["mass2"] - m3["mass2"]) / abs(m3["mass2"])
            same = rel < 1e-5
            a.set_title(f"位数 ne = {ne}　"
                        f"gcd(1,{ne})={math.gcd(1,ne)} / gcd(3,{ne})={math.gcd(3,ne)}\n"
                        f"m=1 と m=3 は{'同類' if same else '別類'}"
                        f"（実測差 {rel:.1e}）",
                        fontsize=10, pad=8,
                        color="crimson" if not same else "black")
        a.set_xlabel("巻き数 m")
        a.margins(y=0.22)
        a.grid(alpha=0.3)
    axes[0].set_ylabel("質量²（分散補償・海中）")
    handles = [plt.Line2D([], [], marker="o", ls="", ms=10, mec="k",
                          color=c, label=f"gcd(m, ne) = {g}")
               for g, c in sorted(cmap.items())]
    axes[-1].legend(handles=handles, fontsize=8, loc="lower right")
    save(fig, "figA8_divisor_class_register_order_v1.png",
         "図A8　約数類定理はレジスタ位数に追随する——"
         "同じ m=1 と m=3 が位数 8・16・32 で14桁一致し、位数 12 でだけ 60% 違う")


if __name__ == "__main__":
    print("=== 論文A 主張図の作成 ===")
    for f in (figA1, figA2, figA3, figA4, figA5, figA6, figA7, figA8):
        f()
    print("完了：8 枚")
