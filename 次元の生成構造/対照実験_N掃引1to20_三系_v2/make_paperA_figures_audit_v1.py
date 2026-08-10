#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""全480枚監査から出た発見の主張図を作る v1（figA16〜figA19・計4枚）

契機: `分析記録_全図面監査_v2.md`。走行図480枚を全数目視した結果、
本文の主張1・2・4・5・6・7 に訂正が必要と判明した。その根拠を図にする。

  figA16 飽和時刻の冪則        τ_sat ∝ P_F^(-1.09)（P_F = n_F·δ² ＝ フェルミオンパワー）
  figA17 二つの法則の対比      空間の誕生は複素振幅の和／物質の飽和はパワー
  figA18 選択性の窓            排他比は δ≈1e-3 で最大、飽和で6桁崩落
  figA19 六判定量のトレードオフ 上4行が増えるほど下2行が消える

判定の定義は母体 `run_tb_nsweep_1to20_v1.fig_matrix` から**逐語で引き継ぐ**:
  構成できた / 空間が生まれた(space_born) / 物質が生まれた(matter_born)
  / 時間が生まれた(time_born) / 3次元が確定(align_med > 0.8)
  / 凝縮体(cond_closure_med < 1e-5)
独自の閾値は導入しない。

飽和時刻 τ_sat は本器で定義する新しい量である:
  τ_sat = 物質分率 f_seed が 0.4 を最初に超えた更新回数（超えなければ None）
  0.4 は「飽和値 ≈ 0.5 の 8 割」であり、走行前ではなく本監査の結果から選んだ
  事後的な閾値である。感度は §出力 JSON に 0.3/0.4/0.45 の 3 通りを併記して示す。

入力はすべて既存の正本（read-only）。新規走行はしない。
出力: figA16〜A19 と `result_paperA_audit_figs_v1.json`

使い方: python3 make_paperA_figures_audit_v1.py
"""
from __future__ import annotations
import hashlib
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
T_LONG = 42000
DELTAS = [1e-15, 1e-08, 1e-04, 1e-03, 1e-02,
          0.03162277660168379, 0.04357, 0.1]
MODES = ["neutral", "electron", "fermion_family", "boson_family", "mixed"]
LABEL = {"neutral": "中性型 (1,0)", "electron": "電子型 (1,3)",
         "fermion_family": "F5 フェルミオン5セル", "boson_family": "B3 ボゾン3セル",
         "mixed": "混合 8セル"}
COLOR = {"neutral": "tab:blue", "electron": "tab:orange",
         "fermion_family": "tab:green", "boson_family": "tab:red",
         "mixed": "tab:purple"}
N_CELLS = {"neutral": 1, "electron": 1, "fermion_family": 5,
           "boson_family": 3, "mixed": 8}
N_FERM = {"neutral": 1, "electron": 1, "fermion_family": 5,
          "boson_family": 0, "mixed": 5}
SAT_LEVEL = 0.4                      # 飽和判定（本器で宣言・事後選択であることを明記）
SAT_LEVELS_SENS = (0.30, 0.40, 0.45)  # 感度確認用


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def npz_path(mode: str, delta: float, T: int = T_LONG, suffix: str = "") -> Path:
    tag = f"_rep-{suffix}" if suffix else ""
    if not suffix and T == T_LONG and abs(delta - 0.01) < 1e-18 \
            and mode in ("mixed", "neutral"):
        return HERE / f"nsweep_{mode}_T{T}_N12_v2.npz"
    return HERE / f"nsweep_{mode}_T{T}_d{delta:g}{tag}_N12_v2.npz"


def result_path(mode: str, delta: float) -> Path:
    if abs(delta - 0.01) < 1e-18 and mode in ("mixed", "neutral"):
        return HERE / f"result_nsweep_{mode}_T{T_LONG}_v2.json"
    return HERE / f"result_nsweep_{mode}_T{T_LONG}_d{delta:g}_v2.json"


def tau_sat(path: Path, level: float = SAT_LEVEL):
    """物質分率 f_seed が level を最初に超えた更新回数（母体と同じ 0基点+1 規約）。"""
    if not path.exists():
        return None
    fs = np.nan_to_num(np.load(path)["m_f_seed"], nan=0.0)
    idx = np.flatnonzero(fs > level)
    return int(idx[0]) + 1 if len(idx) else None


def collect_tau_sat() -> tuple[list, dict]:
    """T=42000 の 40 条件 ＋ 長時間走行（T=300000）から τ_sat を集める。"""
    pts, sens = [], {}
    for mode in MODES:
        for d in DELTAS:
            t = tau_sat(npz_path(mode, d))
            src = "T42000"
            if t is None:                      # 窓の外なら長時間走行を探す
                for suf in (f"s4-n12-{'f5' if mode=='fermion_family' else mode}"
                            f"-t300000-l1",):
                    p = npz_path(mode, d, 300000, suf)
                    if p.exists():
                        t = tau_sat(p)
                        src = "T300000"
            if t is None:
                continue
            pts.append({"mode": mode, "delta": d, "n_cells": N_CELLS[mode],
                        "n_ferm": N_FERM[mode], "P_F": N_FERM[mode] * d * d,
                        "A_coh": N_CELLS[mode] * d, "tau_sat": t, "source": src})
            for lv in SAT_LEVELS_SENS:
                t2 = tau_sat(npz_path(mode, d)) if src == "T42000" else None
                if src == "T300000":
                    p = npz_path(mode, d, 300000,
                                 f"s4-n12-{'f5' if mode=='fermion_family' else mode}-t300000-l1")
                    t2 = tau_sat(p, lv) if p.exists() else None
                else:
                    t2 = tau_sat(npz_path(mode, d), lv)
                sens.setdefault(f"{lv}", []).append(
                    {"mode": mode, "delta": d, "tau_sat": t2})
    return pts, sens


def linfit(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    b = (n * (x * y).sum() - x.sum() * y.sum()) / (n * (x * x).sum() - x.sum() ** 2)
    a = (y.sum() - b * x.sum()) / n
    pred = a + b * x
    ss_res = float(((y - pred) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    return {"a": float(a), "b": float(b), "n": n,
            "r2": float(1 - ss_res / ss_tot) if ss_tot > 0 else float("nan")}


def save(fig, name, title, tight=(0, 0, 1, 0.94)):
    fig.suptitle(title, fontsize=12)
    fig.tight_layout(rect=tight)
    fig.savefig(HERE / name, dpi=150)
    plt.close(fig)
    print(f"  → {name}")


# ===================== figA16 飽和時刻の冪則 =====================
def figA16(pts, out):
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    a = axes[0]
    for mode in MODES:
        q = [p for p in pts if p["mode"] == mode]
        if not q:
            continue
        a.loglog([p["delta"] for p in q], [p["tau_sat"] for p in q], "o-",
                 ms=8, color=COLOR[mode], label=LABEL[mode], alpha=0.9)
        for p in q:
            if p["source"] == "T300000":
                a.plot(p["delta"], p["tau_sat"], "*", ms=20, mfc="none",
                       mec="crimson", mew=1.8, zorder=5)
    a.axhline(T_LONG, color="k", ls="--", lw=1.2)
    a.text(1.2e-2, T_LONG * 1.15, "観測窓 T=42000", fontsize=9)
    a.set_xlabel("シード強度 δ")
    a.set_ylabel("物質が飽和する更新回数 τ_sat")
    a.set_title("★＝長時間走行（T=300000）で初めて観測できた点\n"
                "窓の上に出た点は「飽和しない」のではなく「窓の外で飽和する」", fontsize=10)
    a.legend(fontsize=8)
    a.grid(alpha=0.3, which="both")
    b = axes[1]
    x = [math.log(p["P_F"]) for p in pts if p["P_F"] > 0]
    y = [math.log(p["tau_sat"]) for p in pts if p["P_F"] > 0]
    fit = linfit(x, y)
    out["fit_PF"] = fit
    for mode in MODES:
        q = [p for p in pts if p["mode"] == mode and p["P_F"] > 0]
        if q:
            b.plot([math.log(p["P_F"]) for p in q],
                   [math.log(p["tau_sat"]) for p in q], "o", ms=10,
                   color=COLOR[mode], mec="k", mew=0.5, label=LABEL[mode])
    gx = np.linspace(min(x), max(x), 100)
    b.plot(gx, fit["a"] + fit["b"] * gx, "-", color="k", lw=1.6)
    b.set_xlabel("ln P_F　（P_F = シードのフェルミオンパワー = n_F·δ²）")
    b.set_ylabel("ln τ_sat")
    b.set_title(f"τ_sat ∝ P_F^({fit['b']:.3f})　R² = {fit['r2']:.4f}　"
                f"（{fit['n']} 点）", fontsize=11)
    b.grid(alpha=0.3)
    b.text(0.03, 0.06,
           "B3（ボゾン3セル）は P_F = 0 なので\nこの図に乗らない——実際どの強度でも永久に飽和しない",
           transform=b.transAxes, fontsize=9,
           bbox=dict(fc="mistyrose", ec="gray"))
    b.legend(fontsize=8, loc="upper right")
    save(fig, "figA16_saturation_power_law_v1.png",
         "図A16　物質が飽和するまでの時間は、シードのフェルミオンパワーの逆数に比例する"
         "——強度の「閾値」は存在しない")


# ===================== figA17 二つの法則 =====================
def figA17(pts, out):
    claims = json.loads((HERE / "result_paperA_claims_v1.json").read_text())
    c5 = claims["checks"]["claim5_fit"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    a = axes[0]
    for mode in MODES:
        q = [p for p in c5["points"] if p["mode"] == mode]
        if q:
            a.plot([math.log(p["A_coh"]) for p in q],
                   [p["tau_space"] for p in q], "o", ms=10, color=COLOR[mode],
                   mec="k", mew=0.5, label=LABEL[mode])
    fa = c5["fit_A_coh"]
    gx = np.linspace(min(math.log(p["A_coh"]) for p in c5["points"]),
                     max(math.log(p["A_coh"]) for p in c5["points"]), 100)
    a.plot(gx, fa["a"] + fa["b"] * gx, "-", color="k", lw=1.6)
    a.set_xlabel("ln A　（A = 複素振幅の和 = n·δ）")
    a.set_ylabel("空間が生まれる更新回数 τ_space")
    a.set_title(f"【線形部が決める】τ_space = {fa['a']:.3f} {fa['b']:+.3f}·ln A\n"
                f"R² = {fa['r2']:.6f}", fontsize=11, color="tab:blue")
    a.legend(fontsize=8)
    a.grid(alpha=0.3)
    a.text(0.03, 0.06, "線形部の回転生成子は\n全セルの複素和の偏角だけで決まる",
           transform=a.transAxes, fontsize=9,
           bbox=dict(fc="aliceblue", ec="tab:blue"))
    b = axes[1]
    fit = out["fit_PF"]
    for mode in MODES:
        q = [p for p in pts if p["mode"] == mode and p["P_F"] > 0]
        if q:
            b.plot([math.log(p["P_F"]) for p in q],
                   [math.log(p["tau_sat"]) for p in q], "s", ms=10,
                   color=COLOR[mode], mec="k", mew=0.5, label=LABEL[mode])
    xs = [math.log(p["P_F"]) for p in pts if p["P_F"] > 0]
    gx = np.linspace(min(xs), max(xs), 100)
    b.plot(gx, fit["a"] + fit["b"] * gx, "-", color="k", lw=1.6)
    b.set_xlabel("ln P_F　（P_F = パワー = n_F·δ²）")
    b.set_ylabel("ln τ_sat（物質が飽和する更新回数）")
    b.set_title(f"【非線形部が決める】τ_sat ∝ P_F^({fit['b']:.3f})\n"
                f"R² = {fit['r2']:.4f}", fontsize=11, color="tab:red")
    b.legend(fontsize=8)
    b.grid(alpha=0.3)
    b.text(0.03, 0.06, "非線形部の3次式は\n読出し R（パワー比）に比例する",
           transform=b.transAxes, fontsize=9,
           bbox=dict(fc="mistyrose", ec="tab:red"))
    save(fig, "figA17_two_laws_v1.png",
         "図A17　同じ系に性質の違う二つの法則がある"
         "——見る事象が線形部と非線形部のどちらに属するかで、効く量が変わる")


# ===================== figA18 選択性の窓 =====================
def figA18(out):
    fig, ax = plt.subplots(figsize=(10.5, 5.0))
    rows = []
    for mode in MODES:
        xs, ys = [], []
        for d in DELTAS:
            p = result_path(mode, d)
            if not p.exists():
                continue
            r = json.loads(p.read_text())["N"]["12"]
            v = r.get("excl_ratio_med")
            if v is None or not np.isfinite(v):
                continue
            xs.append(d)
            ys.append(max(v, 1e-16))
            rows.append({"mode": mode, "delta": d, "excl_ratio_med": v})
        if xs:
            ax.loglog(xs, ys, "o-", ms=8, color=COLOR[mode], label=LABEL[mode],
                      alpha=0.9)
    out["exclusion"] = rows
    peak = max((r for r in rows), key=lambda r: r["excl_ratio_med"])
    ax.axvspan(3e-4, 2e-2, color="honeydew", alpha=0.9, zorder=0)
    ax.plot(peak["delta"], peak["excl_ratio_med"], "*", ms=26, color="crimson",
            zorder=6)
    ax.annotate(f"選択性の頂点　δ={peak['delta']:g}\n排他比 {peak['excl_ratio_med']:.2e}"
                f"（1 万倍）",
                (peak["delta"], peak["excl_ratio_med"]),
                textcoords="offset points", xytext=(-235, -18), fontsize=10,
                color="crimson",
                arrowprops=dict(arrowstyle="->", color="crimson"))
    ax.text(2.4e-3, 1e-13, "選択性の窓", fontsize=12, ha="center",
            color="darkgreen", fontweight="bold")
    ax.axhline(1.0, color="k", ls=":", lw=0.9)
    ax.text(3e-15, 2.0, "排他比 = 1（狙ったセルと非狙いセルが同じ）", fontsize=8)
    ax.set_xlabel("シード強度 δ")
    ax.set_ylabel("排他比（狙った相棒セル ÷ 非狙いセル・後半中央値）")
    ax.set_ylim(1e-17, 1e6)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8, loc="center left")
    ax.annotate("飽和すると 128 セル全体に拡散し\n選択性が 6 桁崩落する",
                (0.1, 7.67e-2), textcoords="offset points", xytext=(-60, 95),
                fontsize=9, ha="center",
                arrowprops=dict(arrowstyle="->", color="gray"),
                bbox=dict(fc="mistyrose", ec="gray"))
    ax.text(1e-8, 1e-15, "B3（ボゾン3セル）は全条件で 0\n"
                         "＝ 狙った相棒セルが立たない", fontsize=8, color="tab:red")
    save(fig, "figA18_selectivity_window_v1.png",
         "図A18　シードの住所で生成種を選べるのは中間強度だけ"
         "——選択性は δ≈10⁻³ で最大になり、飽和すると消える")


# ===================== figA19 六判定量のトレードオフ =====================
def figA19(out):
    labels = ["構成できた", "空間が生まれた", "物質が生まれた", "時間が生まれた",
              "3次元が確定 (align>0.8)", "凝縮体 (閉塞<1e-5)"]
    fig, axes = plt.subplots(1, 5, figsize=(17.5, 3.8), sharey=True)
    grid = {}
    for a, mode in zip(axes, MODES):
        Mx = np.zeros((6, len(DELTAS)))
        for j, d in enumerate(DELTAS):
            p = result_path(mode, d)
            if not p.exists():
                continue
            r = json.loads(p.read_text())["N"]["12"]
            Mx[0, j] = 1.0
            Mx[1, j] = 1.0 if r.get("space_born") else 0.0
            Mx[2, j] = 1.0 if r.get("matter_born") else 0.0
            Mx[3, j] = 1.0 if r.get("time_born") else 0.0
            al = r.get("align_med")
            Mx[4, j] = 1.0 if (al is not None and np.isfinite(al) and al > 0.8) else 0.0
            cc = r.get("cond_closure_med")
            Mx[5, j] = 1.0 if (cc is not None and np.isfinite(cc) and cc < 1e-5) else 0.0
        grid[mode] = Mx.tolist()
        a.imshow(Mx, aspect="auto", cmap="RdYlGn", vmin=0, vmax=1)
        for i in range(6):
            for j in range(len(DELTAS)):
                a.text(j, i, "○" if Mx[i, j] > 0.5 else "×", ha="center",
                       va="center", fontsize=11,
                       color="black" if Mx[i, j] > 0.5 else "white")
        a.axhline(3.5, color="k", lw=2.0)
        a.set_xticks(range(len(DELTAS)))
        a.set_xticklabels([f"{d:g}" for d in DELTAS], rotation=60, fontsize=8)
        a.set_title(LABEL[mode], fontsize=10)
        a.set_xlabel("シード強度 δ")
    axes[0].set_yticks(range(6))
    axes[0].set_yticklabels(labels, fontsize=9)
    out["birth_grid"] = grid
    fig.text(0.5, 0.015,
             "太い横線より上の 4 つは強度とともに増えるが、下の 2 つは強度とともに失われる"
             "——単調な「誕生」ではなくトレードオフである",
             ha="center", fontsize=11,
             bbox=dict(fc="lightyellow", ec="gray"))
    save(fig, "figA19_six_criteria_tradeoff_v1.png",
         "図A19　判定量は 6 つある。強いシードは物質と時計を生む代わりに、"
         "3 次元の確定と凝縮体を壊す", tight=(0, 0.10, 1, 0.92))


if __name__ == "__main__":
    print("=== 監査由来の主張図（4 枚）===")
    out = {"generator": {"script": Path(__file__).name,
                         "sha256": sha256(Path(__file__).resolve())},
           "declared": {"sat_level": SAT_LEVEL,
                        "sat_levels_sensitivity": list(SAT_LEVELS_SENS),
                        "birth_criteria": "母体 run_tb_nsweep_1to20_v1.fig_matrix と同一"
                                          "（align_med>0.8 / cond_closure_med<1e-5）"},
           "source": "分析記録_全図面監査_v2.md"}
    pts, sens = collect_tau_sat()
    out["tau_sat_points"] = pts
    out["tau_sat_sensitivity"] = sens
    print(f"  τ_sat を取得できた条件: {len(pts)} 件"
          f"（うち長時間走行由来 {sum(1 for p in pts if p['source']=='T300000')} 件）")
    figA16(pts, out)
    figA17(pts, out)
    figA18(out)
    figA19(out)
    (HERE / "result_paperA_audit_figs_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False, default=str))
    print("完了：4 枚 → result_paperA_audit_figs_v1.json")
