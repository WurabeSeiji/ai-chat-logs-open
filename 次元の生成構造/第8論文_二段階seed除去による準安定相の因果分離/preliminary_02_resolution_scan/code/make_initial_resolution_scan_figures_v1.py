#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第2予備実験 §13 必須図（指定図のみ機械作成）。解釈しない。図中に結果テキストを入れない。

13.1 各run: a_outside(step) 対数縦軸（0値は欠測扱い・描画せず, 回帰線なし）。exec1 を使用。
13.2 3N×5p=15枚: 同N同p で 5 Delta_ref 重ね。
13.3 5p×5Delta_ref=25枚: 同p同Delta_ref で N=5,40,300 重ね。
13.4 13.2/13.3 の固定軸版（step 0〜10000, a 1e-16〜1e-2, 対数）。データ範囲外でも軸不変。
"""
import csv
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CODE = Path(__file__).resolve().parent
P2 = CODE.parent
RAW = P2 / "raw"
FIG = P2 / "figures"; FIG.mkdir(exist_ok=True)
NS = [5, 40, 300]
PS = [0.0, 0.5, 1.0, 1.5, 2.0]
DREFS = [1e-4, 1e-6, 1e-8, 1e-10, 1e-12]
DCOL = {1e-4: "#4C78A8", 1e-6: "#F58518", 1e-8: "#54A24B", 1e-10: "#E45756", 1e-12: "#B279A2"}
NCOL = {5: "#4C78A8", 40: "#F58518", 300: "#E45756"}
FIXED_X = (0, 10000); FIXED_Y = (1e-16, 1e-2)


def rid(n, p, dref, ex=1):
    return f"N{n:05d}_p{p:.1f}_dref{dref:.0e}_exec{ex}"


def base_rid(n, ex=1):
    return f"N{n:05d}_baseline_exec{ex}"


def series(run_id):
    d = RAW / run_id / "timeseries.csv"
    if not d.exists():
        return None
    rows = list(csv.DictReader(open(d)))
    t = np.array([int(r["step"]) for r in rows], float)
    a = np.array([float(r["a_outside"]) for r in rows], float)
    m = a > 0                          # 0値は欠測扱い（描画せず）
    return t[m], a[m]


def fig_per_run(run_id):
    s = series(run_id)
    if s is None:
        return
    t, a = s
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.semilogy(t, a, lw=0.9, color="#333333")
    ax.set_xlabel("step"); ax.set_ylabel("a_outside")
    ax.set_title(run_id); ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout(); fig.savefig(FIG / f"per_run_{run_id}.png", dpi=120); plt.close(fig)


def fig_compare_Np(n, p, fixed=False):
    fig, ax = plt.subplots(figsize=(9, 5))
    for dref in DREFS:
        s = series(rid(n, p, dref))
        if s is None:
            continue
        ax.semilogy(s[0], s[1], lw=0.9, color=DCOL[dref], label=f"Δ_ref={dref:.0e}")
    ax.set_xlabel("step"); ax.set_ylabel("a_outside")
    tag = "fixedaxis" if fixed else "auto"
    ax.set_title(f"N={n} p={p:.1f}  (5 Delta_ref)  [{tag}]"); ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.25)
    if fixed:
        ax.set_xlim(*FIXED_X); ax.set_ylim(*FIXED_Y)
    fig.tight_layout()
    fig.savefig(FIG / f"cmp_Np_N{n:05d}_p{p:.1f}_{tag}.png", dpi=120); plt.close(fig)


def fig_compare_pD(p, dref, fixed=False):
    fig, ax = plt.subplots(figsize=(9, 5))
    for n in NS:
        s = series(rid(n, p, dref))
        if s is None:
            continue
        ax.semilogy(s[0], s[1], lw=0.9, color=NCOL[n], label=f"N={n}")
    ax.set_xlabel("step"); ax.set_ylabel("a_outside")
    tag = "fixedaxis" if fixed else "auto"
    ax.set_title(f"p={p:.1f} Δ_ref={dref:.0e}  (N=5,40,300)  [{tag}]"); ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.25)
    if fixed:
        ax.set_xlim(*FIXED_X); ax.set_ylim(*FIXED_Y)
    fig.tight_layout()
    fig.savefig(FIG / f"cmp_pD_p{p:.1f}_dref{dref:.0e}_{tag}.png", dpi=120); plt.close(fig)


def main():
    # 13.1 各run（ON 全 config + baseline）exec1
    n_pr = 0
    for n in NS:
        for p in PS:
            for dref in DREFS:
                fig_per_run(rid(n, p, dref)); n_pr += 1
        fig_per_run(base_rid(n)); n_pr += 1
    # 13.2 + 13.4
    for n in NS:
        for p in PS:
            fig_compare_Np(n, p, fixed=False); fig_compare_Np(n, p, fixed=True)
    # 13.3 + 13.4
    for p in PS:
        for dref in DREFS:
            fig_compare_pD(p, dref, fixed=False); fig_compare_pD(p, dref, fixed=True)
    print(f"[figs] per_run={n_pr}, cmp_Np={len(NS)*len(PS)}×2, cmp_pD={len(PS)*len(DREFS)}×2 → {FIG}")


if __name__ == "__main__":
    main()
