#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス4：図（事実の記録のみ）。fig1 閉塞率、fig2 重なり欠損、fig3 時計、fig4 最終状態（PR/M と残差）。"""
import os, sys, json, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = ["Hiragino Sans", "Arial Unicode MS", "sans-serif"]
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data"); FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)

def load_ts(tag):
    p = os.path.join(DATA, tag, "timeseries.csv")
    if not os.path.exists(p): return None
    return np.genfromtxt(p, delimiter=",", names=True)

REP = ["hm_N8", "ne_N8", "rb_N8", "hm_N12", "ne_N12", "rb_N12"]
COL = {"hm": "tab:blue", "ne": "tab:red", "rb": "tab:green"}

fig, ax = plt.subplots(figsize=(9, 5.5))
for tag in REP:
    ts = load_ts(tag)
    if ts is None: continue
    m = tag.split("_")[0]
    ax.semilogy(ts["step"], np.maximum(ts["closure_frac"], 1e-18),
                color=COL[m], ls=("-" if "N8" in tag else "--"), label=tag)
ax.set_xlabel("step"); ax.set_ylabel("|ΣZ²|/H（閉塞率）")
ax.set_title("閉塞は力学量：アンカー（hm）は保持、非平衡（ne/rb）は開く")
ax.legend(fontsize=8); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig1_closure_frac.png"), dpi=150); plt.close(fig)

fig, ax = plt.subplots(figsize=(9, 5.5))
for tag in REP:
    ts = load_ts(tag)
    if ts is None: continue
    m = tag.split("_")[0]
    ax.semilogy(ts["step"], np.maximum(ts["overlap_deficit"], 1e-18),
                color=COL[m], ls=("-" if "N8" in tag else "--"), label=tag)
ax.set_xlabel("step"); ax.set_ylabel("1 − |⟨v,Z⟩|²/(‖v‖²‖Z‖²)")
ax.set_title("親からの離脱（重なり欠損）")
ax.legend(fontsize=8); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig2_overlap_deficit.png"), dpi=150); plt.close(fig)

rows = list(csv.DictReader(open(os.path.join(ROOT, "results", "matrix_N_by_method.csv"))))
fig, ax = plt.subplots(figsize=(9, 5))
for m in ("hm", "ne", "rb"):
    xs = [int(r["N"]) for r in rows if r["method"] == m and r["clock_ratio_first200"] not in ("", "nan")]
    ys = [float(r["clock_ratio_first200"]) for r in rows if r["method"] == m and r["clock_ratio_first200"] not in ("", "nan")]
    ax.plot(xs, ys, "o", color=COL[m], label=m)
ax.axhline(1.0, color="k", lw=0.8)
ax.set_xlabel("N"); ax.set_ylabel("実測位相進み / 予測（−Δμ_new）")
ax.set_title("時計の検証（最初の 200 step 平均）")
ax.legend(); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig3_clock.png"), dpi=150); plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for m in ("hm", "ne", "rb"):
    sel = [r for r in rows if r["method"] == m]
    axes[0].plot([int(r["N"]) for r in sel], [float(r["PR_over_M_final"]) for r in sel], "o-", color=COL[m], label=m)
    axes[1].semilogy([int(r["N"]) for r in sel], [max(float(r["final_residual_new_over_r2"]), 1e-18) for r in sel], "o-", color=COL[m], label=m)
axes[0].set_xlabel("N"); axes[0].set_ylabel("PR/M（最終）"); axes[0].set_title("最終状態の局在"); axes[0].legend(); axes[0].grid(alpha=0.3)
axes[1].set_xlabel("N"); axes[1].set_ylabel("res_new/r²（最終）"); axes[1].set_title("新アンカーへの着地（残差）"); axes[1].legend(); axes[1].grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig4_final_states.png"), dpi=150); plt.close(fig)
print("PASS4 OK（figures/fig1..fig4）")
