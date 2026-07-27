#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第8論文 第1予備実験 図（§10 の fig1〜4）＋集計表（§11）。解釈なし。指示外の図は作らない。

fig1: f_outside_parent（A/B/D）
fig2: q3, q4（A/B/D）
fig3: rank_Q（A/B/D）
fig4: 準安定域の代表観測量（f_outside_parent）B vs D（絶対step併記）
集計: summary/preliminary_seed_ablation_summary.csv
共通横軸 絶対step 0..55000, 5000刻み, crossing 不動。
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CODE = Path(__file__).resolve().parent
PAPER8 = CODE.parent
FIG = PAPER8 / "figures"; FIG.mkdir(exist_ok=True)
NS = [5, 40, 300]
XMAX = 55000
FILES = {"A": "condition_A_no_seed", "B": "condition_B_initial_only", "D": "condition_D_existing_two_seed"}
COL = {"A": "#7F7F7F", "B": "#F58518", "D": "#E45756"}
REP_OBS = "f_outside_parent"     # 第7論文の代表観測量（分裂量）


def load(n, c):
    p = PAPER8 / "raw" / f"N{n:05d}" / f"{FILES[c]}.csv"
    if not p.exists():
        return None
    return list(csv.DictReader(open(p)))


def diag(n, c):
    p = PAPER8 / "diagnostics" / f"N{n:05d}_condition_{c}.json"
    return json.load(open(p)) if p.exists() else {}


def col(rows, name):
    return np.array([float(r[name]) for r in rows])


def setx(ax):
    ax.set_xlim(0, XMAX); ax.set_xticks(np.arange(0, XMAX + 1, 5000))


def ready(n):
    return all(load(n, c) is not None for c in "ABD")


def make_all():
    avail = [n for n in NS if ready(n)]
    for n in avail:
        cr = diag(n, "B").get("crossing_step")
        t1 = diag(n, "D").get("injection", {}).get("injected_at")
        rows = {c: load(n, c) for c in "ABD"}
        t = {c: col(rows[c], "time") for c in "ABD"}
        # fig1 f
        fig, ax = plt.subplots(figsize=(10, 5))
        for c in "ABD":
            ax.plot(t[c], col(rows[c], "f_outside_parent"), lw=0.9, color=COL[c], label=f"cond {c}")
        if cr is not None:
            ax.axvline(cr, color="k", ls=":", lw=0.8)
        setx(ax); ax.set_ylim(0, 1); ax.set_xlabel("step (absolute)"); ax.set_ylabel("f_outside_parent = 1 - E_P1")
        ax.set_title(f"N={n} fig01: f_outside_parent (A/B/D)"); ax.legend()
        fig.tight_layout(); fig.savefig(FIG / f"fig01_f_compare_N{n:05d}.png", dpi=130); plt.close(fig)
        # fig2 q3,q4
        fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
        for c in "ABD":
            a1.plot(t[c], col(rows[c], "q3"), lw=0.9, color=COL[c], label=f"cond {c}")
            a2.plot(t[c], col(rows[c], "q4"), lw=0.9, color=COL[c])
        for a in (a1, a2):
            if cr is not None: a.axvline(cr, color="k", ls=":", lw=0.8)
            setx(a)
        a1.set_ylabel("q3"); a2.set_ylabel("q4"); a2.set_xlabel("step (absolute)"); a1.legend()
        a1.set_title(f"N={n} fig02: q3, q4 (A/B/D)")
        fig.tight_layout(); fig.savefig(FIG / f"fig02_q3q4_compare_N{n:05d}.png", dpi=130); plt.close(fig)
        # fig3 rank_Q
        fig, ax = plt.subplots(figsize=(10, 5))
        for c in "ABD":
            ax.plot(t[c], col(rows[c], "rank_Q"), lw=0.9, color=COL[c], label=f"cond {c}")
        if cr is not None: ax.axvline(cr, color="k", ls=":", lw=0.8)
        setx(ax); ax.set_ylim(-0.2, 4.5); ax.set_yticks([0, 2, 4]); ax.set_xlabel("step (absolute)"); ax.set_ylabel("rank_Q")
        ax.set_title(f"N={n} fig03: rank_Q (A/B/D)"); ax.legend()
        fig.tight_layout(); fig.savefig(FIG / f"fig03_rankQ_compare_N{n:05d}.png", dpi=130); plt.close(fig)
        # fig4 metastable B vs D（絶対step併記）
        fig, ax = plt.subplots(figsize=(10, 5))
        for c in "BD":
            ax.plot(t[c], col(rows[c], REP_OBS), lw=0.9, color=COL[c], label=f"cond {c}")
        if t1 is not None:
            ax.axvline(t1, color="g", ls="--", lw=0.9, label=f"t1 (injection) = {t1}")
        if cr is not None: ax.axvline(cr, color="k", ls=":", lw=0.8)
        setx(ax); ax.set_xlabel("step (absolute)"); ax.set_ylabel(REP_OBS)
        ax.set_title(f"N={n} fig04: metastable observable ({REP_OBS}), B vs D"); ax.legend()
        fig.tight_layout(); fig.savefig(FIG / f"fig04_metastable_B_vs_D_N{n:05d}.png", dpi=130); plt.close(fig)
    make_summary(avail)
    print(f"[figs] fig01〜04 を {FIG} に出力（N={avail}）")


def make_summary(avail):
    (PAPER8 / "summary").mkdir(exist_ok=True)
    with open(PAPER8 / "summary" / "preliminary_seed_ablation_summary.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["N", "condition", "initial_seed_enabled", "metastable_seed_enabled",
                    "crossing_detected", "crossing_step", "metastable_start_detected", "metastable_start_step",
                    "max_f", "final_f", "max_q3", "max_q4", "final_q3", "final_q4",
                    "max_rank_Q", "final_rank_Q", "mean_metastable_observable", "std_metastable_observable",
                    "max_norm_error", "max_zero_square_error", "max_projection_closure_error"])
        for n in avail:
            for c in "ABD":
                rows = load(n, c); d = diag(n, c)
                f = col(rows, "f_outside_parent"); q3 = col(rows, "q3"); q4 = col(rows, "q4")
                rq = np.array([int(r["rank_Q"]) for r in rows]); tt = col(rows, "time")
                ms = d.get("metastable_start_step")
                # 準安定域観測量: metastable_start_step 以降の代表観測量(f)。検出なしは NaN
                if ms is not None:
                    win = f[tt >= ms]
                    mean_ms = float(np.mean(win)) if len(win) else float("nan")
                    std_ms = float(np.std(win)) if len(win) else float("nan")
                else:
                    mean_ms = std_ms = float("nan")
                cr = d.get("crossing_step")
                w.writerow([n, c, int(d.get("initial_seed_enabled", False)), int(d.get("metastable_seed_enabled", False)),
                            int(cr is not None), (cr if cr is not None else ""),
                            int(ms is not None), (ms if ms is not None else ""),
                            "%.6e" % np.max(f), "%.6e" % f[-1], "%.6e" % np.max(q3), "%.6e" % np.max(q4),
                            "%.6e" % q3[-1], "%.6e" % q4[-1], int(np.max(rq)), int(rq[-1]),
                            "%.6e" % mean_ms, "%.6e" % std_ms,
                            "%.3e" % d.get("max_norm_error", float("nan")),
                            "%.3e" % d.get("max_zero_square_abs", float("nan")),
                            "%.3e" % d.get("max_projection_closure_error", float("nan"))])
    print("[summary] summary/preliminary_seed_ablation_summary.csv")


if __name__ == "__main__":
    make_all()
