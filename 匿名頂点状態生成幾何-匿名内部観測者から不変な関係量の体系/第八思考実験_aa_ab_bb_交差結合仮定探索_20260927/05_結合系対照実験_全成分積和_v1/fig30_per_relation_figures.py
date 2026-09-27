#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""結合系の観測 row データから、関係ごとの図を描く。

  python3 fig30_per_relation_figures.py <pattern_id>

図を描く部分は、第六思考実験の postprocess_strict_charged8_5cases.py の plot_case と
sample_raw を写したものである。元のプログラムは paper6_figure_programs/ に無変更で置いてある。
元のプログラムと違う点
  1. 入力は、結合系の row を第六の行の形式に並べ替えたもの（fig_common.paper6_rows）。
  2. 保存先。
  3. figure05（最初の macrostep の work state）は描かない。
     結合系の実行では microstep の途中の状態を保存していないため。
図に描く点は、第六と同じく全行から約 8000 点を等間隔に取り出したもの。
"""
from __future__ import annotations
import json
import math
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fig_common import case_id, figure_dir, paper6_rows, read_unified_rows, relation_conditions

HSTEP = 2 * math.pi / 4000


def sample_raw(g, total_rows, limit=8000):
    """第六の sample_raw と同じ取り出し方（step が stride の倍数の行と、最後の行）。"""
    stride = max(1, total_rows // limit)
    st = g[:, 0].astype(np.int64)
    m = (st % stride) == 0
    out = [g[m]]
    last = g[-1:]
    if int(out[-1][-1, 0]) != int(last[0, 0]):
        out.append(last)
    return np.concatenate(out), stride


def save(fig, figdir, name):
    fig.tight_layout()
    fig.savefig(figdir / (name + ".svg"))
    fig.savefig(figdir / (name + ".png"), dpi=180)
    plt.close(fig)


def plot_case(cid, g_all, anabase, figdir):
    figdir.mkdir(parents=True, exist_ok=True)
    ref = pd.read_csv(anabase / cid / "charged_binary_analytic_reference_v1_raw.csv",
                      usecols=["r", "t", "phi", "x_rel", "y_rel"])
    rr = ref["r"].to_numpy(); tr = ref["t"].to_numpy(); pr = ref["phi"].to_numpy()
    xr = ref["x_rel"].to_numpy(); yr = ref["y_rel"].to_numpy()
    ir = np.linspace(0, len(ref) - 1, min(8000, len(ref)), dtype=int)
    g, stride = sample_raw(g_all, g_all.shape[0])
    st = g[:, 0]; r = g[:, 1]; t = g[:, 11]; x = g[:, 12]; y = g[:, 13]

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot(xr[ir], yr[ir], label="Analytic reference", linewidth=1.1)
    ax.plot(x, y, "--", label="Strict state generator", linewidth=.9)
    ax.set_xlabel("x / M"); ax.set_ylabel("y / M"); ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=.25); ax.legend(); ax.set_title(cid + ": orbit overlay")
    save(fig, figdir, "figure01_orbit_overlay")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(tr[ir], rr[ir], label="Analytic reference")
    mt = np.isfinite(t)
    ax.plot(t[mt], r[mt], "--", label="Generator (finite t readout)")
    ax.set_xlabel("t / M"); ax.set_ylabel("r / M"); ax.grid(True, alpha=.25); ax.legend()
    ax.set_title(cid + ": radius vs time")
    save(fig, figdir, "figure02_radius_vs_time")

    mt = np.isfinite(t) & (t <= tr[-1])
    fig, ax = plt.subplots(figsize=(8, 5))
    if np.any(mt):
        ri = np.interp(t[mt], tr, rr); pi = np.interp(t[mt], tr, pr)
        er = np.abs(r[mt] - ri); exy = np.hypot(x[mt] - ri * np.cos(pi), y[mt] - ri * np.sin(pi))
        ax.semilogy(t[mt], np.maximum(er, 1e-18), label="|Delta r| at same t")
        ax.semilogy(t[mt], np.maximum(exy, 1e-18), label="||Delta(x,y)|| at same t")
    ax.set_xlabel("t / M"); ax.set_ylabel("absolute residual"); ax.grid(True, alpha=.25); ax.legend()
    ax.set_title(cid + ": finite-time residuals")
    save(fig, figdir, "figure03_reference_residuals")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(st, g[:, 8], label="N"); ax.plot(st, g[:, 9], label="C"); ax.plot(st, g[:, 10], label="D")
    ax.set_xlabel("macro step"); ax.set_ylabel("state value"); ax.grid(True, alpha=.25); ax.legend()
    ax.set_title(cid + ": identity-propagated states")
    save(fig, figdir, "figure04_identity_states")

    ph = st * HSTEP; mp = ph <= pr[-1] + 1e-12
    ri = np.interp(ph[mp], pr, rr)
    er = np.abs(r[mp] - ri); exy = np.hypot(x[mp] - ri * np.cos(ph[mp]), y[mp] - ri * np.sin(ph[mp]))
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.semilogy(st[mp] / 4000, np.maximum(er, 1e-18), label="|Delta r| at same phase")
    ax.semilogy(st[mp] / 4000, np.maximum(exy, 1e-18), label="||Delta(x,y)|| at same phase")
    ax.set_xlabel("orbital cycles"); ax.set_ylabel("absolute residual"); ax.grid(True, alpha=.25); ax.legend()
    ax.set_title(cid + ": full-run phase-domain residuals")
    save(fig, figdir, "figure06_full_phase_residuals")

    return {"case": cid, "rows": int(g_all.shape[0]), "rows_drawn": int(g.shape[0]), "stride": int(stride),
            "P_first": float(g_all[0, 1]), "P_last": float(g_all[-1, 1]),
            "points_within_reference_time_range": int(mt.sum()),
            "points_within_reference_phase_range": int(mp.sum())}


def main() -> None:
    pattern_id = sys.argv[1]
    base = figure_dir(pattern_id)
    data, columns = read_unified_rows(pattern_id)
    records = []
    for rel, cid, la, lb, c0, d0 in relation_conditions(pattern_id):
        g_all = paper6_rows(data, columns, rel)
        rec = plot_case(cid, g_all, base / "analytic_reference", base / cid)
        rec.update({"relation": rel, "C": c0, "D": d0})
        records.append(rec)
        print(cid, rec, flush=True)
    (base / "per_relation_figure_record.json").write_text(json.dumps(records, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
