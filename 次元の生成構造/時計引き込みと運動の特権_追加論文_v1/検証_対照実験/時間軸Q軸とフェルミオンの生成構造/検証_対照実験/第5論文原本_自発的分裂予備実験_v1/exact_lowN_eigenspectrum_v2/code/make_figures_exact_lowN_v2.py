#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第7論文 第1段階修正版(v2) N=5 図（§10 の13図）。解釈なし。N=5のみ。"""
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CODE = Path(__file__).resolve().parent
BASE = CODE.parent


def rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def pivot(path, key, val, tcol="time"):
    d = defaultdict(dict); keys = set()
    for r in rows(path):
        t = int(float(r[tcol])); k = int(float(r[key]))
        d[t][k] = float(r[val]); keys.add(k)
    ts = sorted(d); keys = sorted(keys)
    M = np.full((len(ts), len(keys)), np.nan)
    ki = {k: i for i, k in enumerate(keys)}
    for i, t in enumerate(ts):
        for k, v in d[t].items():
            M[i, ki[k]] = v
    return np.array(ts), keys, M


def make(n=5):
    raw = BASE / "raw" / f"N{n:05d}"
    fd = BASE / "figures" / f"N{n:05d}"; fd.mkdir(parents=True, exist_ok=True)
    cr = json.load(open(BASE / "diagnostics" / f"N{n:05d}.json"))["crossing"]
    ev = raw / "eigenvalues.csv"

    # 図1: 全 μ_j 線形・対数
    t, idx, MU = pivot(ev, "eigen_index", "eigenvalue")
    for logy, tag in [(False, "01_mu_linear"), (True, "01_mu_abslog")]:
        fig, ax = plt.subplots(figsize=(10, 6))
        for j in range(MU.shape[1]):
            (ax.semilogy if logy else ax.plot)(t, np.abs(MU[:, j]) if logy else MU[:, j], lw=0.6)
        ax.axvline(cr, color="k", ls=":", lw=0.8); ax.set_xlabel("time"); ax.set_ylabel("mu_j" + (" |.|log" if logy else ""))
        ax.set_title(f"N={n} Fig{tag}: all eigenvalues mu_j=eig(iK)")
        fig.tight_layout(); fig.savefig(fd / f"fig{tag}.png", dpi=130); plt.close(fig)

    # 図2: |σ_j|/σ1 線形・対数
    _, _, R = pivot(ev, "eigen_index", "sigma_over_sigma1")
    for logy, tag in [(False, "02_ratio_linear"), (True, "02_ratio_log")]:
        fig, ax = plt.subplots(figsize=(10, 6))
        for j in range(R.shape[1]):
            (ax.semilogy if logy else ax.plot)(t, np.clip(R[:, j], 1e-20, None) if logy else R[:, j], lw=0.6)
        ax.axvline(cr, color="k", ls=":", lw=0.8); ax.set_xlabel("time"); ax.set_ylabel("|sigma_j|/sigma_1")
        ax.set_title(f"N={n} Fig{tag}"); fig.tight_layout(); fig.savefig(fd / f"fig{tag}.png", dpi=130); plt.close(fig)

    # 図3,4: N, N²倍
    for col, tag in [("N_sigma_over_sigma1", "03_N_ratio"), ("N2_sigma_over_sigma1", "04_N2_ratio")]:
        _, _, RR = pivot(ev, "eigen_index", col)
        fig, ax = plt.subplots(figsize=(10, 6))
        for j in range(RR.shape[1]):
            ax.plot(t, RR[:, j], lw=0.6)
        ax.axvline(cr, color="k", ls=":", lw=0.8); ax.set_xlabel("time"); ax.set_ylabel(tag)
        ax.set_title(f"N={n} Fig{tag}"); fig.tight_layout(); fig.savefig(fd / f"fig{tag}.png", dpi=130); plt.close(fig)

    # 図5: 数値床比 σ/(eps|K|)
    _, _, FL = pivot(ev, "eigen_index", "abs_sigma_over_eps_normK")
    fig, ax = plt.subplots(figsize=(10, 6))
    for j in range(FL.shape[1]):
        ax.semilogy(t, np.clip(FL[:, j], 1e-2, None), lw=0.6)
    ax.axhline(1000, color="r", ls="--", lw=0.8, label="1000 (rot primary)")
    ax.axhline(100, color="orange", ls="--", lw=0.8, label="100 (alt)")
    ax.axvline(cr, color="k", ls=":", lw=0.8); ax.legend(); ax.set_xlabel("time")
    ax.set_ylabel("|sigma_j|/(eps |K|)"); ax.set_title(f"N={n} Fig05: distance from numerical floor")
    fig.tight_layout(); fig.savefig(fd / "fig05_floor_ratio.png", dpi=130); plt.close(fig)

    # 図6: 正負対誤差 / 図13: 平面間直交 / (diagnostics)
    dg = rows(raw / "diagnostics_timeseries.csv")
    dt = np.array([int(float(r["time"])) for r in dg])
    pe = np.array([float(r["pair_error"]) for r in dg])
    ip = np.array([float(r["max_interplane_nondeg"]) for r in dg])
    ci = np.array([float(r["max_cluster_interplane"]) for r in dg])
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.semilogy(dt, np.clip(pe, 1e-20, None), lw=0.8); ax.axvline(cr, color="k", ls=":", lw=0.8)
    ax.set_xlabel("time"); ax.set_ylabel("pair error max|mu_k+mu_{M-1-k}|")
    ax.set_title(f"N={n} Fig06: positive-negative pair error")
    fig.tight_layout(); fig.savefig(fd / "fig06_pair_error.png", dpi=130); plt.close(fig)
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.semilogy(dt, np.clip(ip, 1e-20, None), lw=0.8, label="non-degenerate interplane")
    ax.semilogy(dt, np.clip(ci, 1e-20, None), lw=0.8, ls="--", label="cluster interplane")
    ax.axhline(1e-12, color="r", ls=":", lw=0.8, label="1e-12 acceptance")
    ax.axvline(cr, color="k", ls=":", lw=0.8); ax.legend(); ax.set_xlabel("time")
    ax.set_ylabel("max ||B_j^T B_k||_2"); ax.set_title(f"N={n} Fig13: max interplane orthogonality error")
    fig.tight_layout(); fig.savefig(fd / "fig13_interplane_error.png", dpi=130); plt.close(fig)

    # 図7: 縮退クラスタ時間推移（mult>1 クラスタ数、クラスタ数）
    cl = rows(raw / "clusters.csv")
    bytime = defaultdict(lambda: {"n": 0, "deg": 0})
    for r in cl:
        tt = int(float(r["time"])); bytime[tt]["n"] += 1
        if int(r["mult"]) > 1:
            bytime[tt]["deg"] += 1
    ts = sorted(bytime)
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(ts, [bytime[x]["n"] for x in ts], lw=0.9, label="n_clusters")
    ax.plot(ts, [bytime[x]["deg"] for x in ts], lw=0.9, label="n_degenerate_clusters(mult>1)")
    ax.axvline(cr, color="k", ls=":", lw=0.8); ax.legend(); ax.set_xlabel("time")
    ax.set_ylabel("count"); ax.set_title(f"N={n} Fig07: degenerate clusters over time")
    fig.tight_layout(); fig.savefig(fd / "fig07_degenerate_clusters.png", dpi=130); plt.close(fig)

    # 図8: q1..q4
    q = rows(raw / "q_svd.csv")
    qt = np.array([int(float(r["time"])) for r in q])
    Q = np.array([[float(r["q1"]), float(r["q2"]), float(r["q3"]), float(r["q4"])] for r in q])
    fig, ax = plt.subplots(figsize=(10, 6))
    for j, lab in enumerate(["q1", "q2", "q3", "q4"]):
        ax.plot(qt, Q[:, j], lw=0.8, label=lab)
    ax.axvline(cr, color="k", ls=":", lw=0.8); ax.legend(); ax.set_xlabel("time")
    ax.set_ylabel("sv of Q=[B0|Bdom]"); ax.set_title(f"N={n} Fig08: four-basis singular values")
    fig.tight_layout(); fig.savefig(fd / "fig08_q_svd.png", dpi=130); plt.close(fig)

    # 図9: 初期床から非零化した branch の σ/σ1 推移
    brtime = defaultdict(dict); brflag = {}
    for r in cl:
        b = int(r["cluster_branch"]); tt = int(float(r["time"]))
        brtime[b][tt] = float(r["sigma_over_sigma1"]); brflag[b] = int(r["initial_floor_flag"])
    fig, ax = plt.subplots(figsize=(10, 6))
    grew = [b for b in brtime if brflag.get(b, 0) == 1]
    for b in (grew or list(brtime)):
        xs = sorted(brtime[b]); ax.plot(xs, [brtime[b][x] for x in xs], lw=0.7)
    ax.axvline(cr, color="k", ls=":", lw=0.8); ax.set_xlabel("time"); ax.set_ylabel("sigma/sigma1")
    ax.set_title(f"N={n} Fig09: branches that appeared after initial (initial_floor_flag=1): {len(grew)}")
    fig.tight_layout(); fig.savefig(fd / "fig09_floor_branches.png", dpi=130); plt.close(fig)

    # 図10: 対象 branch ごとの δ
    fig, ax = plt.subplots(figsize=(10, 6))
    dbt = defaultdict(dict)
    for r in cl:
        dbt[int(r["cluster_branch"])][int(float(r["time"]))] = float(r["delta_C"])
    for b in dbt:
        xs = sorted(dbt[b]); ax.plot(xs, [dbt[b][x] for x in xs], lw=0.6)
    ax.axvline(cr, color="k", ls=":", lw=0.8); ax.set_xlabel("time"); ax.set_ylabel("delta_C per branch")
    ax.set_title(f"N={n} Fig10: delta per cluster branch"); fig.tight_layout()
    fig.savefig(fd / "fig10_delta_branches.png", dpi=130); plt.close(fig)

    # 図11: 対象別残差行列 全特異値
    tg = rows(raw / "delta_targets.csv")
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    for ax, ttype in zip(axes.flat, ["A_dominant", "C_occupied", "D_degenerate", "all_diagnostic"]):
        dd = defaultdict(dict)
        for r in tg:
            if r["target_type"] == ttype:
                dd[int(float(r["time"]))][int(r["k"])] = float(r["singular_value"])
        ts2 = sorted(dd)
        if ts2:
            kmax = max(max(dd[x]) for x in ts2) + 1
            for k in range(kmax):
                ax.semilogy(ts2, [np.clip(dd[x].get(k, np.nan), 1e-20, None) for x in ts2], lw=0.6)
        ax.axvline(cr, color="k", ls=":", lw=0.8); ax.set_title(ttype); ax.set_xlabel("time")
        ax.set_ylabel("singular values")
    fig.suptitle(f"N={n} Fig11: target residual (I-Pi0)B_target singular values")
    fig.tight_layout(); fig.savefig(fd / "fig11_target_residual_svd.png", dpi=130); plt.close(fig)

    # 図12: クラスタ占有 E_C 線形・対数
    for logy, tag in [(False, "12_occ_linear"), (True, "12_occ_log")]:
        oc = defaultdict(dict)
        for r in cl:
            oc[int(r["cluster_branch"])][int(float(r["time"]))] = float(r["occupation_fraction"])
        fig, ax = plt.subplots(figsize=(10, 6))
        for b in oc:
            xs = sorted(oc[b]); ys = [oc[b][x] for x in xs]
            (ax.semilogy if logy else ax.plot)(xs, np.clip(ys, 1e-20, None) if logy else ys, lw=0.6)
        ax.axvline(cr, color="k", ls=":", lw=0.8); ax.set_xlabel("time"); ax.set_ylabel("E_C / |Z|^2")
        ax.set_title(f"N={n} Fig{tag}: cluster occupation")
        fig.tight_layout(); fig.savefig(fd / f"fig{tag}.png", dpi=130); plt.close(fig)

    print(f"[v2 figs] N={n}: §10 の図を {fd} に出力")


if __name__ == "__main__":
    make(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
