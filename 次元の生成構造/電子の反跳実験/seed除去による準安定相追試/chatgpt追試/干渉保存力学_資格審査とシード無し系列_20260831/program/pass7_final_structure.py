#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス7（本フレーム追加分析・読出し専用）：最終状態の構造検証。
仮説（星型アンカー）：等振幅の頂点星（1 頂点の N−1 辺、z² 整列、符号は Z₂ 枝）は新力学の厳密平衡で
  μ_star=(N−2)a²=N(N−2)r̄²/2、PR=N−1。各走行の最終 step 状態から検証する。"""
import os, csv, json
import numpy as np
import sys
sys.path.insert(0, os.path.dirname(__file__))
from common import edges, adjacency
from interference_dynamics import hermitian_H
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
def fmt(x, d=4):
    return "—" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{d}f}"
def fe(x):
    return "—" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.2e}"
rows = []
for tag in sorted(os.listdir(DATA)):
    if tag == "reference": continue
    fp = os.path.join(DATA, tag, "states_treatment.npz")
    cj = os.path.join(DATA, tag, "parent_checks.json")
    if not (os.path.exists(fp) and os.path.exists(cj)): continue
    c = json.load(open(cj)); N = c["N"]; E = edges(N); M = len(E); A = adjacency(N)
    Z = np.load(fp)["Z"][-1]
    H = float(np.vdot(Z, Z).real); r2 = H/M
    Hz = hermitian_H(Z, A) @ Z
    mu = float((np.vdot(Z, Hz)/np.vdot(Z, Z)).real)
    res = float(np.linalg.norm(Hz - mu*Z)/np.linalg.norm(Z))/r2
    W = np.zeros(N)
    for k, (i, j) in enumerate(E):
        W[i] += abs(Z[k])**2; W[j] += abs(Z[k])**2
    istar = int(np.argmax(W)); star = [k for k, e in enumerate(E) if istar in e]
    star_frac = float(sum(abs(Z[k])**2 for k in star)/H)
    amps = np.array([abs(Z[k]) for k in star])
    amp_unif = float(np.ptp(amps)/amps.mean()) if amps.mean() > 0 else float("nan")
    phs = np.array([np.angle(Z[k]) for k in star])
    g = np.angle(np.sum(np.exp(2j*phs)))/2.0
    n_plus = int(np.sum(np.cos(phs - g) > 0)); n_minus = len(phs) - n_plus
    z2_align = float(np.abs(np.sum(Z[star]**2))/max(np.sum(np.abs(Z[star])**2), 1e-300))
    rows.append(dict(tag=tag, N=N, method=c["method"], star_vertex=istar, star_frac=star_frac,
                     amp_uniformity=amp_unif, n_plus=n_plus, n_minus=n_minus,
                     z2_alignment=z2_align, mu_final_over_r2=mu/r2, star_mu_pred=N*(N-2)/2.0,
                     mu_dev=abs(mu/r2 - N*(N-2)/2.0)/(N*(N-2)/2.0),
                     residual_new_over_r2=res, landed=bool(res < 1e-6)))
with open(os.path.join(ROOT, "results", "final_state_structure.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
lines = ["# 最終状態の構造（星型アンカー検証・本フレーム追加分析）", "",
         "| tag | 星頂点 | 星純度 | 振幅均等(ptp/mean) | ±内訳 | z²整列 | μ_fin/r̄² | N(N−2)/2 | 相対差 | res_new/r̄² | 着地 |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
for r in rows:
    lines.append("| " + " | ".join([r["tag"], str(r["star_vertex"]), fmt(r["star_frac"]),
        fmt(r["amp_uniformity"], 3), str(r["n_plus"]) + "/" + str(r["n_minus"]), fmt(r["z2_alignment"], 6),
        fmt(r["mu_final_over_r2"], 3), fmt(r["star_mu_pred"], 1), fe(r["mu_dev"]),
        fe(r["residual_new_over_r2"]), "YES" if r["landed"] else "no"]) + " |")
with open(os.path.join(ROOT, "results", "final_state_structure.md"), "w") as f:
    f.write("\n".join(lines) + "\n")
print(f"PASS7 OK（{len(rows)} 走行、results/final_state_structure.md）")
