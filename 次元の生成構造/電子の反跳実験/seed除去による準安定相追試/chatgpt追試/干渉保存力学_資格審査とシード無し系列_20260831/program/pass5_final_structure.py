#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス5：最終状態の構造検証（走行後分析）。
仮説（星型アンカー）：等振幅・位相整列の頂点星（1 頂点の N−1 辺）は新力学の厳密平衡で
  μ_star = (N−2)a² = N(N−2)r̄²/2（a² = M r̄²/(N−1)）、PR = N−1、|Σz²|/H = 1。
検証：各走行の最終状態について 最大頂点重み割合・星純度・星上の振幅均等性・位相整列・μ_final/r̄² と N(N−2)/2 の差。
追加：hm の成長率実測（ln closure_frac の傾き/2 ≒ λ_f/2 …閉塞は振幅の 2 次なので傾き=2×成長率）と走行前 ρ 予測の照合。"""
import os, csv, json, math
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
import sys
sys.path.insert(0, os.path.dirname(__file__))
from common import edges

rows = []
for tag in sorted(os.listdir(DATA)):
    fp = os.path.join(DATA, tag, "final_state.npz")
    cj = os.path.join(DATA, tag, "parent_checks.json")
    if not os.path.exists(fp): continue
    c = json.load(open(cj))
    N = c["N"]; E = edges(N); M = len(E)
    Z = np.load(fp)["Z"]
    H = float(np.vdot(Z, Z).real); r2 = H/M
    W = np.zeros(N)
    for k, (i, j) in enumerate(E):
        W[i] += abs(Z[k])**2; W[j] += abs(Z[k])**2
    istar = int(np.argmax(W))
    star_idx = [k for k, e in enumerate(E) if istar in e]
    star_frac = float(sum(abs(Z[k])**2 for k in star_idx)/H)
    amps = np.array([abs(Z[k]) for k in star_idx])
    amp_unif = float(np.ptp(amps)/amps.mean()) if amps.mean() > 0 else float("nan")
    phs = np.array([np.angle(Z[k]) for k in star_idx])
    zsum = np.abs(np.exp(1j*phs).mean())
    # 符号内訳：z² の共通軸 g に対して cos(φ−g) の符号で ± を数える（Z2^M 枝）
    g = np.angle(np.sum(np.exp(2j*phs)))/2.0
    n_plus = int(np.sum(np.cos(phs - g) > 0)); n_minus = len(phs) - n_plus
    z2_align = float(np.abs(np.sum(Z[star_idx]**2))/max(np.sum(np.abs(Z[star_idx])**2), 1e-300))
    mu_over_r2 = None
    sj = json.load(open(os.path.join(DATA, tag, "summary.json")))
    mu_over_r2 = sj["final_mu_new_over_r2"]
    star_mu = N*(N-2)/2.0
    ts = np.genfromtxt(os.path.join(DATA, tag, "timeseries.csv"), delimiter=",", names=True)
    y = ts["closure_frac"]; mask = (y > 1e-12) & (y < 1e-4)
    lam_meas = None
    if mask.sum() > 100:
        sl = np.polyfit(ts["step"][mask], np.log(y[mask]), 1)[0]
        lam_meas = float(sl/2.0)   # 閉塞率は振幅 2 次 → 状態の成長率は傾き/2
    lam_pred = c.get("pred_lambda_f")
    rows.append(dict(tag=tag, N=N, method=c["method"],
                     star_vertex=istar, star_frac=star_frac, amp_uniformity=amp_unif,
                     phase_alignment=float(zsum), n_plus=n_plus, n_minus=n_minus,
                     sign_imbalance_pred=float(abs(n_plus - n_minus)/max(len(phs), 1)),
                     z2_alignment=z2_align,
                     mu_final_over_r2=mu_over_r2, star_mu_pred=star_mu,
                     mu_dev=abs(mu_over_r2 - star_mu)/star_mu if star_mu else None,
                     closure_final=sj["closure_final"],
                     lambda_pred=lam_pred, lambda_measured=lam_meas,
                     lambda_ratio=(lam_meas/(lam_pred/2) if (lam_meas and lam_pred and lam_pred > 1e-6) else None),
                     pred_disp1=c.get("pred_disp1"), disp1_measured=sj.get("disp1_measured")))
with open(os.path.join(ROOT, "results", "final_state_structure.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
def fmt(x, d=4):
    return "—" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{d}f}"
def fe(x):
    return "—" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.2e}"
lines = ["# 最終状態の構造（星型アンカー検証）", "",
         "| tag | 星頂点 | 星純度 | 振幅均等(ptp/mean) | 位相整列 | ±内訳 | z²整列 | μ_fin/r̄² | N(N−2)/2 | 相対差 | 閉塞末 | λ実測 | λ予測/2 | 比 | disp1予/実 |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
for r in rows:
    lines.append("| " + " | ".join([r["tag"], str(r["star_vertex"]), fmt(r["star_frac"], 4), fmt(r["amp_uniformity"], 3),
                                    fmt(r["phase_alignment"], 4),
                                    str(r["n_plus"]) + "/" + str(r["n_minus"]), fmt(r["z2_alignment"], 6),
                                    fmt(r["mu_final_over_r2"], 3), fmt(r["star_mu_pred"], 1),
                                    fe(r["mu_dev"]), fe(r["closure_final"]),
                                    fe(r["lambda_measured"]), fe(r["lambda_pred"]/2 if r["lambda_pred"] else None),
                                    fmt(r["lambda_ratio"], 3),
                                    (fe(r["pred_disp1"]) + "/" + fe(r["disp1_measured"]))]) + " |")
with open(os.path.join(ROOT, "results", "final_state_structure.md"), "w") as f:
    f.write("\n".join(lines) + "\n")
print("PASS5 OK（results/final_state_structure.md）")
for r in rows:
    if r["method"] == "hm":
        print(f"{r['tag']}: 星純度={fmt(r['star_frac'],3)} 位相整列={fmt(r['phase_alignment'],3)} "
              f"μ_fin/r̄²={fmt(r['mu_final_over_r2'],3)} vs {r['star_mu_pred']:.1f} (dev={fe(r['mu_dev'])}) "
              f"λ実測={fe(r['lambda_measured'])} λ予測/2={fe(r['lambda_pred']/2 if r['lambda_pred'] else None)} 比={fmt(r['lambda_ratio'],3)}")
