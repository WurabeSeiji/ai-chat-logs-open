#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス3：全走行の判定（予測 vs 実測）と行列表。図なし。
実測分類（機械適用・裁量なし）：
  held      : overlap_deficit_max < 1e-6 かつ closure_max < 1e-10
  inflating : H_perp_frac_max > 0.05
  drifting  : それ以外
照合規則：eq_neutral→held、eq_inflating→inflating、non_equilibrium→（drifting または inflating）。
時計：dphi_mean_first200 / clock_pred_dphi（equilibrium は 1.000 のはず）。
着地：final_residual_new_over_r2 < 1e-6 なら「新アンカーに着地」。μ_final/μ_parent の最良有理近似（分母≤12）。"""
import os, csv, json
from fractions import Fraction
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")

def fe(x):
    return "—" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.1e}"
def ff(x, d=4):
    return "—" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{d}f}"

rows = []
for tag in sorted(os.listdir(DATA)):
    sj = os.path.join(DATA, tag, "summary.json")
    cj = os.path.join(DATA, tag, "parent_checks.json")
    if not (os.path.exists(sj) and os.path.exists(cj)):
        continue
    s = json.load(open(sj)); c = json.load(open(cj))
    if s["overlap_deficit_max"] < 1e-6 and s["closure_max"] < 1e-10:
        measured = "held"
    elif s["Hperp_frac_max"] > 0.05:
        measured = "inflating"
    else:
        measured = "drifting"
    pred = c["pred_kind"]
    match = (pred == "eq_neutral" and measured == "held") or \
            (pred == "eq_inflating" and measured == "inflating") or \
            (pred == "non_equilibrium" and measured in ("drifting", "inflating"))
    ratio = s["final_mu_new"]/c["mu_new"] if abs(c["mu_new"]) > 1e-12 else float("nan")
    fr = Fraction(ratio).limit_denominator(12) if np.isfinite(ratio) else None
    clock = s["dphi_mean_first200"]/s["clock_pred_dphi"] if abs(s["clock_pred_dphi"]) > 1e-12 else float("nan")
    rows.append(dict(tag=tag, N=c["N"], method=c["method"],
                     pred_kind=pred, measured_kind=measured, match=match,
                     res_new_over_r2_parent=c["residual_new_over_r2"],
                     pred_rho_minus_1=c.get("pred_rho_minus_1"), pred_disp1=c.get("pred_disp1"),
                     disp1_measured=s.get("disp1_measured"),
                     unitarity=s["unitarity_max_rel_drift"], closure_max=s["closure_max"],
                     closure_final=s["closure_final"], Hperp_frac_max=s["Hperp_frac_max"],
                     overlap_deficit_max=s["overlap_deficit_max"], PR_over_M_final=s["PR_over_M_final"],
                     clock_ratio_first200=clock,
                     growth_slope=(s["growth_fit"] or {}).get("slope_ln_Hperp_per_step"),
                     growth_R2=(s["growth_fit"] or {}).get("R2"),
                     final_residual_new_over_r2=s["final_residual_new_over_r2"],
                     landed_on_anchor=bool(s["final_residual_new_over_r2"] < 1e-6),
                     mu_ratio_final=ratio,
                     mu_ratio_rational=(str(fr.numerator) + "/" + str(fr.denominator)) if fr else None,
                     mu_ratio_rational_err=(abs(ratio - float(fr)) if fr else None)))
os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
keys = list(rows[0].keys())
with open(os.path.join(ROOT, "results", "matrix_N_by_method.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
nm = sum(1 for r in rows if r["match"])
lines = ["# 干渉保存力学・シード無し系列：予測 vs 実測（機械判定）", "",
         f"予測一致 {nm}/{len(rows)}", "",
         "| tag | 予測 | 実測 | 一致 | res_new/r²(親) | ρ−1 | 閉塞max | H⊥率max | 重なり欠損max | PR/M末 | 時計比 | 着地 | μ比 | 有理近似 |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
for r in rows:
    lines.append("| " + " | ".join([
        r["tag"], r["pred_kind"], r["measured_kind"], "OK" if r["match"] else "NG",
        fe(r["res_new_over_r2_parent"]), fe(r["pred_rho_minus_1"]),
        fe(r["closure_max"]), fe(r["Hperp_frac_max"]), fe(r["overlap_deficit_max"]),
        ff(r["PR_over_M_final"], 3), ff(r["clock_ratio_first200"], 6),
        "YES" if r["landed_on_anchor"] else "no",
        ff(r["mu_ratio_final"], 4), r["mu_ratio_rational"] or "—"]) + " |")
with open(os.path.join(ROOT, "results", "matrix_N_by_method.md"), "w") as f:
    f.write("\n".join(lines) + "\n")
print(f"PASS3 OK: 予測一致 {nm}/{len(rows)}（results/matrix_N_by_method.md）")
