#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第7論文 A3補：長時間緩和の確定。

A3-1〜3で観測された冪的分離と有限振幅閾値の正体を、50000ステップの長時間追跡で確定する。
(1) 基準軌道（無摂動）：第二準安定状態の帯占有・非支配占有が長時間でどうなるか。
(2) 最大増幅方向 η_opt に沿った有限振幅 ε 走査：閾値超えの励起が持続するか緩和するか、
    最終状態が単一モード（有効ランク2）か第二有限占有平面（ランク>2）を持つか。

原本エンジン不変更 import。判定は結果を修正せずそのまま報告。
使い方: python3 run_longtime_relaxation_a3_v1.py 10 40 --steps 50000
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_transverse_stability_v1 import reconstruct_metastable, evolve_step, retract, parent_plane_f
from run_halfband_stability_a2_v1 import (
    precompute_baseline, plane_decomp, band_bases, basis_from, occ, phase_align,
)
from run_transient_nonlinear_a3_v1 import a3_2

RESULT_DIR = HERE / "longtime_relaxation_a3_result_v1"


def plane_occ_spectrum(sys_lr, Z, topk=6):
    planes, smax = plane_decomp(sys_lr, Z)
    occs = []
    for (s, vr, vi) in planes:
        occs.append((float(s), occ(basis_from([(s, vr, vi)]), Z)))
    occs.sort(key=lambda x: -x[1])
    return smax, occs[:topk]


def eff_rank(frames):
    X = np.column_stack([c for Z in frames for c in (Z.real, Z.imag)])
    s = np.linalg.svd(X, compute_uv=False); lam = s ** 2
    return float((lam.sum() ** 2) / np.sum(lam ** 2))


def track(rec, Z0, wp0, steps, sample_every, p, q, late_frac=0.16):
    sys_lr = rec["sys"]
    Z = Z0.copy(); wp = wp0.copy()
    series = []
    late = []
    for k in range(steps + 1):
        if k % sample_every == 0:
            planes, smax = plane_decomp(sys_lr, Z)
            Bd, Bh, sb = band_bases(planes, smax)
            eh = occ(Bh, Z); ed = occ(Bd, Z)
            series.append({"k": k, "f": parent_plane_f(Z, p, q), "EH": eh, "Edom": ed,
                           "Eker": float(1 - eh - ed)})
        if k >= steps * (1 - late_frac) and k % (sample_every) == 0:
            late.append(Z.copy())
        if k < steps:
            Z, wp = evolve_step(sys_lr, Z, wp)
    smax, occs = plane_occ_spectrum(sys_lr, Z)
    er = eff_rank(late)
    top = occs[0][1]; second = occs[1][1] if len(occs) > 1 else 0.0
    return {"series": series, "final_sigma_max": smax,
            "final_plane_spectrum": occs, "eff_rank_late": er,
            "top_occ": top, "second_occ": second,
            "second_over_top": float(second / max(top, 1e-300)),
            "EH_final": series[-1]["EH"], "f_final": series[-1]["f"]}


def decay_rate(series):
    ks = np.array([s["k"] for s in series], float)
    eh = np.array([max(s["EH"], 1e-300) for s in series])
    m = (ks > 0) & (eh > 1e-200)
    if m.sum() < 4:
        return None
    sl = np.polyfit(ks[m], np.log(eh[m]), 1)[0]
    return float(sl)


def run_N(n, steps, eps_list):
    rec = reconstruct_metastable(n)
    sys_lr = rec["sys"]; p, q = rec["p"], rec["q"]
    se = max(200, steps // 250)
    report = {"n": n, "m": sys_lr.m, "t_pert": rec["t_pert"], "steps": steps}

    # (1) 基準軌道
    base = track(rec, rec["Z_pert"], rec["wp_at_pert"], steps, se, p, q)
    report["baseline"] = {
        "EH_decay_rate_per_step": decay_rate(base["series"]),
        "EH_start": base["series"][0]["EH"], "EH_final": base["EH_final"],
        "f_start": base["series"][0]["f"], "f_final": base["f_final"],
        "eff_rank_late": base["eff_rank_late"],
        "final_top_occ": base["top_occ"], "final_second_over_top": base["second_over_top"],
        "EH_series": [(s["k"], s["EH"], s["f"], s["Edom"]) for s in base["series"]][::max(1, len(base["series"])//15)],
    }

    # (2) η_opt に沿った有限振幅
    cache = precompute_baseline(rec, 2000, 50)
    a32, _ = a3_2(rec, cache, [100, 500, 1000, 2000])
    eta_opt = a32[2000]["eta_opt"]
    scans = []
    for eps in eps_list:
        Z0 = retract(rec["Z_pert"] + eps * eta_opt)
        r = track(rec, Z0, rec["wp_at_pert"], steps, se, p, q)
        scans.append({
            "eps": eps,
            "EH_max": float(max(s["EH"] for s in r["series"])),
            "EH_final": r["EH_final"], "eff_rank_late": r["eff_rank_late"],
            "final_top_occ": r["top_occ"], "final_second_over_top": r["second_over_top"],
            "f_final": r["f_final"],
        })
    report["opt_direction_scan"] = scans

    # 分類：第二有限占有平面が漸近的に存在するか
    asymp_second = max(s["final_second_over_top"] for s in scans)
    asymp_rank = max(s["eff_rank_late"] for s in scans)
    if asymp_second > 1e-3 and asymp_rank > 2.5:
        M = "M2or M3_finite_second_mode"
    elif any(s["EH_max"] > 100 * report["baseline"]["EH_start"] for s in scans):
        M = "M1_nonnormal_transient_relaxes_to_single_mode"
    else:
        M = "M0_asymptotic_single_mode"
    report["M_classification"] = M
    report["asymptotic_second_over_top_max"] = float(asymp_second)
    report["asymptotic_eff_rank_max"] = float(asymp_rank)

    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    with open(RESULT_DIR / f"relax_N{n:05d}.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, ensure_ascii=False)
    return report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ns", type=int, nargs="+")
    ap.add_argument("--steps", type=int, default=50000)
    ap.add_argument("--eps", type=float, nargs="+", default=[1e-4, 1e-3, 1e-2, 3e-2, 1e-1])
    args = ap.parse_args()
    for n in args.ns:
        r = run_N(n, args.steps, args.eps)
        b = r["baseline"]
        print(f"\n===== N={n} (M={r['m']}) {args.steps}step 長時間緩和 =====")
        print(f"[基準軌道] 帯占有 EH: {b['EH_start']:.2e}→{b['EH_final']:.2e} "
              f"(減衰率{b['EH_decay_rate_per_step']:.2e}/step)  "
              f"f: {b['f_start']:.3f}→{b['f_final']:.3f}(プラトー)")
        print(f"           最終 支配占有={b['final_top_occ']:.6f} 第二/支配={b['final_second_over_top']:.1e} "
              f"有効ランク後期={b['eff_rank_late']:.3f}")
        print(f"[η_opt有限振幅走査] 最終状態の正体:")
        print(f"  {'eps':>8} {'EH_max':>9} {'EH_final':>9} {'支配占有':>9} {'第二/支配':>9} {'有効ランク':>8}")
        for s in r["opt_direction_scan"]:
            print(f"  {s['eps']:>8.0e} {s['EH_max']:>9.2e} {s['EH_final']:>9.2e} {s['final_top_occ']:>9.6f} "
                  f"{s['final_second_over_top']:>9.1e} {s['eff_rank_late']:>8.3f}")
        print(f"  → 分類: {r['M_classification']}  (漸近第二/支配max={r['asymptotic_second_over_top_max']:.1e}, "
              f"有効ランクmax={r['asymptotic_eff_rank_max']:.3f})")


if __name__ == "__main__":
    main()
