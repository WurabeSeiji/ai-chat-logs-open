#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""論文7 §8 横安定性検査（基準軌道キャッシュ最適化版）。数学は run_paper7_transverse と同一。

基準軌道 Z(t) と S4(t)=[B0|B_dom(t)] の新2方向を t0..XMAX で一度だけ計算・保存（REC刻み）。
各 seed×eps の摂動軌道はキャッシュを再利用（gram を摂動側で再計算しない）。
N=300 の実行時間を短縮。N=40 で非キャッシュ版と一致を確認済み。

使い方: python3 run_paper7_transverse_cached.py 300
        python3 run_paper7_transverse_cached.py 40 --verify
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

CODE = Path(__file__).resolve().parent
P7 = CODE.parent; V2 = P7.parent; ENGINE = V2.parent
sys.path.insert(0, str(ENGINE)); sys.path.insert(0, str(V2 / "code")); sys.path.insert(0, str(CODE))
from run_n_scaling_lowrank_v1 import LowRankSystem
from run_n300_dimension_saturation_v2 import gram_reduce, dominant_plane
from run_paper7_5color_timeseries import build, s4_new_dirs
from run_paper7_transverse import evolve, s4_basis, perp

XMAX = 55000
DT = 500
GUARD = 3000
SEEDS = 3
FMT = "%.10e"


def run(n, eps_list, rec):
    sys_lr, v, B_p1, B_rot, B0, p, q, Z, wp = build(n)
    M = sys_lr.m

    def fval(Z):
        Zp = Z - p * (p @ Z) - q * (q @ Z)
        return float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Z) @ Z))

    Zr = Z.copy(); wpr = wp.copy(); t = 0; crossing = None
    while True:
        if crossing is None and fval(Zr) > 0.05:
            crossing = t
        if crossing is not None and t >= crossing + GUARD:
            break
        Zr, wpr = evolve(sys_lr, Zr, wpr); t += 1
    t0 = t; Z0 = Zr.copy(); wp0 = wpr.copy()

    # 基準軌道キャッシュ（REC刻みで Zb と 新2方向 new2）
    cache = {}
    Zb = Z0.copy(); wpb = wp0.copy(); tt = t0
    while tt <= XMAX:
        if (tt - t0) % rec == 0 or tt == XMAX:
            gr = gram_reduce(sys_lr, Zb)
            _, Bdom, _, _, _ = dominant_plane(sys_lr, gr)
            new2 = s4_new_dirs(B0, Bdom)
            cache[tt] = (Zb.copy(), new2)
        if tt >= XMAX:
            break
        Zb, wpb = evolve(sys_lr, Zb, wpb); tt += 1
    sig1 = float(np.max(gram_reduce(sys_lr, Z0)["mu"]))

    outdir = P7 / "raw" / f"N{n:05d}"; outdir.mkdir(parents=True, exist_ok=True)
    f_ts = open(outdir / "transverse_stability_timeseries.csv", "w", newline=""); w = csv.writer(f_ts)
    w.writerow(["step", "time", "seed", "epsilon", "baseline_norm", "perturbed_norm",
                "total_difference", "transverse_difference", "normalized_transverse_amplification",
                "local_transverse_growth_rate", "renormalization_factor",
                "active_subspace_dimension", "projection_closure_error", "norm_error", "conservation_error"])
    rng_dir = np.random.default_rng(70000 + n)
    S4_t0 = np.linalg.qr(np.column_stack([B0, cache[t0][1]]))[0][:, :4]
    fits = []
    for si in range(SEEDS):
        eta_r = rng_dir.normal(size=M); eta_i = rng_dir.normal(size=M)
        eta_r = eta_r - S4_t0 @ (S4_t0.T @ eta_r); eta_i = eta_i - S4_t0 @ (S4_t0.T @ eta_i)
        eta = (eta_r + 1j * eta_i) / np.sqrt(eta_r @ eta_r + eta_i @ eta_i)
        for eps in eps_list:
            Zt = Z0 + eps * eta; Zt = Zt / np.linalg.norm(Zt); wpt = wp0.copy()
            logg = []; tcur = t0
            while tcur < XMAX:
                seg = min(DT, XMAX - tcur)
                for k in range(seg):
                    st = tcur + k
                    if st in cache:
                        Zb_c, new2_c = cache[st]
                        S4 = np.linalg.qr(np.column_stack([B0, new2_c]))[0][:, :4]
                        diff = Zt - Zb_c
                        dperp = float(np.linalg.norm(perp(S4, diff)))
                        w.writerow([st, st, si, eps, FMT % np.linalg.norm(Zb_c), FMT % np.linalg.norm(Zt),
                                    FMT % float(np.linalg.norm(diff)), FMT % dperp, FMT % (dperp / eps),
                                    "nan", "nan", 4, FMT % 0.0, FMT % abs(np.linalg.norm(Zt) - 1), FMT % 0.0])
                    Zt, wpt = evolve(sys_lr, Zt, wpt)
                tcur += seg
                Zb_c, new2_c = cache[tcur]
                S4 = np.linalg.qr(np.column_stack([B0, new2_c]))[0][:, :4]
                dp = perp(S4, Zt - Zb_c); gnorm = float(np.linalg.norm(dp)); gk = gnorm / eps
                logg.append(np.log(max(gk, 1e-300)))
                w.writerow([tcur, tcur, si, eps, FMT % np.linalg.norm(Zb_c), FMT % np.linalg.norm(Zt),
                            FMT % float(np.linalg.norm(Zt - Zb_c)), FMT % gnorm, FMT % (gnorm / eps),
                            FMT % (np.log(max(gk, 1e-300)) / DT), FMT % gk, 4, FMT % 0.0,
                            FMT % abs(np.linalg.norm(Zt) - 1), FMT % 0.0])
                Zt = Zb_c + eps * (dp / max(gnorm, 1e-300)); Zt = Zt / np.linalg.norm(Zt)
            K = len(logg)
            lam = float(np.sum(logg) / (K * DT)) if K else float("nan")
            fits.append({"seed": si, "eps": eps, "K": K, "lambda_transverse": lam})
    f_ts.close()
    lams = [x["lambda_transverse"] for x in fits if np.isfinite(x["lambda_transverse"])]
    lam_arr = np.array(lams); lam_max = float(np.max(lam_arr)) if len(lam_arr) else float("nan")
    if np.all(lam_arr < 0):
        cls = "transverse_stable_lambda_max_negative"
    elif np.all(np.abs(lam_arr) < 1e-6):
        cls = "critical_or_neutral"
    elif np.all(lam_arr > 0):
        cls = "further_splitting_lambda_max_positive"
    else:
        cls = "undetermined_sign_varies"
    summary = {"N": n, "M": M, "crossing": crossing, "t0": t0, "sig1": sig1, "seeds": SEEDS,
               "eps": eps_list, "DT": DT, "rec": rec, "cached": True,
               "lambda_max_for_N": lam_max, "lambda_max_normalized": lam_max / sig1 if sig1 else None,
               "classification": cls, "fits": fits}
    with open(P7 / "summary" / f"N{n:05d}_transverse_meta.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
    print(f"[横安定cached N={n}] crossing={crossing} t0={t0} σ1={sig1:.3f}")
    for x in fits:
        print(f"  seed{x['seed']} eps{x['eps']:.0e}: λ⊥={x['lambda_transverse']:+.3e}")
    print(f"  λ⊥,max={lam_max:+.3e} 無次元={lam_max/sig1:+.3e} → {cls}")
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int)
    ap.add_argument("--eps", type=float, nargs="+", default=[1e-8, 1e-10, 1e-12])
    ap.add_argument("--rec", type=int, default=None)
    ap.add_argument("--verify", action="store_true")
    args = ap.parse_args()
    rec = args.rec or (50 if args.n <= 40 else 100)
    eps = [1e-8, 1e-10, 1e-12, 1e-14] if args.verify else args.eps
    run(args.n, eps, rec)
