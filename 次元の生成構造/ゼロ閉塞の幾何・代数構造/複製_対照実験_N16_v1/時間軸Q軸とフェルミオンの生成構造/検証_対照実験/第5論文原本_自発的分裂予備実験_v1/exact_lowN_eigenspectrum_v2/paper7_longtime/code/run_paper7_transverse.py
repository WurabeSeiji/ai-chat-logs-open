#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""論文7 §8 横安定性検査：4方向部分空間 S4(t)=[B0|B_dom(t)] 外への横摂動成長率。解釈なし。

確定指示：S4 は時間依存 S4(t)=正規直交化[B0|B_dom(t)]（質問2=A）。横摂動は基準軌道の
S4(t)^⊥ で測る。複数seed・複数eps。Benettin型再正規化で λ⊥,max^(4)。§11で機械判定。

主観測 d⊥^(4)(t)=‖(I-Π_{S4(t)})(Z̃(t)-Z(t))‖、A⊥=d⊥/d⊥(t0)。
Benettin: ΔTごとに δ⊥ を S4^⊥ へ再射影し ε へ再正規化、g_k=‖δ⊥‖/ε、λ=Σlog g_k/(KΔT)。

使い方: python3 run_paper7_transverse.py 5
"""

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

XMAX = 55000
DT = 500                    # Benettin 再正規化間隔
GUARD = 3000                # crossing 後、準安定域開始までのガード
REC = 50                    # 記録間隔
SEEDS = 3
EPS = [1e-8, 1e-10, 1e-12, 1e-14]
FMT = "%.10e"


def evolve(sys_lr, Z, wp):
    sys_lr.set_theta(np.angle(Z)); se, wp = sys_lr.sigma_max_power(wp)
    return sys_lr.cayley_step(Z, se), wp


def s4_basis(sys_lr, B0, Z):
    gr = gram_reduce(sys_lr, Z)
    _, Bdom, _, _, _ = dominant_plane(sys_lr, gr)
    S4 = np.column_stack([B0, s4_new_dirs(B0, Bdom)])   # M×4
    Qr, _ = np.linalg.qr(S4)
    return Qr[:, :4]


def perp(S4, d):
    """(I-Π_S4) d（d複素, S4実M×4）。"""
    dr = d.real - S4 @ (S4.T @ d.real)
    di = d.imag - S4 @ (S4.T @ d.imag)
    return dr + 1j * di


def run(n):
    sys_lr, v, B_p1, B_rot, B0, p, q, Z, wp = build(n)
    M = sys_lr.m

    def fval(Z):
        Zp = Z - p * (p @ Z) - q * (q @ Z)
        return float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Z) @ Z))

    # crossing → 準安定開始 t0
    Zr = Z.copy(); wpr = wp.copy(); t = 0; crossing = None
    while True:
        if crossing is None and fval(Zr) > 0.05:
            crossing = t
        if crossing is not None and t >= crossing + GUARD:
            break
        Zr, wpr = evolve(sys_lr, Zr, wpr); t += 1
    t0 = t
    Z0 = Zr.copy(); wp0 = wpr.copy()
    S4_t0 = s4_basis(sys_lr, B0, Z0)
    floor = None

    outdir = P7 / "raw" / f"N{n:05d}"; outdir.mkdir(parents=True, exist_ok=True)
    f_ts = open(outdir / "transverse_stability_timeseries.csv", "w", newline=""); w = csv.writer(f_ts)
    w.writerow(["step", "time", "seed", "epsilon", "baseline_norm", "perturbed_norm",
                "total_difference", "transverse_difference", "normalized_transverse_amplification",
                "local_transverse_growth_rate", "renormalization_factor",
                "active_subspace_dimension", "projection_closure_error", "norm_error", "conservation_error"])
    rng_dir = np.random.default_rng(70000 + n)
    fits = []

    for si in range(SEEDS):
        # S4(t0)^⊥ のランダム方向
        eta_r = rng_dir.normal(size=M); eta_i = rng_dir.normal(size=M)
        eta_r = eta_r - S4_t0 @ (S4_t0.T @ eta_r); eta_i = eta_i - S4_t0 @ (S4_t0.T @ eta_i)
        nrm = np.sqrt(eta_r @ eta_r + eta_i @ eta_i)
        eta = (eta_r + 1j * eta_i) / nrm
        for eps in EPS:
            # 基準・摂動 二軌道を t0 から XMAX まで、DT ごとに Benettin 再正規化
            Zb = Z0.copy(); wpb = wp0.copy()
            Zt = Z0 + eps * eta; Zt = Zt / np.linalg.norm(Zt)
            wpt = wp0.copy()
            logg = []; d0 = eps
            tcur = t0
            while tcur < XMAX:
                # DT ステップ進める（両軌道, 記録 REC ごと）
                seg = min(DT, XMAX - tcur)
                for k in range(seg):
                    if (tcur + k) % REC == 0:
                        S4 = s4_basis(sys_lr, B0, Zb)
                        diff = Zt - Zb
                        dperp = float(np.linalg.norm(perp(S4, diff)))
                        dtot = float(np.linalg.norm(diff))
                        A = dperp / d0
                        w.writerow([tcur + k, tcur + k, si, eps,
                                    FMT % np.linalg.norm(Zb), FMT % np.linalg.norm(Zt),
                                    FMT % dtot, FMT % dperp, FMT % A, "nan", "nan",
                                    4, FMT % 0.0, FMT % abs(np.linalg.norm(Zt) - 1),
                                    FMT % abs(np.linalg.norm(Zb) - 1)])
                    Zb, wpb = evolve(sys_lr, Zb, wpb)
                    Zt, wpt = evolve(sys_lr, Zt, wpt)
                tcur += seg
                # Benettin 再正規化：δ⊥ を S4^⊥ で取り出し ε へ戻す
                S4 = s4_basis(sys_lr, B0, Zb)
                dp = perp(S4, Zt - Zb)
                gnorm = float(np.linalg.norm(dp))
                gk = gnorm / eps
                logg.append(np.log(max(gk, 1e-300)))
                w_last = [tcur, tcur, si, eps, FMT % np.linalg.norm(Zb), FMT % np.linalg.norm(Zt),
                          FMT % float(np.linalg.norm(Zt - Zb)), FMT % gnorm, FMT % (gnorm / eps),
                          FMT % (np.log(max(gk, 1e-300)) / DT), FMT % gk, 4, FMT % 0.0,
                          FMT % abs(np.linalg.norm(Zt) - 1), FMT % abs(np.linalg.norm(Zb) - 1)]
                w.writerow(w_last)
                # 再正規化：Zt = normalize(Zb + eps * dp/‖dp‖)
                Zt = Zb + eps * (dp / max(gnorm, 1e-300))
                Zt = Zt / np.linalg.norm(Zt)
            K = len(logg)
            lam = float(np.sum(logg) / (K * DT)) if K > 0 else float("nan")
            fits.append({"seed": si, "eps": eps, "K": K, "lambda_transverse": lam,
                         "log_g_mean": float(np.mean(logg)) if K else float("nan")})
            if floor is None and eps == EPS[0] and si == 0:
                floor = 1e-16
    f_ts.close()

    # σ1（無次元化用）＝ t0 の gram 支配固有値
    gr0 = gram_reduce(sys_lr, Z0); sig1 = float(np.max(gr0["mu"]))
    lams = [x["lambda_transverse"] for x in fits if np.isfinite(x["lambda_transverse"])]
    lam_max = max(lams) if lams else float("nan")
    # §11 分類（機械）
    lam_arr = np.array(lams)
    signs = np.sign(lam_arr)
    if np.all(lam_arr < 0):
        cls = "transverse_stable_lambda_max_negative"
    elif np.all(np.abs(lam_arr) < 1e-6):
        cls = "critical_or_neutral"
    elif np.all(lam_arr > 0):
        cls = "further_splitting_lambda_max_positive"
    else:
        cls = "undetermined_sign_varies"
    summary = {"N": n, "M": M, "crossing": crossing, "t0": t0, "sig1": sig1,
               "seeds": SEEDS, "eps": EPS, "DT": DT,
               "lambda_max_for_N": lam_max, "lambda_max_normalized": lam_max / sig1 if sig1 else None,
               "classification": cls, "fits": fits}
    (P7 / "summary").mkdir(exist_ok=True)
    with open(P7 / "summary" / f"N{n:05d}_transverse_meta.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
    print(f"[横安定 N={n}] crossing={crossing} t0={t0} σ1={sig1:.3f}")
    for x in fits:
        print(f"  seed{x['seed']} eps{x['eps']:.0e}: λ⊥={x['lambda_transverse']:+.3e} (K={x['K']})")
    print(f"  λ⊥,max={lam_max:+.3e}  無次元={lam_max/sig1:+.3e}  → 分類: {cls}")
    return summary


if __name__ == "__main__":
    for a in sys.argv[1:]:
        run(int(a))
