#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第7論文 予備実験：第二準安定状態の零閉鎖層内・横安定性（段階1）。

第六論文の力学（原本、SHA固定・不変更）で第二準安定状態を作り、そこから
二軌道（基準／摂動）を同じ発展則で走らせ、既存活性空間 S2 の外側へ向かう
横方向成長率 λ⊥^(2) を測る。

摂動は必ず零閉鎖多様体 M0={Z: Z^†Z=1, Z^TZ=0} 上に保つ：
  摂動方向 η ∈ S2^⊥ ∩ T_Z M0（接空間条件 Z^Tη=0, Re(Z^†η)=0）を
  制約の交差（零空間）として直接構成し、polar retraction で M0 へ戻す。
一般乱数＋複素正規化は使わない（別閉鎖層へ移す誤認を避ける）。

対照：A零摂動（丸め床）、B活性空間内接方向 η∥∈S2、C活性空間外 η⊥（主対象）。
主観測量：位相整列後の全差分 / 固定包絡外差分 / 新規占有量 ΔE_new。
回帰試験：基準軌道の分裂量 f が第六論文 fcurve とビット一致。

まず A0（N=5 実装検証）＋ N=5 A1 のみ実行する。
使い方: python3 run_transverse_stability_v1.py 5
"""

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_n_scaling_lowrank_v1 import (
    LowRankSystem, make_parent, zero_closure_kernel_seed, GAMMA,
)

RESULT_DIR = HERE / "transverse_stability_result_v1"
REF_FCURVE = HERE / "metastable_series_result_v1"

DELTA = 1e-15
GUARD = 500
LEARN = 1000
VALID = 1000
TWIN_STEPS = 2000


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# ---------- 力学（原本ループの忠実複製、power反復σ） ----------
def evolve_step(sys_lr, Z, wp):
    sys_lr.set_theta(np.angle(Z))
    sig, wp = sys_lr.sigma_max_power(wp)
    return sys_lr.cayley_step(Z, sig), wp


def parent_plane_f(Z, p, q):
    Zp = Z - p * (p @ Z) - q * (q @ Z)
    return float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Z) @ Z))


def reconstruct_metastable(n, seed=0, tol=1e-12):
    """第六論文と同一の第二準安定状態を再構成し、窓と摂動開始点を返す。"""
    sys_lr = LowRankSystem(n)
    rng = np.random.default_rng(40260722 + 1000 * n + seed)
    v, residual, sig = make_parent(sys_lr, rng, iters=1200, tol=tol)
    g = zero_closure_kernel_seed(sys_lr, rng)
    Z = v + DELTA * g
    Z = Z / np.linalg.norm(Z)
    p = v.real / np.linalg.norm(v.real)
    q = v.imag - (v.imag @ p) * p
    q = q / np.linalg.norm(q)
    wp = rng.normal(size=sys_lr.m)

    crossing = None
    t_pert = None
    Zs = {}          # t -> Z（交差後を保存）
    fser = []
    t = 0
    # 交差検出→交差後 GUARD+LEARN+VALID まで走行し軌道保存
    while True:
        f = parent_plane_f(Z, p, q)
        fser.append((t, f))
        if crossing is None and f > 0.05:
            crossing = t
            t_pert = crossing + GUARD + LEARN + VALID
        if crossing is not None and t >= crossing:
            Zs[t] = Z.copy()
        if t_pert is not None and t >= t_pert:
            break
        Z, wp = evolve_step(sys_lr, Z, wp)
        t += 1
    return {
        "sys": sys_lr, "v": v, "p": p, "q": q, "parent_residual": residual,
        "crossing": crossing, "t_pert": t_pert, "Zs": Zs, "wp_at_pert": wp.copy(),
        "Z_pert": Zs[t_pert], "fser": fser,
    }


# ---------- 固定包絡 S2（学習/検証再構成誤差） ----------
def build_S2(Zs, t0_learn, t0_valid, n_learn, n_valid, sys_m):
    def frame_matrix(rng_ts):
        cols = []
        for t in rng_ts:
            if t in Zs:
                cols.append(Zs[t].real)
                cols.append(Zs[t].imag)
        return np.column_stack(cols)  # M x 2L
    learn_ts = range(t0_learn, t0_learn + n_learn)
    valid_ts = range(t0_valid, t0_valid + n_valid)
    Xl = frame_matrix(learn_ts)
    U, S, _ = np.linalg.svd(Xl, full_matrices=False)  # U: M x k
    # 検証窓の状態列
    Vcols = frame_matrix(valid_ts)          # M x 2Lv
    tot = np.sum(Vcols ** 2)
    # r ごとの検証残差比
    resid = []
    for r in range(0, min(len(S), sys_m) + 1):
        if r == 0:
            resid.append(1.0)
        else:
            Br = U[:, :r]
            proj = Br @ (Br.T @ Vcols)
            resid.append(float(np.sum((Vcols - proj) ** 2) / tot))
    resid = np.array(resid)
    # 主基準：改善が小さくなり残差が十分小になる最小r
    r_star = len(S)
    for r in range(1, len(S) + 1):
        improve = (resid[r - 1] - resid[r]) / max(resid[r - 1], 1e-300)
        if resid[r] < 1e-8 and improve < 1e-2:
            r_star = r
            break
    # 副基準
    s2 = S ** 2
    cum = np.cumsum(s2) / np.sum(s2)
    r_cum = int(np.searchsorted(cum, 0.999999) + 1)
    r_rel = int(np.sum(S / S[0] > 1e-10))
    return {
        "U": U, "S": S, "r_star": r_star, "r_cum": r_cum, "r_rel": r_rel,
        "valid_resid_curve": resid[:min(len(resid), 40)].tolist(),
        "learn_range": [t0_learn, t0_learn + n_learn],
        "valid_range": [t0_valid, t0_valid + n_valid],
    }


def projector_apply(B, Z):
    """複素 Z を実基底 B（M×r, 直交列）の張る空間へ射影。"""
    return B @ (B.T @ Z.real) + 1j * (B @ (B.T @ Z.imag))


def outside_norm2(B, Z):
    Zo = Z - projector_apply(B, Z)
    return float(np.real(np.conj(Zo) @ Zo))


# ---------- 零閉鎖接空間 ∩ S2直交補 への摂動 ----------
def tangent_perturbation(Z, B, rng, inplane=False):
    """η ∈ (S2 or S2^⊥) ∩ T_Z M0, |η|=1 を制約零空間として構成。

    実表示 η=(ηr,ηi)∈R^{2M}。制約：
      接空間 Re(Z^†η)=0, Re(Z^Tη)=0, Im(Z^Tη)=0
      空間側 B^Tηr=0, B^Tηi=0 （inplane=True なら S2 内＝直交補を反転）
    """
    M = len(Z)
    Zr, Zi = Z.real, Z.imag
    eta_r = rng.normal(size=M)
    eta_i = rng.normal(size=M)

    def proj_space(vr, vi):
        # S2 成分の除去（inplane なら S2 成分のみ残す）
        pr = B @ (B.T @ vr)
        pi = B @ (B.T @ vi)
        if inplane:
            return pr, pi
        return vr - pr, vi - pi

    eta_r, eta_i = proj_space(eta_r, eta_i)

    # 接空間の3勾配 (ηr,ηi)-表示
    gA = (Zr.copy(), Zi.copy())          # Re(Z^†η)
    gB = (Zr.copy(), -Zi.copy())         # Re(Z^Tη)
    gC = (Zi.copy(), Zr.copy())          # Im(Z^Tη)
    grads = [gA, gB, gC]
    # 勾配も同じ空間へ射影してから直交化し、η から除去
    basis = []
    for (gr, gi) in grads:
        gr2, gi2 = proj_space(gr.copy(), gi.copy())
        for (br, bi) in basis:
            d = gr2 @ br + gi2 @ bi
            gr2 = gr2 - d * br
            gi2 = gi2 - d * bi
        nrm = np.sqrt(gr2 @ gr2 + gi2 @ gi2)
        if nrm > 1e-12:
            basis.append((gr2 / nrm, gi2 / nrm))
    for (br, bi) in basis:
        d = eta_r @ br + eta_i @ bi
        eta_r = eta_r - d * br
        eta_i = eta_i - d * bi

    nrm = np.sqrt(eta_r @ eta_r + eta_i @ eta_i)
    eta = (eta_r + 1j * eta_i) / nrm
    return eta


def retract(W):
    """W(≈M0近傍) を polar retraction で M0 上へ。frame [x,y] を直交正規化。"""
    x = np.sqrt(2.0) * W.real
    y = np.sqrt(2.0) * W.imag
    F = np.column_stack([x, y])            # M x 2
    G = F.T @ F                            # 2x2 SPD
    ev, evec = np.linalg.eigh(G)
    Ginvsqrt = evec @ np.diag(1.0 / np.sqrt(ev)) @ evec.T
    Fo = F @ Ginvsqrt
    Z = (Fo[:, 0] + 1j * Fo[:, 1]) / np.sqrt(2.0)
    return Z


def closure_errors(Z):
    return (abs(float(np.real(np.conj(Z) @ Z)) - 1.0),
            abs(complex(Z @ Z)))


# ---------- 観測量 ----------
def phase_align_diff(Zt, Ztil):
    ip = np.vdot(Zt, Ztil)                 # Z^† Ztil
    ph = ip / abs(ip) if abs(ip) > 0 else 1.0
    Zal = Zt * ph                          # e^{iφ*} Z, φ*=arg(Z^†Ztil)
    d = Ztil - Zal
    return float(np.linalg.norm(d)), Zal


def twin_orbit(rec, B, eta, eps, steps, record_every=1):
    sys_lr = rec["sys"]
    p, q = rec["p"], rec["q"]
    Z0 = rec["Z_pert"]
    W = Z0 + eps * eta
    Ztil0 = retract(W)
    ce = closure_errors(Ztil0)
    Zb = Z0.copy()
    Zt = Ztil0.copy()
    wpb = rec["wp_at_pert"].copy()
    wpt = rec["wp_at_pert"].copy()
    rows = []
    for k in range(steps + 1):
        # 観測
        total_pa, Zal = phase_align_diff(Zb, Zt)
        diff = Zt - Zal
        out_diff = np.sqrt(max(0.0, outside_norm2(B, diff)))
        base_out = outside_norm2(B, Zb)
        pert_out = outside_norm2(B, Zt)
        dE_new = pert_out - base_out
        if k % record_every == 0:
            rows.append((k, parent_plane_f(Zb, p, q), parent_plane_f(Zt, p, q),
                         float(np.linalg.norm(Zt - Zb)), total_pa, out_diff,
                         np.sqrt(max(0.0, base_out)), np.sqrt(max(0.0, pert_out)),
                         dE_new,
                         abs(float(np.real(np.conj(Zt) @ Zt)) - 1.0),
                         abs(complex(Zt @ Zt))))
        if k < steps:
            Zb, wpb = evolve_step(sys_lr, Zb, wpb)
            Zt, wpt = evolve_step(sys_lr, Zt, wpt)
    return rows, ce


def fit_lambda(taus, d, floor):
    """log d の線形回帰。数値床100倍以上・R²≥0.99・連続200step以上・最長区間。"""
    d = np.asarray(d)
    taus = np.asarray(taus, float)
    mask = d > max(100 * floor, 1e-300)
    if mask.sum() < 200:
        return None
    idx = np.where(mask)[0]
    # 連続区間のうち最長
    best = None
    start = idx[0]
    prev = idx[0]
    segs = []
    for i in idx[1:]:
        if i == prev + 1:
            prev = i
        else:
            segs.append((start, prev)); start = i; prev = i
    segs.append((start, prev))
    for (a, b) in segs:
        if b - a + 1 < 200:
            continue
        x = taus[a:b + 1]; y = np.log(d[a:b + 1])
        A = np.vstack([x, np.ones_like(x)]).T
        sol, res, *_ = np.linalg.lstsq(A, y, rcond=None)
        slope = sol[0]
        yhat = A @ sol
        ss = np.sum((y - yhat) ** 2); st = np.sum((y - y.mean()) ** 2)
        r2 = 1 - ss / st if st > 0 else 0.0
        length = b - a + 1
        if r2 >= 0.99 and (best is None or length > best["length"]):
            best = {"slope": float(slope), "r2": float(r2),
                    "start": int(a), "end": int(b), "length": int(length)}
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int)
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--eps", type=float, nargs="+", default=[1e-8, 1e-10, 1e-12])
    ap.add_argument("--steps", type=int, default=TWIN_STEPS)
    args = ap.parse_args()
    n = args.n
    RESULT_DIR.mkdir(parents=True, exist_ok=True)

    report = {"n": n, "m": n * (n - 1) // 2,
              "engine_sha256": sha256(HERE / "run_n_scaling_lowrank_v1.py"),
              "params": {"delta": DELTA, "guard": GUARD, "learn": LEARN,
                         "valid": VALID, "twin_steps": args.steps,
                         "seeds": args.seeds, "eps": args.eps,
                         "sigma_method": "sigma_max_power(iters=3)"}}

    # --- 第二準安定状態の再構成 ---
    rec = reconstruct_metastable(n)
    report["crossing"] = rec["crossing"]
    report["t_pert"] = rec["t_pert"]
    report["parent_residual"] = rec["parent_residual"]

    # 回帰試験：f が第六論文 fcurve とビット一致
    ref = REF_FCURVE / f"fcurve_N{n:05d}_delta1e-15_seed0.csv"
    reg = {"reference": str(ref), "checked": False}
    if ref.exists():
        d = {}
        with open(ref) as fh:
            for r in csv.DictReader(fh):
                d[int(r["tau"])] = float(r["f"])
        dev = max((abs(f - d[t]) for (t, f) in rec["fser"] if t in d), default=None)
        reg = {"reference": str(ref), "checked": True,
               "max_f_dev": (None if dev is None else float(dev)),
               "passed": bool(dev is not None and dev < 1e-12)}
    report["regression_vs_paper6"] = reg

    # --- 固定包絡 S2 ---
    t0_learn = rec["crossing"] + GUARD
    t0_valid = rec["crossing"] + GUARD + LEARN
    s2 = build_S2(rec["Zs"], t0_learn, t0_valid, LEARN, VALID, rec["sys"].m)
    B = s2["U"][:, :s2["r_star"]]
    report["S2"] = {k: v for k, v in s2.items() if k not in ("U", "S")}
    report["S2"]["singular_values_head"] = s2["S"][:12].tolist()

    # --- A0: 摂動構成・retraction・零閉鎖検証（seed0, eps=1e-8, C=外側） ---
    rng0 = np.random.default_rng(70260725 + n)
    eta_c = tangent_perturbation(rec["Z_pert"], B, rng0, inplane=False)
    # 接空間・空間制約の残差
    Z0 = rec["Z_pert"]
    tang = {
        "Re_Zdag_eta": float(np.real(np.vdot(Z0, eta_c))),
        "abs_ZT_eta": float(abs(Z0 @ eta_c)),
        "S2perp_leak": float(np.linalg.norm(B.T @ eta_c.real) + np.linalg.norm(B.T @ eta_c.imag)),
        "eta_norm": float(np.linalg.norm(eta_c)),
    }
    Wtest = Z0 + 1e-8 * eta_c
    Ztil_test = retract(Wtest)
    ce_pre = closure_errors(Wtest / np.linalg.norm(Wtest))
    ce_post = closure_errors(Ztil_test)
    report["A0"] = {
        "tangent_constraint_residuals": tang,
        "closure_err_after_naive_normalize": {"norm": ce_pre[0], "ZTZ": ce_pre[1]},
        "closure_err_after_retraction": {"norm": ce_post[0], "ZTZ": ce_post[1]},
        "Z_pert_closure": {"norm": closure_errors(Z0)[0], "ZTZ": closure_errors(Z0)[1]},
    }

    # --- A1: transverse(C) と inplane(B) と zero(A) ---
    fits = []
    ts_saved = []
    floor_hint = None
    for kind, inplane, is_zero in [("transverse", False, False),
                                   ("inplane", True, False),
                                   ("zero", False, True)]:
        for si in range(args.seeds):
            rng = np.random.default_rng(70260725 + 1000 * n + 10 * si + (0 if not inplane else 1))
            eta = tangent_perturbation(Z0, B, rng, inplane=inplane)
            for eps in ([0.0] if is_zero else args.eps):
                rows, ce = twin_orbit(rec, B, eta, eps, args.steps,
                                      record_every=(1 if n <= 20 else 10))
                taus = [r[0] for r in rows]
                dperp = [r[5] for r in rows]           # 位相整列後・固定包絡外差分
                if kind == "zero":
                    floor_hint = float(np.median(dperp[len(dperp) // 2:]) + 1e-300)
                tag = f"N{n:05d}_{kind}_seed{si}_eps{eps:.0e}"
                with open(RESULT_DIR / f"ts_{tag}.csv", "w", newline="") as fh:
                    w = csv.writer(fh)
                    w.writerow(["time", "splitting_fraction", "pert_splitting_fraction",
                                "total_difference", "phasealigned_total", "transverse_difference",
                                "baseline_residual", "perturbed_residual", "dE_new",
                                "norm_error", "conservation_error"])
                    w.writerows(rows)
                fit = fit_lambda(taus, dperp, floor_hint or 1e-15)
                fits.append({"kind": kind, "seed": si, "eps": eps,
                             "retract_closure": {"norm": ce[0], "ZTZ": ce[1]},
                             "dperp_final": float(dperp[-1]),
                             "dperp_start": float(dperp[0]),
                             "lambda_fit": fit})
                if len(ts_saved) < 3 and kind == "transverse":
                    ts_saved.append(tag)
    report["numerical_floor"] = floor_hint
    report["fits"] = fits

    with open(RESULT_DIR / f"summary_N{n:05d}.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, ensure_ascii=False)

    # コンソール要約
    print(json.dumps({
        "n": n, "crossing": rec["crossing"], "t_pert": rec["t_pert"],
        "regression": reg, "S2_dims": {"r_star": s2["r_star"], "r_cum": s2["r_cum"], "r_rel": s2["r_rel"]},
        "A0_closure_after_retraction": report["A0"]["closure_err_after_retraction"],
        "A0_tangent_resid": tang,
        "floor": floor_hint,
        "transverse_lambdas": [f["lambda_fit"] for f in fits if f["kind"] == "transverse"][:3],
    }, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
