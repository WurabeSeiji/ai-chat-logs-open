#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第7論文 予備実験A3：冪的分離の正体（A3-1）／最大有限時間増幅方向（A3-2）／
有限振幅走査による非線形遷移（A3-3）。

A2で半周波帯は無限小ランダム摂動に指数不安定でなく d_H~t^0.74 の冪的分離のみを示した。
A3はこれを分離判定する：
  A3-1 冪分離が漸近中立/過渡増幅/歳差せん断/見かけの冪則のどれか（局所傾き・瞬時率・
       AIC/BICモデル比較・支配平面整列前後）。
  A3-2 帯内の最大有限時間増幅方向 η_opt を有限差分で制限伝播作用素の最大特異ベクトルとして求める。
  A3-3 η_opt に沿った有限振幅走査で非線形遷移（別準安定/第三有限占有モード）を探す。
分類 M0(漸近中立)〜M4(歳差せん断)。正成長を前提に調整しない。

原本エンジン不変更 import、A0/A1/A2機構を再利用。
使い方: python3 run_transient_nonlinear_a3_v1.py 10 --part all
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_transverse_stability_v1 import (
    reconstruct_metastable, evolve_step, retract, closure_errors, parent_plane_f,
)
from run_halfband_stability_a2_v1 import (
    precompute_baseline, plane_decomp, band_bases, basis_from, occ, proj_overlap,
    phase_align, tangent_grads, BAND_LO, BAND_HI, DOM_THRESH,
)

RESULT_DIR = HERE / "transient_nonlinear_a3_result_v1"


# ---------- 共通：摂動軌道の伝播（射影は基準キャッシュを使用） ----------
def propagate(rec, eta, eps, steps, sample_set):
    sys_lr = rec["sys"]
    Zt = retract(rec["Z_pert"] + eps * eta) if eps > 0 else rec["Z_pert"].copy()
    wpt = rec["wp_at_pert"].copy()
    out = {}
    for k in range(steps + 1):
        if k in sample_set:
            out[k] = Zt.copy()
        if k < steps:
            Zt, wpt = evolve_step(sys_lr, Zt, wpt)
    return out


def band_dom_align_diff(Zt, Zb, Bh, Bd):
    """位相整列後の帯内差分 d_H、全差分、支配平面整列後の残差 d_H(平面整列)。"""
    Zta = phase_align(Zt, Zb)
    diff = Zta - Zb
    dH = float(np.sqrt(occ(Bh, diff))) if Bh is not None else 0.0
    dtot = float(np.linalg.norm(diff))
    return dH, dtot, Zta, diff


# ---------- admissible 帯接方向基底（H_0 ∩ T_Z M0, 支配除去, 正規直交） ----------
def admissible_band_basis(Z, Bh0, Bd0):
    grads = []
    for (gr, gi) in tangent_grads(Z):
        gr2, gi2 = gr.copy(), gi.copy()
        if Bd0 is not None:
            gr2 = gr2 - Bd0 @ (Bd0.T @ gr2); gi2 = gi2 - Bd0 @ (Bd0.T @ gi2)
        for (br, bi) in grads:
            d = gr2 @ br + gi2 @ bi
            gr2 = gr2 - d * br; gi2 = gi2 - d * bi
        nrm = np.sqrt(gr2 @ gr2 + gi2 @ gi2)
        if nrm > 1e-12:
            grads.append((gr2 / nrm, gi2 / nrm))

    def to_tangent_band(cr, ci):
        if Bd0 is not None:
            cr = cr - Bd0 @ (Bd0.T @ cr); ci = ci - Bd0 @ (Bd0.T @ ci)
        for (br, bi) in grads:
            d = cr @ br + ci @ bi
            cr = cr - d * br; ci = ci - d * bi
        return cr, ci

    raw = []
    d0 = Bh0.shape[1]
    for i in range(d0):
        for part in (0, 1):
            cr = Bh0[:, i].copy() if part == 0 else np.zeros_like(Bh0[:, i])
            ci = np.zeros_like(Bh0[:, i]) if part == 0 else Bh0[:, i].copy()
            cr, ci = to_tangent_band(cr, ci)
            raw.append((cr, ci))
    # Gram-Schmidt 直交化（2M実表示）
    basis = []
    for (cr, ci) in raw:
        for (br, bi) in basis:
            d = cr @ br + ci @ bi
            cr = cr - d * br; ci = ci - d * bi
        nrm = np.sqrt(cr @ cr + ci @ ci)
        if nrm > 1e-9:
            basis.append((cr / nrm, ci / nrm))
    etas = [(br + 1j * bi) for (br, bi) in basis]
    return etas


# ---------- A3-2：最大有限時間増幅方向 ----------
def a3_2(rec, cache, Ts, eps=1e-9):
    Z0 = rec["Z_pert"]
    Bh0, Bd0 = cache[0]["Bh"], cache[0]["Bd"]
    etas = admissible_band_basis(Z0, Bh0, Bd0)
    K = len(etas)
    maxT = max(Ts)
    sample_set = set(Ts)
    # 各基底方向を伝播（線形域 eps）
    resp = {T: [] for T in Ts}
    for eta in etas:
        traj = propagate(rec, eta, eps, maxT, sample_set)
        for T in Ts:
            Zt = traj[T]; Zb = cache[T]["Zb"]; Bh = cache[T]["Bh"]
            Zta = phase_align(Zt, Zb)
            r = (Zta - Zb) / eps
            # 帯へ制限（B_H(T)）
            aH = np.concatenate([Bh.T @ r.real, Bh.T @ r.imag]) if Bh is not None else np.zeros(1)
            r_out = r.copy()
            if Bh is not None:
                r_out = r_out - (Bh @ (Bh.T @ r.real) + 1j * Bh @ (Bh.T @ r.imag))
            resp[T].append((aH, float(np.linalg.norm(r)), float(np.linalg.norm(r_out))))
    out = {}
    for T in Ts:
        A = np.column_stack([resp[T][k][0] for k in range(K)])  # (2 d_T) x K
        U, S, Vt = np.linalg.svd(A, full_matrices=False)
        # η_opt = Σ V[k,0] η_k
        v1 = Vt[0]
        eta_opt = sum(v1[k] * etas[k] for k in range(K))
        eta_opt = eta_opt / np.linalg.norm(eta_opt)
        v2 = Vt[1] if len(Vt) > 1 else v1
        eta_2 = sum(v2[k] * etas[k] for k in range(K))
        eta_2 = eta_2 / np.linalg.norm(eta_2)
        # η_opt の帯内保持率（応答全体のうち帯内割合）＝ σ1 は帯内応答なので、全応答ノルムと比較
        # 全応答（帯内＋帯外）を再計算
        traj = propagate(rec, eta_opt, eps, T, {T})
        Zt = traj[T]; Zb = cache[T]["Zb"]; Bh = cache[T]["Bh"]
        Zta = phase_align(Zt, Zb); r = (Zta - Zb) / eps
        tot = float(np.linalg.norm(r))
        inb = float(np.sqrt(occ(Bh, r))) if Bh is not None else 0.0
        out[T] = {
            "sigma1": float(S[0]), "sigma2": float(S[1] if len(S) > 1 else 0.0),
            "K": K, "eta_opt": eta_opt, "eta_2": eta_2,
            "opt_total_resp": tot, "opt_inband_resp": inb,
            "opt_band_fraction": (inb / tot) if tot > 0 else 0.0,
        }
    return out, etas


# ---------- A3-1：冪分離の正体 ----------
def _pow(t, a, al, t0): return a * (t + t0) ** al
def _exp(t, a, lam): return a * np.exp(lam * t)
def _lin(t, a, b): return a + b * t
def _logm(t, a, b, t0): return a + b * np.log(t + t0)
def _sat(t, a, t0): return a * t / (t + t0)


def aic_bic(y, yhat, k):
    n = len(y); rss = float(np.sum((y - yhat) ** 2))
    if rss <= 0: rss = 1e-300
    aic = n * np.log(rss / n) + 2 * k
    bic = n * np.log(rss / n) + k * np.log(n)
    return aic, bic, rss


def a3_1(rec, cache, sample_keys, seeds=3, eps=1e-10):
    Z0 = rec["Z_pert"]
    Bh0, Bd0 = cache[0]["Bh"], cache[0]["Bd"]
    from run_halfband_stability_a2_v1 import make_perturbation
    curves = []
    for si in range(seeds):
        rng = np.random.default_rng(80260725 + 10 * si)
        eta = make_perturbation(Z0, "H", Bh0, Bd0, rng)
        maxT = max(sample_keys)
        traj = propagate(rec, eta, eps, maxT, set(sample_keys))
        ts, dH, dH_planealigned, domang = [], [], [], []
        for k in sample_keys:
            if k == 0:
                continue
            Zt = traj[k]; Zb = cache[k]["Zb"]; Bh = cache[k]["Bh"]; Bd = cache[k]["Bd"]
            dh, dtot, Zta, diff = band_dom_align_diff(Zt, Zb, Bh, Bd)
            ts.append(k); dH.append(dh)
            # 支配平面整列：摂動軌道の支配平面を基準の支配平面へ最適整列した後の帯内残差。
            # 実装：diff から支配平面成分(基準Bd)を除去した帯内成分（歳差成分の粗い除去）
            diff2 = diff.copy()
            if Bd is not None:
                diff2 = diff2 - (Bd @ (Bd.T @ diff.real) + 1j * Bd @ (Bd.T @ diff.imag))
            dH_planealigned.append(float(np.sqrt(occ(Bh, diff2))) if Bh is not None else 0.0)
            domang.append(proj_overlap(cache[0]["Bd"], Bd))
        curves.append({"seed": si, "ts": ts, "dH": dH,
                       "dH_planealigned": dH_planealigned, "dom_overlap_vs_t0": domang})

    # 代表 seed0 で解析
    c = curves[0]
    t = np.array(c["ts"], float); d = np.array(c["dH"], float)
    m = (t > 0) & (d > 0)
    t, d = t[m], d[m]
    ld, lt = np.log(d), np.log(t)
    # 局所対数傾き（移動窓）
    w = max(5, len(t) // 15)
    alpha_local = []
    for i in range(len(t)):
        a, b = max(0, i - w), min(len(t), i + w + 1)
        sl = np.polyfit(lt[a:b], ld[a:b], 1)[0]
        alpha_local.append(float(sl))
    # 局所指数率 dlogd/dt
    lam_local = np.gradient(ld, t)
    # モデル比較（AIC/BIC, d 空間）
    models = {}
    fit_range = (t >= t.max() * 0.02)   # 過渡初期を除く
    tt, dd = t[fit_range], d[fit_range]
    try:
        for name, fn, p0, k in [
            ("power", _pow, (dd[0], 0.7, 1.0), 3),
            ("exp", _exp, (dd[0], 1e-4), 2),
            ("linear", _lin, (dd[0], (dd[-1]-dd[0])/(tt[-1]-tt[0])), 2),
            ("log", _logm, (dd[0], 1.0, 1.0), 3),
            ("saturation", _sat, (dd[-1]*2, tt[-1]), 2),
        ]:
            try:
                popt, _ = curve_fit(fn, tt, dd, p0=p0, maxfev=20000)
                yhat = fn(tt, *popt)
                aic, bic, rss = aic_bic(dd, yhat, k)
                models[name] = {"params": [float(x) for x in popt], "aic": float(aic),
                                "bic": float(bic), "rss": float(rss)}
            except Exception as e:
                models[name] = {"error": str(e)}
    except Exception as e:
        models["_error"] = str(e)
    best = min((m for m in models if "aic" in models[m]), key=lambda k_: models[k_]["bic"], default=None)

    # 支配平面整列後に冪が残るか
    d2 = np.array(c["dH_planealigned"], float)[m]
    m2 = d2 > 0
    alpha_pa = float(np.polyfit(lt[m2], np.log(d2[m2]), 1)[0]) if m2.sum() > 5 else None
    ratio_pa = float(np.median(d2[m2] / d[m2])) if m2.sum() > 5 else None

    return {
        "curves": curves,
        "alpha_local_final": float(np.median(alpha_local[-max(3, w):])),
        "alpha_local_series": [float(x) for x in alpha_local],
        "lam_local_final": float(np.median(lam_local[-max(3, w):])),
        "lam_local_early": float(np.median(lam_local[:max(3, w)])),
        "model_comparison": models, "best_model_bic": best,
        "alpha_after_plane_align": alpha_pa,
        "planealign_residual_ratio": ratio_pa,
        "ts": [int(x) for x in t], "dH": [float(x) for x in d],
    }


# ---------- A3-3：有限振幅走査 ----------
def a3_3(rec, cache, eta_opt, eta_2, etas, steps, eps_list, seeds_band=2):
    Z0 = rec["Z_pert"]
    sample_every = 100
    sample_set = set(range(0, steps + 1, sample_every))
    Bh0, Bd0 = cache[0]["Bh"], cache[0]["Bd"]
    rng = np.random.default_rng(90260725)
    eta_rand = etas[rng.integers(len(etas))]
    # 帯外方向・支配平面内方向
    from run_halfband_stability_a2_v1 import make_perturbation
    eta_out = make_perturbation(Z0, "R", Bh0, Bd0, np.random.default_rng(91260725))
    eta_dom = make_perturbation(Z0, "D", Bh0, Bd0, np.random.default_rng(92260725))
    dirs = {"opt": eta_opt, "second": eta_2, "rand_band": eta_rand,
            "out_band": eta_out, "dom": eta_dom}
    results = {}
    for dname, eta in dirs.items():
        rows = []
        for eps in eps_list:
            traj = propagate(rec, eta, eps, steps, sample_set)
            # 実効初期距離（retraction後・位相整列）
            Z0b = cache[0]["Zb"]
            eff0 = float(np.linalg.norm(phase_align(traj[0], Z0b) - Z0b))
            ks = sorted(traj.keys())
            EH = []; Edom = []; Eker = []; dist = []; effr_cols = []
            for k in ks:
                Zt = traj[k]; Zb = cache[k]["Zb"]; Bh = cache[k]["Bh"]; Bd = cache[k]["Bd"]
                Zta = phase_align(Zt, Zb)
                EH.append(occ(Bh, Zta)); Edom.append(occ(Bd, Zta))
                Eker.append(float(1 - occ(Bd, Zta) - occ(Bh, Zta)))
                dist.append(float(np.linalg.norm(Zta - Zb)))
                effr_cols.append(Zta.real); effr_cols.append(Zta.imag)
            ce = closure_errors(traj[ks[-1]])
            X = np.column_stack(effr_cols[-40:] if len(effr_cols) >= 40 else effr_cols)
            s = np.linalg.svd(X, compute_uv=False); lam = s ** 2
            eff_rank = float((lam.sum() ** 2) / np.sum(lam ** 2))
            EHb = [occ(cache[k]["Bh"], cache[k]["Zb"]) for k in ks]
            rows.append({
                "eps": eps, "eff_init_dist": eff0,
                "EH_final": float(EH[-1]), "EH_max": float(np.max(EH)),
                "EH_base_mean": float(np.mean(EHb)), "EH_base_max": float(np.max(EHb)),
                "Edom_final": float(Edom[-1]), "Eker_final": float(Eker[-1]),
                "dist_final": float(dist[-1]), "dist_max": float(np.max(dist)),
                "eff_rank_end": eff_rank,
                "closure_norm": ce[0], "closure_ZTZ": ce[1],
                "EH_series": [float(x) for x in EH][::max(1, len(EH)//20)],
                "dist_series": [float(x) for x in dist][::max(1, len(dist)//20)],
            })
        results[dname] = rows
    return results


def classify_M(a31, a33):
    lam_fin = a31["lam_local_final"]; lam_early = a31["lam_local_early"]
    alpha_pa = a31["alpha_after_plane_align"]
    ratio_pa = a31["planealign_residual_ratio"]
    best = a31["best_model_bic"]
    # 有限振幅遷移の有無：opt方向で EH_max が基準を十分超え、非比例、かつ**有効ランクが増加**
    # （過渡的な EH_max の跳躍だけでは遷移としない。第二有限占有モードには有効ランク>2.5 が必要。
    #   漸近的持続性は run_longtime_relaxation_a3 が確定する。）
    trans = False
    opt = a33["opt"]
    base_max = opt[0]["EH_base_max"]
    for r in opt:
        if (r["EH_max"] > 5 * base_max and r["dist_final"] > 10 * r["eff_init_dist"]
                and r["eff_rank_end"] > 2.5):
            trans = True
    # 分類
    if a31.get("planealign_residual_ratio") is not None and ratio_pa < 0.1:
        M = "M4_precession_shear"
    elif trans:
        M = "M2_finite_amplitude_threshold"
    elif best == "saturation" or (lam_fin < lam_early * 0.5 and lam_fin < 1e-4):
        M = "M0_asymptotic_neutral" if not trans else "M2_finite_amplitude_threshold"
    elif best in ("power", "log") and abs(lam_fin) < 1e-4:
        M = "M1_nonnormal_transient"
    elif best == "power" and lam_fin > 1e-4:
        M = "M3_longtime_power_instability"
    else:
        M = "M1_nonnormal_transient"
    return M, {"lam_local_early": lam_early, "lam_local_final": lam_fin,
               "best_model": best, "alpha_after_plane_align": alpha_pa,
               "planealign_residual_ratio": ratio_pa, "finite_amp_transition": trans}


def run_N(n, parts, long_steps, seeds, a33_steps):
    rec = reconstruct_metastable(n)
    horizon = max(long_steps, a33_steps, 2000)
    se = 50 if n <= 10 else (100 if n <= 20 else 200)
    print(f"[N={n}] 基準軌道キャッシュ（{horizon}step, every {se}） …")
    cache = precompute_baseline(rec, horizon, se)
    keys = sorted(cache.keys())
    report = {"n": n, "m": rec["sys"].m, "t_pert": rec["t_pert"],
              "band_def": {"lo": BAND_LO, "hi": BAND_HI},
              "sigma_dom": cache[0]["smax"],
              "band_dim": int(cache[0]["Bh"].shape[1]) if cache[0]["Bh"] is not None else 0}

    a31 = a32 = a33 = None
    Ts = [100, 500, 1000, 2000]
    if parts in ("all", "a31", "a31a32"):
        long_keys = [k for k in keys if k <= long_steps]
        a31 = a3_1(rec, cache, long_keys, seeds=seeds)
        report["A3_1"] = {k: v for k, v in a31.items() if k != "curves"}
    if parts in ("all", "a32", "a31a32", "a32a33"):
        a32, etas = a3_2(rec, cache, Ts)
        report["A3_2"] = {str(T): {kk: vv for kk, vv in a32[T].items()
                                   if kk not in ("eta_opt", "eta_2")} for T in Ts}
    if parts in ("all", "a33", "a32a33"):
        if a32 is None:
            a32, etas = a3_2(rec, cache, Ts)
        eps_list = [1e-10, 1e-8, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 3e-2, 1e-1]
        a33 = a3_3(rec, cache, a32[max(Ts)]["eta_opt"], a32[max(Ts)]["eta_2"], etas,
                   a33_steps, eps_list)
        report["A3_3"] = a33
    if a31 is not None and a33 is not None:
        M, minfo = classify_M(a31, a33)
        report["M_classification"] = M
        report["M_info"] = minfo

    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    with open(RESULT_DIR / f"a3_summary_N{n:05d}.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, ensure_ascii=False)
    return report, a31, a32, a33


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ns", type=int, nargs="+")
    ap.add_argument("--part", default="all",
                    choices=["all", "a31", "a32", "a33", "a31a32", "a32a33"])
    ap.add_argument("--long_steps", type=int, default=20000)
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--a33_steps", type=int, default=10000)
    args = ap.parse_args()
    for n in args.ns:
        rep, a31, a32, a33 = run_N(n, args.part, args.long_steps, args.seeds, args.a33_steps)
        print(f"\n===== N={n} (M={rep['m']}) 帯dim={rep['band_dim']} σ_dom={rep['sigma_dom']:.2f} =====")
        if a31:
            print("[A3-1 冪分離の正体]")
            print(f"  局所冪α(終)={a31['alpha_local_final']:.3f}  "
                  f"局所指数率λ(早)={a31['lam_local_early']:.2e}→(終)={a31['lam_local_final']:.2e}  "
                  f"(真冪ならλ~α/t→0)")
            mc = a31["model_comparison"]
            order = sorted((k for k in mc if 'bic' in mc[k]), key=lambda k: mc[k]['bic'])
            print("  BIC順:", ", ".join(f"{k}={mc[k]['bic']:.0f}" for k in order))
            print(f"  最良モデル(BIC)={a31['best_model_bic']}  "
                  f"支配平面整列後α={a31['alpha_after_plane_align']}  "
                  f"整列後残差比={a31['planealign_residual_ratio']}")
        if a32:
            print("[A3-2 最大有限時間増幅]")
            for T in [100, 500, 1000, 2000]:
                r = a32[T]
                print(f"  T={T:5d}: σ1={r['sigma1']:.3f} σ2={r['sigma2']:.3f} "
                      f"σ1/σ2={r['sigma1']/max(r['sigma2'],1e-30):.2f} "
                      f"opt帯内保持={r['opt_band_fraction']:.3f}")
        if a33:
            print("[A3-3 有限振幅走査 opt方向]")
            base_max = a33['opt'][0]['EH_base_max']
            print(f"  基準帯占有max={base_max:.2e}")
            print(f"  {'eps':>8} {'実効初期':>9} {'EH_max':>9} {'EH_fin':>9} {'dist_fin':>9} {'effRank':>7} {'閉包誤差':>9}")
            for r in a33['opt']:
                print(f"  {r['eps']:>8.0e} {r['eff_init_dist']:>9.2e} {r['EH_max']:>9.2e} "
                      f"{r['EH_final']:>9.2e} {r['dist_final']:>9.2e} {r['eff_rank_end']:>7.2f} {r['closure_norm']:>9.1e}")
        if "M_classification" in rep:
            print(f"  → 分類: {rep['M_classification']}  {rep['M_info']}")


if __name__ == "__main__":
    main()
