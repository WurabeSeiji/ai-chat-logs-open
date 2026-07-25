#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第7論文 予備実験A2：半周波低占有帯の選択的横不安定性。

A1で確定した構造（支配1枚 σ_dom＝ノルム99.9%、σ≈σ_dom/2 の準縮退帯へ~1e-4）を踏まえ、
半周波帯 H(t)={0.4σ_dom≤σ≤0.6σ_dom の回転部分空間} へ加えた零閉鎖接方向摂動が、
基準軌道の自然な微小占有を超えて選択的に増幅するかを判定する。

仕様（確定指示）の要点：
- 摂動対象は帯（個別平面は追跡不能）。帯射影 Π_H^(0)(t) は基準軌道から一度だけ構成し、
  基準・摂動の両軌道へ適用する（摂動軌道から再定義しない＝成長方向を観測空間に吸収しない）。
- 主観測量は帯内状態差 d_H(t)=|Π_H^(0)(t)(Z̃_t−Z_t^(0))|、位相整列後。d_H~e^{λ_H t}。
- 対照群 H(帯内)／D(支配平面内)／R(帯・支配外接方向)／Z(零摂動)。
- 判定 H1(選択的不安定)〜H5(帯外拡散) を同一コードで。正成長を前提に調整しない。

原本エンジン run_n_scaling_lowrank_v1.py（SHA固定）は不変更で import。
A0/A1機構（tangent_perturbation・retract・reconstruct_metastable）は横安定性版から import。

使い方:
  A2-0: python3 run_halfband_stability_a2_v1.py 10 --seeds 1 --eps 1e-8 --a2_0
  A2-1: python3 run_halfband_stability_a2_v1.py 10 20 40 --seeds 3 --eps 1e-8 1e-10 1e-12
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_n_scaling_lowrank_v1 import LowRankSystem  # noqa: F401  (原本エンジン)
from run_transverse_stability_v1 import (
    reconstruct_metastable, evolve_step, retract, closure_errors,
    parent_plane_f, GUARD, LEARN, VALID, TWIN_STEPS,
)

RESULT_DIR = HERE / "halfband_stability_a2_result_v1"

BAND_LO, BAND_HI = 0.40, 0.60      # 半周波帯 σ/σ_dom 範囲
DOM_THRESH = 0.75                  # σ>DOM*σ_dom を支配帯
SIG_REL = 1e-6


# ---------- 生成子の平面分解（密 eig） ----------
def plane_decomp(sys_lr, Z):
    sys_lr.set_theta(np.angle(Z))
    M = sys_lr.m
    K = np.column_stack([sys_lr.kmatvec(np.eye(M)[:, j]) for j in range(M)])
    w, V = np.linalg.eig(K)
    sig = w.imag
    smax = float(sig.max())
    thr = SIG_REL * smax
    planes = []
    for i in range(M):
        if sig[i] > thr:
            planes.append((float(sig[i]), V[:, i].real.copy(), V[:, i].imag.copy()))
    return planes, smax


def basis_from(planes):
    if not planes:
        return None
    cols = []
    for (_, vr, vi) in planes:
        cols.append(vr); cols.append(vi)
    Q, R = np.linalg.qr(np.column_stack(cols))
    keep = np.abs(np.diag(R)) > 1e-8
    return Q[:, keep]


def band_bases(planes, smax):
    dom = [p for p in planes if p[0] > DOM_THRESH * smax]
    band = [p for p in planes if BAND_LO * smax <= p[0] <= BAND_HI * smax]
    Bd = basis_from(dom)
    Bh = basis_from(band)
    sig_band = [p[0] for p in band]
    return Bd, Bh, sig_band


def occ(B, Z):
    if B is None:
        return 0.0
    return float(np.sum((B.T @ Z.real) ** 2) + np.sum((B.T @ Z.imag) ** 2))


def proj_overlap(Ba, Bb):
    if Ba is None or Bb is None:
        return 0.0
    F = Ba.T @ Bb
    return float(np.sum(F ** 2) / min(Ba.shape[1], Bb.shape[1]))


# ---------- 基準軌道と帯射影のキャッシュ（一度だけ） ----------
def precompute_baseline(rec, steps, sample_every):
    """基準軌道 Z^(0)_t を走らせ、抽出時刻で帯射影 Π_H,Π_dom を構成しキャッシュ。"""
    sys_lr = rec["sys"]
    Zb = rec["Z_pert"].copy()
    wpb = rec["wp_at_pert"].copy()
    cache = {}          # k -> dict(Zb, Bh, Bd, sig_band, smax)
    for k in range(steps + 1):
        if k % sample_every == 0:
            planes, smax = plane_decomp(sys_lr, Zb)
            Bd, Bh, sig_band = band_bases(planes, smax)
            cache[k] = {"Zb": Zb.copy(), "Bh": Bh, "Bd": Bd,
                        "sig_band": sig_band, "smax": smax}
        if k < steps:
            Zb, wpb = evolve_step(sys_lr, Zb, wpb)
    return cache


# ---------- 摂動構成（帯内 H / 支配内 D / 帯支配外 R） ----------
def tangent_grads(Z):
    Zr, Zi = Z.real, Z.imag
    return [(Zr.copy(), Zi.copy()), (Zr.copy(), -Zi.copy()), (Zi.copy(), Zr.copy())]


def make_perturbation(Z, kind, Bh, Bd, rng):
    """η(複素) を構成。kind: H(帯内), D(支配内), R(帯・支配外)。
    接空間 T_Z M0（3勾配除去）＋ kind に応じた空間制約。|η|=1。"""
    M = len(Z)
    if kind == "H":
        assert Bh is not None
        cr = Bh @ rng.normal(size=Bh.shape[1])
        ci = Bh @ rng.normal(size=Bh.shape[1])
        remove = [Bd]                          # 支配平面成分を除去
    elif kind == "D":
        assert Bd is not None
        cr = Bd @ rng.normal(size=Bd.shape[1])
        ci = Bd @ rng.normal(size=Bd.shape[1])
        remove = []
    elif kind == "R":
        cr = rng.normal(size=M); ci = rng.normal(size=M)
        remove = [Bd, Bh]                       # 支配・帯を除去
    else:
        raise ValueError(kind)

    def rm_space(vr, vi):
        for B in remove:
            if B is not None:
                vr = vr - B @ (B.T @ vr)
                vi = vi - B @ (B.T @ vi)
        return vr, vi

    cr, ci = rm_space(cr, ci)
    # 接空間3勾配を同じ除去空間へ射影して直交化、η から除去
    basis = []
    for (gr, gi) in tangent_grads(Z):
        gr2, gi2 = rm_space(gr.copy(), gi.copy())
        for (br, bi) in basis:
            d = gr2 @ br + gi2 @ bi
            gr2 = gr2 - d * br; gi2 = gi2 - d * bi
        nrm = np.sqrt(gr2 @ gr2 + gi2 @ gi2)
        if nrm > 1e-12:
            basis.append((gr2 / nrm, gi2 / nrm))
    for (br, bi) in basis:
        d = cr @ br + ci @ bi
        cr = cr - d * br; ci = ci - d * bi
    nrm = np.sqrt(cr @ cr + ci @ ci)
    return (cr + 1j * ci) / nrm


def phase_align(Zt, Zb):
    ip = np.vdot(Zt, Zb)                        # Zt^† Zb
    ph = ip / abs(ip) if abs(ip) > 0 else 1.0
    return Zt * ph


# ---------- 双子軌道（摂動軌道は基準キャッシュの射影で観測） ----------
def run_pair(rec, cache, sample_keys, kind, eta, eps, steps, sample_every):
    sys_lr = rec["sys"]
    p, q = rec["p"], rec["q"]
    Z0 = rec["Z_pert"]
    Zt = retract(Z0 + eps * eta) if eps > 0 else Z0.copy()
    ce0 = closure_errors(Zt)
    wpt = rec["wp_at_pert"].copy()
    rows = []
    pert_frames = {}
    Bd0 = cache[0]["Bd"]
    for k in range(steps + 1):
        if k % sample_every == 0:
            c = cache[k]
            Zb = c["Zb"]; Bh = c["Bh"]; Bd = c["Bd"]
            Zta = phase_align(Zt, Zb)
            diff = Zta - Zb
            EH0 = occ(Bh, Zb); EHe = occ(Bh, Zta); dEH = EHe - EH0
            dH = float(np.sqrt(occ(Bh, diff)))
            # 帯外・支配外の差分
            dperp = diff.copy()
            if Bh is not None:
                dperp = dperp - (Bh @ (Bh.T @ dperp.real) + 1j * Bh @ (Bh.T @ dperp.imag))
            if Bd is not None:
                dperp = dperp - (Bd @ (Bd.T @ dperp.real) + 1j * Bd @ (Bd.T @ dperp.imag))
            dHbar = float(np.linalg.norm(dperp))
            Edom_b = occ(Bd, Zb); Edom_t = occ(Bd, Zta)
            # 支配平面 主角（基準 t vs t0）＝歳差、および 基準 vs 摂動 支配平面
            prec = proj_overlap(Bd0, Bd)
            ce = closure_errors(Zt)
            smax = c["smax"]
            sb = c["sig_band"]
            rows.append((k,
                         parent_plane_f(Zb, p, q), parent_plane_f(Zta, p, q),
                         EH0, EHe, dEH, dH, dHbar,
                         Edom_b, Edom_t, float(1 - Edom_b - EH0),  # 核占有(近似)
                         smax, (min(sb) if sb else 0.0), (max(sb) if sb else 0.0),
                         (np.mean(sb) if sb else 0.0), (max(sb) / smax if sb else 0.0),
                         prec,
                         float(np.linalg.norm(diff)),
                         ce[0], ce[1]))
            pert_frames[k] = Zta.copy()
        if k < steps:
            Zt, wpt = evolve_step(sys_lr, Zt, wpt)
    return rows, ce0, pert_frames


def eff_rank_window(frames, keys):
    cols = []
    for k in keys:
        Z = frames[k]
        cols.append(Z.real); cols.append(Z.imag)
    X = np.column_stack(cols)
    s = np.linalg.svd(X, compute_uv=False)
    lam = s ** 2
    return float((lam.sum() ** 2) / np.sum(lam ** 2))


def fit_power_exp(taus, d, floor):
    """dH(t) を冪 d~t^α と指数 d~e^{λt} の両方でフィットし判別。marginal(冪α≲1)か
    exponential(指数がR²で勝ち成長係数が大)かを返す。"""
    taus = np.asarray(taus, float); d = np.asarray(d, float)
    m = (taus > 0) & (d > max(floor, 1e-300))
    if m.sum() < 5:
        return None
    x, y = taus[m], d[m]
    lx, ly = np.log(x), np.log(y)
    pa = np.polyfit(lx, ly, 1); pyh = np.polyval(pa, lx)
    pr2 = 1 - np.sum((ly - pyh) ** 2) / max(np.sum((ly - ly.mean()) ** 2), 1e-300)
    ea = np.polyfit(x, ly, 1); eyh = np.polyval(ea, x)
    er2 = 1 - np.sum((ly - eyh) ** 2) / max(np.sum((ly - ly.mean()) ** 2), 1e-300)
    gf = float(d[m][-1] / max(d[m][0], 1e-300))
    return {"alpha": float(pa[0]), "pow_r2": float(pr2),
            "lambda": float(ea[0]), "exp_r2": float(er2),
            "growth_factor": gf, "is_exponential": bool(er2 > pr2 + 0.02 and pa[0] > 1.3)}


def fit_lambda(taus, d, floor, min_steps=200):
    """log d の最長連続区間回帰。床100倍・初期10倍以上・R²≥0.99・min_steps以上。"""
    taus = np.asarray(taus, float); d = np.asarray(d, float)
    d0 = d[0] if d[0] > 0 else floor
    thr = max(100 * floor, 10 * d0, 1e-300)
    mask = d > thr
    if mask.sum() < 3:
        return None
    idx = np.where(mask)[0]
    segs = []; start = prev = idx[0]
    for i in idx[1:]:
        if taus[i] - taus[prev] <= (taus[1] - taus[0]) * 1.5:
            prev = i
        else:
            segs.append((start, prev)); start = i; prev = i
    segs.append((start, prev))
    best = None
    for (a, b) in segs:
        if taus[b] - taus[a] < min_steps or (b - a + 1) < 4:
            continue
        x = taus[a:b + 1]; y = np.log(d[a:b + 1])
        A = np.vstack([x, np.ones_like(x)]).T
        sol, *_ = np.linalg.lstsq(A, y, rcond=None)
        yhat = A @ sol
        ss = np.sum((y - yhat) ** 2); st = np.sum((y - y.mean()) ** 2)
        r2 = 1 - ss / st if st > 0 else 0.0
        length = taus[b] - taus[a]
        if r2 >= 0.99 and (best is None or length > best["span"]):
            best = {"slope": float(sol[0]), "r2": float(r2),
                    "span": float(length), "n": int(b - a + 1)}
    return best


CSV_HEADER = ["time", "f_base", "f_pert", "EH_base", "EH_pert", "dEH", "dH", "dHbar",
              "Edom_base", "Edom_pert", "Eker_base", "sigma_dom", "sig_band_min",
              "sig_band_max", "sig_band_mean", "sigH_over_dom", "dom_precession_overlap",
              "total_diff", "norm_err", "ZTZ_err"]


def classify(lam):
    """H1〜H5 を同一コードで判定。指数不安定は「指数フィットが冪に勝ち α>1.3」かつ
    「対Rで選択的」かつ「d_H優勢」の全てを要求。冪/線形(α≲1)は marginal=H2 とする。"""
    def gf(kind):   # 成長係数の中央値（eps=1e-8のみ、線形域）
        v = [r["pe"]["growth_factor"] for r in lam.get(kind, [])
             if r.get("pe") and r["eps"] == 1e-8]
        return float(np.median(v)) if v else None
    def alpha(kind):
        v = [r["pe"]["alpha"] for r in lam.get(kind, []) if r.get("pe") and r["eps"] == 1e-8]
        return float(np.median(v)) if v else None
    def any_exp(kind):
        return any(r.get("pe") and r["pe"]["is_exponential"] for r in lam.get(kind, []))
    # 帯選択比 dH(H)/dH(R)（eps=1e-8, seed対応平均）
    dHH = np.median([r["dH_final"] for r in lam.get("H", []) if r["eps"] == 1e-8]) if lam.get("H") else None
    dHR = np.median([r["dH_final"] for r in lam.get("R", []) if r["eps"] == 1e-8]) if lam.get("R") else None
    sel = (dHH / dHR) if (dHH and dHR) else None
    ratios = [r["dH_final"] / max(r["dHbar_final"], 1e-300) for r in lam.get("H", []) if r["eps"] == 1e-8]
    dH_dom = (np.median(ratios) > 1.0) if ratios else False
    info = {"alpha_H": alpha("H"), "alpha_R": alpha("R"), "gf_H": gf("H"), "gf_R": gf("R"),
            "band_selectivity_H_over_R": sel, "dH_dominant": bool(dH_dom),
            "H_exponential": any_exp("H")}
    aH = alpha("H")
    if aH is None:
        return "insufficient", info
    # 指数不安定判定：H群が指数的 かつ 対Rで選択的(>2x) かつ d_H優勢
    if any_exp("H") and sel is not None and sel > 2.0 and dH_dom:
        cls = "H1_selective_instability"
    elif not dH_dom and sel is not None and sel < 1.0:
        cls = "H5_band_out_diffusion"
    elif any_exp("H") and any_exp("R"):
        cls = "H4_nonselective_mixing"
    else:
        # 冪/線形(α≲1.3)・成長係数有界 → 中立/marginal
        cls = "H2_marginal_no_selective_instability"
    return cls, info


def run_N(n, seeds, eps_list, steps, sample_every, a2_0=False):
    rec = reconstruct_metastable(n)
    cache = precompute_baseline(rec, steps, sample_every)
    sample_keys = sorted(cache.keys())
    Z0 = rec["Z_pert"]
    Bh0, Bd0 = cache[0]["Bh"], cache[0]["Bd"]
    smax0 = cache[0]["smax"]

    # 基準の帯占有自然変動（床）
    base_EH = np.array([cache[k]["Bh"] is not None and occ(cache[k]["Bh"], cache[k]["Zb"]) or 0.0
                        for k in sample_keys])
    base_dom_prec = np.array([proj_overlap(Bd0, cache[k]["Bd"]) for k in sample_keys])

    report = {
        "n": n, "m": rec["sys"].m, "crossing": rec["crossing"], "t_pert": rec["t_pert"],
        "band_def": {"lo": BAND_LO, "hi": BAND_HI, "dom_thresh": DOM_THRESH},
        "sigma_dom_at_t0": smax0,
        "band_dim_at_t0": int(Bh0.shape[1]) if Bh0 is not None else 0,
        "band_sigma_over_dom_at_t0": (max(cache[0]["sig_band"]) / smax0) if cache[0]["sig_band"] else None,
        "baseline_EH_mean": float(np.mean(base_EH)),
        "baseline_EH_std": float(np.std(base_EH)),
        "baseline_dom_precession_end": float(base_dom_prec[-1]),
        "steps": steps, "sample_every": sample_every,
    }

    groups = ["H", "D", "R"] if a2_0 else ["H", "D", "R"]
    lam = {g: [] for g in groups}
    zero_floor = None

    # Z: 零摂動床
    rows0, _, frames0 = run_pair(rec, cache, sample_keys, "H", None, 0.0, steps, sample_every)
    zero_floor = float(np.median([r[6] for r in rows0][len(rows0)//2:]) + 1e-300)  # dH median tail
    write_csv(n, "Z", 0, 0.0, rows0)
    report["zero_floor_dH"] = zero_floor
    report["baseline_eff_rank"] = eff_rank_window(frames0, [k for k in sample_keys if k >= sample_keys[len(sample_keys)//2]])

    eps_use = [eps_list[0]] if a2_0 else eps_list
    seed_use = 1 if a2_0 else seeds
    for kind in groups:
        for si in range(seed_use):
            rng = np.random.default_rng(80260725 + 1000 * n + 10 * si + {"H": 0, "D": 1, "R": 2}[kind])
            eta = make_perturbation(Z0, kind, Bh0, Bd0, rng)
            # 摂動の帯保持・接空間・支配漏れ
            diag = {
                "in_band_frac": float(occ(Bh0, eta)) if Bh0 is not None else None,
                "dom_leak": float(occ(Bd0, eta)) if Bd0 is not None else None,
                "Re_Zdag_eta": float(np.real(np.vdot(Z0, eta))),
                "abs_ZT_eta": float(abs(Z0 @ eta)),
            }
            for eps in eps_use:
                rows, ce0, frames = run_pair(rec, cache, sample_keys, kind, eta, eps, steps, sample_every)
                write_csv(n, kind, si, eps, rows)
                taus = [r[0] for r in rows]
                dH = [r[6] for r in rows]; dHbar = [r[7] for r in rows]
                er = eff_rank_window(frames, [k for k in sample_keys if k >= sample_keys[len(sample_keys)//2]])
                lam[kind].append({
                    "seed": si, "eps": eps, "pert_diag": diag,
                    "retract_closure": {"norm": ce0[0], "ZTZ": ce0[1]},
                    "dH_final": float(dH[-1]), "dHbar_final": float(dHbar[-1]),
                    "EH_final": float(rows[-1][4]), "dEH_final": float(rows[-1][5]),
                    "eff_rank": er,
                    "pe": fit_power_exp(taus, dH, zero_floor),
                    "lam_dH": fit_lambda(taus, dH, zero_floor),
                    "lam_dHbar": fit_lambda(taus, dHbar, zero_floor),
                })
    report["groups"] = lam
    cls, cinfo = classify(lam)
    report["classification"] = cls
    report["classification_info"] = cinfo
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    with open(RESULT_DIR / f"a2_summary_N{n:05d}.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, ensure_ascii=False)
    return report


def write_csv(n, kind, seed, eps, rows):
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    tag = f"N{n:05d}_{kind}_seed{seed}_eps{eps:.0e}"
    with open(RESULT_DIR / f"a2_ts_{tag}.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(CSV_HEADER); w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ns", type=int, nargs="+")
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--eps", type=float, nargs="+", default=[1e-8, 1e-10, 1e-12])
    ap.add_argument("--steps", type=int, default=TWIN_STEPS)
    ap.add_argument("--a2_0", action="store_true", help="A2-0: N1個・eps1個・seed1・機能確認")
    args = ap.parse_args()
    for n in args.ns:
        se = 5 if n <= 10 else (10 if n <= 20 else 20)
        r = run_N(n, args.seeds, args.eps, args.steps, se, a2_0=args.a2_0)
        print(f"\n=== N={n} (M={r['m']}) t_pert={r['t_pert']} 帯dim={r['band_dim_at_t0']} "
              f"σH/σdom={r['band_sigma_over_dom_at_t0']:.3f} ===")
        print(f"基準帯占有 EH={r['baseline_EH_mean']:.3e}±{r['baseline_EH_std']:.1e}  "
              f"零摂動床 d_H={r['zero_floor_dH']:.2e}  基準有効ランク={r['baseline_eff_rank']:.2f}  "
              f"支配歳差(終)={r['baseline_dom_precession_end']:.3f}")
        for kind in ["H", "D", "R"]:
            for rr in r["groups"].get(kind, []):
                if rr["eps"] != 1e-8:
                    continue
                pe = rr["pe"]
                pes = (f"α={pe['alpha']:+.2f}(R²{pe['pow_r2']:.2f}) λ={pe['lambda']:+.1e}(R²{pe['exp_r2']:.2f}) "
                       f"成長{pe['growth_factor']:.0f}x {'指数' if pe['is_exponential'] else '冪/線形'}") if pe else "なし"
                print(f"  {kind} s{rr['seed']} e1e-8: {pes}  dH={rr['dH_final']:.2e} dHbar={rr['dHbar_final']:.2e} "
                      f"eff={rr['eff_rank']:.2f} 帯保持={rr['pert_diag']['in_band_frac']:.3f}")
        ci = r["classification_info"]
        print(f"  → 分類: {r['classification']}")
        print(f"     α_H={ci['alpha_H']} α_R={ci['alpha_R']} 成長H={ci['gf_H']:.0f}x 成長R={ci['gf_R']:.0f}x "
              f"帯選択比H/R={ci['band_selectivity_H_over_R']:.1f} dH優勢={ci['dH_dominant']} H指数={ci['H_exponential']}")


if __name__ == "__main__":
    main()
