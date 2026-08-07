#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""検証プログラム1番: 三者同時検定（模型内完結形）——辞書仮説 vs 成分仮説

設計: 種の三読出し（時計質量 Mt=⟨δω⟩・線幅 σ_ω・巻き m）を独立に振った
8種を用意し、基準源との重力結合 E を測る。

計測原理（初版の反証から・分解定理の帰結）: 初版は基準源を中性(m=0)とした
ため、中性種(m=0)との対で巻き整合コヒーレント（ゲージ的）チャネルが開き、
重力（盲目）チャネルの質量計測が汚染された（E(B,ref)=−13.3, E(D,ref)=−20.8
の異常深化）。**重力質量の計測は、コヒーレントチャネルを閉じた対＝基準巻きを
全種から到達不能にして行う**。本版は m_ref=5（0,1 から倍加・交差頂点とも
到達不能。(1,5)盲目は G9c で実証済。m=2 からは3次経路のみ=弱汚染を注記）。
盲目（重力）チャネルが読む「質量」が
  辞書仮説   : Mt 単独（三読出しは独立な辞書項目）
  成分仮説   : 不変量 M_inv = √(Mt² + α·σ² + β·m²)（二次形式・(t,R,Q)成分）
のどちらかを、残差比較＋ホールドアウト予測（種Hはフィットに使わない）で判別。

判定（事前固定）:
 (T1a) 成分仮説の優越: rms残差(M2: Mt,σ) < 0.5 × rms残差(辞書: Mt単独)
       かつ α>0。
 (T1b) Q成分: M3(Mt,σ,m) が荷電種 F,G の残差を M2 比で半減、かつ β>0。
 (T1c) ホールドアウト: フィット外の種 H の E を M3 が相対誤差<10%で予測。
使い方: python3 run_t1_threeway_invariant_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_t1", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
base = ex.base

def main():
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    x = np.arange(n); eta = np.arange(ne)
    k = np.rint(np.fft.fftfreq(n, d=1.0 / n)).astype(int)
    L = np.exp(-((np.abs(k) / 3.0) ** 4))
    Wf = ((k % 2) == 0).astype(float) * (1.0 - L); Wb = 1.0 - Wf
    GOLD = 0.6180339887498949
    sea = np.zeros(n, complex)
    for kk in (1, 3, 5, 7, 9, 11):
        sea += (0.2 / np.sqrt(6)) * np.exp(2j * np.pi * kk * x / n
                                           + 2j * np.pi * ((kk * GOLD) % 1.0))

    def step(a2, b2):
        Fa = np.fft.fft(a2, axis=0); Fb = np.fft.fft(b2, axis=0)
        f = (np.sum(np.abs(np.fft.ifft(Fa * Wf[:, None], axis=0)) ** 2, axis=1)
             + np.sum(np.abs(np.fft.ifft(Fb * Wf[:, None], axis=0)) ** 2, axis=1))
        bo = (np.sum(np.abs(np.fft.ifft(Fa * Wb[:, None], axis=0)) ** 2, axis=1)
              + np.sum(np.abs(np.fft.ifft(Fb * Wb[:, None], axis=0)) ** 2, axis=1))
        th = np.arctan2(np.sqrt(f), np.sqrt(bo + 1e-300))
        c, s_ = np.cos(th)[:, None], np.sin(th)[:, None]
        a2, b2 = c * a2 - s_ * b2, s_ * a2 + c * b2
        phi = 2.0 * (np.sin(th) ** 2)[:, None] * np.imag(np.conj(b2) * a2)
        cp, sp_ = np.cos(phi), np.sin(phi)
        return cp * a2 - sp_ * b2, sp_ * a2 + cp * b2

    def ladder_prof(center, fsrc, sig_k):
        prof = np.zeros(n); AMP = 0.05 * np.sqrt(n) * 0.1
        for parity, wgt in (("even", np.sqrt(fsrc)), ("odd", np.sqrt(1 - fsrc))):
            ks = np.arange(4 if parity == "even" else 5, n // 2, 2)
            Wk = np.exp(-0.5 * (ks / sig_k) ** 2)
            sub = np.zeros(n)
            for kk, w in zip(ks, Wk):
                sub += w * np.cos(2 * np.pi * kk * (x - center) / n)
            prof += wgt * sub / np.sqrt(np.sum(sub ** 2))
        return AMP * prof

    def lump(m, fsrc, sig_k, center):
        return (ladder_prof(center, fsrc, sig_k)[:, None]
                * np.exp(2j * np.pi * m * eta / ne)[None, :])

    def run_tau(l2, Tburn=500, Tavg=200):
        a2 = (sea[:, None] * np.ones((1, ne)) + l2).astype(complex)
        b2 = -1j * a2
        acc = np.zeros(n)
        for j in range(Tburn + Tavg):
            ap = a2.copy(); a2, b2 = step(a2, b2)
            if j >= Tburn:
                acc += np.angle(np.einsum("xe,xe->x", np.conj(ap), a2))
        return acc / Tavg

    SPECIES = {  # name: (m, fsrc, sig_k)
        "A": (0, 0.40, 32.0), "B": (0, 0.60, 32.0), "C": (0, 0.80, 32.0),
        "D": (0, 0.60, 24.0), "E": (0, 0.60, 40.0),
        "F": (1, 0.60, 32.0), "G": (2, 0.60, 32.0),
        "H": (1, 0.45, 28.0),  # ホールドアウト（三成分すべて基準から外す）
    }
    REF = (5, 0.60, 32.0)
    SEPS = [24, 52, 81]; cA = 100
    tau0 = run_tau(np.zeros((n, ne), complex))

    # 三読出し（solo・種は cA に置く）
    reads = {}
    solo_at = {}
    for nm, (m, f, sg) in SPECIES.items():
        tau_s = run_tau(lump(m, f, sg, cA))
        solo_at[nm] = tau_s
        prof = ladder_prof(cA, f, sg)
        w = prof ** 2 / np.sum(prof ** 2)
        dω = tau_s - tau0
        Mt = float(np.sum(w * dω))
        sg_ω = float(np.sqrt(max(np.sum(w * dω ** 2) - Mt ** 2, 0.0)))
        reads[nm] = {"Mt": Mt, "sigma": sg_ω, "m": m}
        print(f"種{nm} (m={m}, f={f}, σ_k={sg}): Mt={Mt:+.4e}  σ_ω={sg_ω:.4e}")

    ref_solo = {d: run_tau(lump(*REF, (cA + d) % n)) for d in SEPS}

    def E_with_ref(nm):
        m, f, sg = SPECIES[nm]
        Es = [float(np.sum(run_tau(lump(m, f, sg, cA) + lump(*REF, (cA + d) % n))
                           - solo_at[nm] - ref_solo[d] + tau0)) for d in SEPS]
        return float(np.mean(Es))

    E = {nm: E_with_ref(nm) for nm in SPECIES}
    for nm in SPECIES:
        print(f"E({nm},ref) = {E[nm]:+.4e}")

    FIT = [nm for nm in SPECIES if nm != "H"]
    y = np.array([E[nm] for nm in FIT])

    def rms_fit(Minv):
        # E = c0 + K·Minv を最小二乗（2係数）で解いた残差rms
        A = np.vstack([np.ones_like(Minv), Minv]).T
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        r = y - A @ coef
        return float(np.sqrt(np.mean(r ** 2))), coef

    Mt = np.array([reads[nm]["Mt"] for nm in FIT])
    Sg = np.array([reads[nm]["sigma"] for nm in FIT])
    Mm = np.array([abs(reads[nm]["m"]) for nm in FIT], dtype=float)

    rms_dict, cd = rms_fit(Mt)
    best2 = (1e9, None, None)
    for al in np.linspace(0.0, 30.0, 601):
        r, c = rms_fit(np.sqrt(Mt ** 2 + al * Sg ** 2))
        if r < best2[0]:
            best2 = (r, al, c)
    rms2, alpha2, c2 = best2
    best3 = (1e9, None, None, None)
    sc = np.mean(Mt) ** 2  # β の無次元化: β·(m/1)²·⟨Mt⟩²
    for al in np.linspace(0.0, 30.0, 121):
        for be in np.linspace(0.0, 2.0, 201):
            r, c = rms_fit(np.sqrt(Mt ** 2 + al * Sg ** 2 + be * sc * Mm ** 2))
            if r < best3[0]:
                best3 = (r, al, be, c)
    rms3, alpha3, beta3, c3 = best3

    # 荷電種 F,G の残差（M2 vs M3）
    def resid(nms, al, be, coef):
        out = {}
        for nm in nms:
            Minv = np.sqrt(reads[nm]["Mt"] ** 2 + al * reads[nm]["sigma"] ** 2
                           + be * sc * abs(reads[nm]["m"]) ** 2)
            out[nm] = float(E[nm] - (coef[0] + coef[1] * Minv))
        return out
    rFG2 = resid(["F", "G"], alpha2, 0.0, c2)
    rFG3 = resid(["F", "G"], alpha3, beta3, c3)
    # ホールドアウト H
    H3 = resid(["H"], alpha3, beta3, c3)["H"]
    predH = E["H"] - H3
    errH = abs(H3 / E["H"])

    print(f"\n辞書(Mt単独): rms={rms_dict:.4e}")
    print(f"M2(Mt,σ): rms={rms2:.4e}  α={alpha2:.2f}")
    print(f"M3(Mt,σ,m): rms={rms3:.4e}  α={alpha3:.2f}  β={beta3:.3f}")
    t1a = (rms2 < 0.5 * rms_dict) and (alpha2 > 0)
    fg2 = np.sqrt(np.mean([v ** 2 for v in rFG2.values()]))
    fg3 = np.sqrt(np.mean([v ** 2 for v in rFG3.values()]))
    t1b = (fg3 < 0.5 * fg2) and (beta3 > 0)
    t1c = errH < 0.10
    print(f"(T1a) 成分仮説優越 rms(M2)<0.5·rms(辞書) かつ α>0: "
          f"{rms2:.3e} vs {0.5*rms_dict:.3e} → {'通過' if t1a else '不成立'}")
    print(f"(T1b) Q成分 荷電種残差半減 かつ β>0: F,G残差 M2={fg2:.3e} M3={fg3:.3e}"
          f" β={beta3:.3f} → {'通過' if t1b else '不成立'}")
    print(f"(T1c) ホールドアウトH予測: E={E['H']:+.4e} 予測={predH:+.4e} "
          f"相対誤差={errH:.3f} → {'通過' if t1c else '不成立'}")
    verdict = ("成分仮説成立: 重力チャネルは Mt 単独（辞書）でなく不変量 "
               "√(Mt²+ασ²+βQ²) を読む" if (t1a and t1b and t1c) else "要精査")
    print(verdict)
    out = {"reads": reads, "E": E,
           "rms_dict": rms_dict, "rms_M2": rms2, "alpha_M2": alpha2,
           "rms_M3": rms3, "alpha_M3": alpha3, "beta_M3": beta3,
           "residFG_M2": rFG2, "residFG_M3": rFG3,
           "holdout_H": {"E": E["H"], "pred": predH, "rel_err": errH},
           "T1a": bool(t1a), "T1b": bool(t1b), "T1c": bool(t1c),
           "verdict": verdict, "runtime_sec": time.time() - t0}
    (HERE / "result_t1_threeway_invariant_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()
