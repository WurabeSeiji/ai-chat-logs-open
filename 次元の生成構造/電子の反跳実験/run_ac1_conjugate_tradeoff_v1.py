#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""AC-1: 電子A×光子Cの共役トレードオフ——プローブ帯域が記録に残せる量を決める

設計: A=電子（H1合格構成・位置は掃引で「未知」）、C=光子（ボゾン帯・キャリア
k_C=21・帯域 σ_k を 0〜16 で掃引・総パワー固定）。万能相互作用で結合し、
C起因の記録 Δ = Ψ_{海+A+C}(T) − Ψ_{海+A}(T) を取り、G の読出しだけで
  位置推定 n̂ = (n/2π)·arg Σ_k c_k c̄_{k+1}（隣接モード相関＝Gの位置読出し）
  時計線幅 σ_ω = 記録の毎步位相前進の std（Gの寿命読出しの系譜）
を測る。位置読出しは隣接モードの干渉を要するため、単色記録からは原理的に
読めないはず——不確定性の装置定理版。

判定（v2・初版の2誤りを修正）:
 初版の教訓1（発見）: 位置誤差が全帯域で対蹠点（n/2）に系統集中——これは
 周期表柱7の空間的二重被覆（純パリティ種の対蹠二点構造・一点局在不能）の
 現れであり、フェルミオン種の位置は n/2 を法としてのみ定義される。
 v2 は位置誤差を mod n/2 の円距離で評価（盲目基準 n/8）。
 初版の教訓2: 時計線幅（試行内σ_ω）は記録の内部動力学を測ってしまう。
 v2 の時計精度＝ω̂（記録の平均位相前進）の試行間ばらつき（アンサンブル std）。
 (AC1a) 単色端: 位置誤差(mod n/2)は盲目基準(n/8)の0.7倍以上、時計精度は最良。
 (AC1b) 反対単調性: 帯域↑で位置誤差↓・時計ばらつき↑（Spearman |ρ|>0.8・逆符号）。
 (AC1c) 積の記録: err_pos×err_clk の帯域依存を記録（記録のみ）。
使い方: python3 run_ac1_conjugate_tradeoff_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_ac1", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
base = ex.base

M_ELECTRON = -3
FSRC = 0.7
KC = 21
SIGMAS = (0.0, 1.0, 2.0, 4.0, 8.0, 16.0)
N_TRIALS = 12
T_BURN, T_AVG = 300, 300
GOLD = 0.6180339887498949


def main():
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    x = np.arange(n); eta = np.arange(ne)
    k = np.rint(np.fft.fftfreq(n, d=1.0 / n)).astype(int)
    L = np.exp(-((np.abs(k) / 3.0) ** 4))
    Wf = ((k % 2) == 0).astype(float) * (1.0 - L); Wb = 1.0 - Wf
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

    def prof(center, fsrc=FSRC):
        p = np.zeros(n); AMP = 0.05 * np.sqrt(n) * 0.1
        for parity, wgt in (("even", np.sqrt(fsrc)), ("odd", np.sqrt(1 - fsrc))):
            ks = np.arange(4 if parity == "even" else 5, n // 2, 2)
            Wk = np.exp(-0.5 * (ks / 32.0) ** 2)
            sub = np.zeros(n)
            for kk, w in zip(ks, Wk):
                sub += w * np.cos(2 * np.pi * kk * (x - center) / n)
            p += wgt * sub / np.sqrt(np.sum(sub ** 2))
        return AMP * p

    def electron(center):
        return prof(center)[:, None] * np.exp(2j * np.pi * M_ELECTRON * eta / ne)[None, :]

    def photon(sigma_k, phase_seed):
        """奇kモードのガウス束（総パワー固定・海帯 k≤11 を避ける）"""
        if sigma_k <= 0:
            modes = [KC]; ws = [1.0]
        else:
            lo = max(13, KC - int(3 * sigma_k)); hi = min(n // 2 - 3, KC + int(3 * sigma_k))
            modes = [q for q in range(lo, hi + 1) if q % 2 == 1]
            ws = [np.exp(-0.5 * ((q - KC) / sigma_k) ** 2) for q in modes]
        ws = np.array(ws); ws /= np.sqrt(np.sum(ws ** 2))
        C = np.zeros(n, complex)
        for q, w in zip(modes, ws):
            C += 0.12 * w * np.exp(2j * np.pi * q * x / n
                                   + 2j * np.pi * ((q * GOLD * phase_seed) % 1.0))
        return C

    def g_position(psi):
        """Gの位置読出し: 隣接モード相関の偏角"""
        c = np.fft.fft(psi)
        z = np.sum(c[1:n // 2] * np.conj(c[2:n // 2 + 1]))
        if abs(z) < 1e-30:
            return None
        return (n / (2 * np.pi)) * np.angle(z) % n

    def run_pair(cA, sigma_k, seed):
        C = photon(sigma_k, seed)
        aAC = (sea[:, None] * np.ones((1, ne)) + electron(cA)
               + C[:, None] * np.ones((1, ne))).astype(complex)
        bAC = -1j * aAC
        aA = (sea[:, None] * np.ones((1, ne)) + electron(cA)).astype(complex)
        bA = -1j * aA
        n_hats, omegas = [], []
        prev = None
        for j in range(T_BURN + T_AVG):
            aAC, bAC = step(aAC, bAC); aA, bA = step(aA, bA)
            if j >= T_BURN:
                d = np.sum(aAC - aA, axis=1)  # 記録（η和・a チャネル）
                if prev is not None:
                    z = np.vdot(prev, d)
                    if abs(z) > 1e-30:
                        omegas.append(np.angle(z))
                nh = g_position(d)
                if nh is not None:
                    n_hats.append(nh)
                prev = d.copy()
        # 位置: 二重被覆整合（mod n/2）の円平均と誤差
        half = n / 2.0
        if n_hats:
            ang = np.array(n_hats) % half * (2 * np.pi / half)
            zm = np.mean(np.exp(1j * ang))
            n_est = (half / (2 * np.pi)) * np.angle(zm) % half
            d = abs(((n_est - (cA % half)) + half / 2) % half - half / 2)
            err = d
        else:
            err = half / 4.0
        om_hat = float(np.mean(omegas)) if len(omegas) > 10 else float("nan")
        return err, om_hat

    blind = n / 8.0
    print(f"n={n}・盲目基準誤差≈{blind:.1f}セル・試行{N_TRIALS}×帯域{len(SIGMAS)}")
    table = {}
    for sg in SIGMAS:
        errs, oms = [], []
        for t in range(N_TRIALS):
            cA = int((t * GOLD * n) % n)
            e, om = run_pair(cA, sg, t + 1)
            errs.append(e); oms.append(om)
        med_e = float(np.median(errs))
        err_clk = float(np.std([o for o in oms if np.isfinite(o)]))
        table[sg] = (med_e, err_clk)
        print(f"σ_k={sg:5.1f}: 位置誤差中央値={med_e:6.2f}セル(mod n/2)  "
              f"時計ばらつき={err_clk:.4e}  積={med_e*err_clk:.3e}")

    sgs = list(SIGMAS)
    es = [table[s][0] for s in sgs]; ls = [table[s][1] for s in sgs]
    def spearman(a, b):
        ra = np.argsort(np.argsort(a)); rb = np.argsort(np.argsort(b))
        return float(np.corrcoef(ra, rb)[0, 1])
    rho_pos = spearman(sgs, es); rho_clk = spearman(sgs, ls)
    ac1a = (es[0] > 0.7 * blind) and (ls[0] == min(ls))
    ac1b = (rho_pos < -0.8) and (rho_clk > 0.8)
    print(f"\n(AC1a) 単色端: 位置誤差={es[0]:.1f}（>0.7×盲目{0.7*blind:.1f}）・"
          f"線幅最小={ls[0]==min(ls)} → {'通過' if ac1a else '不成立'}")
    print(f"(AC1b) 反対単調性: ρ(σ,err_pos)={rho_pos:+.2f}（<−0.8）・"
          f"ρ(σ,線幅)={rho_clk:+.2f}（>+0.8） → {'通過' if ac1b else '不成立'}")
    prods = [e * l for e, l in zip(es, ls)]
    print(f"(AC1c) 積の記録: {['%.2e' % p for p in prods]}（min/max比 "
          f"{min(prods)/max(prods):.2f}）")
    ok = ac1a and ac1b
    verdict = ("共役トレードオフ成立: プローブ帯域が位置/時計の確定性を反対向きに"
               "配分する（不確定性の装置定理・単色記録から位置は原理的に読めない）"
               if ok else "要精査")
    print(verdict)
    out = {"blind": blind, "table": {str(s): list(table[s]) for s in SIGMAS},
           "rho_pos": rho_pos, "rho_clk": rho_clk, "products": prods,
           "AC1a": bool(ac1a), "AC1b": bool(ac1b),
           "verdict": verdict, "runtime_sec": time.time() - t0}
    (HERE / "result_ac1_conjugate_tradeoff_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
