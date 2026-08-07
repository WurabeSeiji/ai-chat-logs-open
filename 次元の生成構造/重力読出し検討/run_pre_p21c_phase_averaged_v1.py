#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""P21c: キャリア位相平均化した二体結合——符号パズルの判別と質量積の再回帰

設計（事前記録）: P21a の符号反転・散らばりはキャリア位相汚染の疑い
（sep=32 は λ=n/3≈170.7 と非整合）。n=512 では奇数Kの整数波長が不在のため、
整列でなく**平均消去**: sep ∈ {24,52,81,109,138,166}（φ_AB=2π·3·sep/512 が
ほぼ均等に 0..2π を掃引）で各ペアの E(sep) を測り、
  Ē = ⟨E⟩_sep（位相平均結合）と A_mod = (max−min)/2（位相変調振幅）
に分解する。判定:
 (i) 符号反転が sep 依存で消えれば（Ē が全ペア同符号）→ P21a の符号は汚染。
     Ē の符号がペア依存で残れば → 本物の構造。
 (ii) log|Ē| の回帰で M_t積 vs 不変量積 vs P奇積 を再判別
     （直交対 S3/S4 の共通相手比を含む）。
使い方: python3 run_pre_p21c_phase_averaged_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from itertools import combinations
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_p21c", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
base = ex.base
K_SEA = 3

def main():
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    x = np.arange(n)
    k = np.rint(np.fft.fftfreq(n, d=1.0 / n)).astype(int)
    L = np.exp(-((np.abs(k) / 3.0) ** 4))
    Wf = ((k % 2) == 0).astype(float) * (1.0 - L); Wb = 1.0 - Wf
    carrier = np.exp(2j * np.pi * K_SEA * x / n)

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

    def ladder_at(fsrc, amp, center, sig_k=32.0):
        prof = np.zeros(n); AMP = amp * np.sqrt(n) * 0.1
        for parity, wgt in (("even", np.sqrt(fsrc)), ("odd", np.sqrt(1 - fsrc))):
            ks = np.arange(4 if parity == "even" else 5, n // 2, 2)
            Wk = np.exp(-0.5 * (ks / sig_k) ** 2)
            sub = np.zeros(n)
            for kk, w in zip(ks, Wk):
                sub += w * np.cos(2 * np.pi * kk * (x - center) / n)
            prof += wgt * sub / np.sqrt(np.sum(sub ** 2))
        return AMP * prof

    def run_tau(lump_prof, Tburn=500, Tavg=200):
        a2 = ((0.2 * carrier + lump_prof)[:, None] * np.ones((1, ne))).astype(complex)
        b2 = -1j * a2
        acc = np.zeros(n)
        for j in range(Tburn + Tavg):
            ap = a2.copy(); a2, b2 = step(a2, b2)
            if j >= Tburn:
                acc += np.angle(np.einsum("xe,xe->x", np.conj(ap), a2))
        return acc / Tavg

    SRCS = {"S1": (0.3, 0.05), "S2": (0.9, 0.05), "S3": (0.3, 0.10),
            "S4": (0.9, 0.10), "S5": (0.6, 0.07)}
    SEPS = [24, 52, 81, 109, 138, 166]
    cA = 100
    tau0 = run_tau(np.zeros(n))
    soloA = {nm: run_tau(ladder_at(f, a, cA)) for nm, (f, a) in SRCS.items()}
    soloB = {}
    for d in SEPS:
        cB = (cA + d) % n
        for nm, (f, a) in SRCS.items():
            soloB[(nm, d)] = run_tau(ladder_at(f, a, cB))
    print(f"solo完了 {time.time()-t0:.0f}s")

    results = {}
    for (nA, nB) in combinations(SRCS.keys(), 2):
        Es = []
        for d in SEPS:
            cB = (cA + d) % n
            t12 = run_tau(ladder_at(*SRCS[nA], cA) + ladder_at(*SRCS[nB], cB))
            E = float(np.sum(t12 - soloA[nA] - soloB[(nB, d)] + tau0))
            Es.append(E)
        Ebar = float(np.mean(Es)); Amod = float((max(Es) - min(Es)) / 2)
        results[f"{nA}+{nB}"] = {"E_by_sep": Es, "E_bar": Ebar, "A_mod": Amod}
        print(f"{nA}+{nB}: Ē={Ebar:+.4e}  変調={Amod:.3e}  符号列=" +
              "".join("+" if e > 0 else "-" for e in Es))

    # 位相平均結合の回帰
    # 単独源の Mt, σ を測る（cA位置・帯内平均）
    dAmask = np.minimum(np.abs(x - cA), n - np.abs(x - cA)) <= 8
    comp = {}
    for nm in SRCS:
        dtt = soloA[nm] - tau0
        comp[nm] = dict(Mt=float(np.mean(np.abs(dtt[dAmask]))),
                        sig=float(np.std(dtt[dAmask])),
                        Podd=SRCS[nm][0] * SRCS[nm][1] ** 2)
    Ebars = np.array([abs(v["E_bar"]) for v in results.values()])
    keys = list(results.keys())
    prods = {
        "M_t積": np.array([comp[k.split("+")[0]]["Mt"] * comp[k.split("+")[1]]["Mt"] for k in keys]),
        "不変量積": np.array([np.sqrt(comp[k.split("+")[0]]["Mt"]**2 + comp[k.split("+")[0]]["sig"]**2)
                             * np.sqrt(comp[k.split("+")[1]]["Mt"]**2 + comp[k.split("+")[1]]["sig"]**2)
                             for k in keys]),
        "P奇積": np.array([comp[k.split("+")[0]]["Podd"] * comp[k.split("+")[1]]["Podd"] for k in keys]),
    }
    print("\n位相平均結合 |Ē| の回帰:")
    for name, P in prods.items():
        ok = (Ebars > 1e-15) & (P > 1e-15)
        pc = np.polyfit(np.log(P[ok]), np.log(Ebars[ok]), 1)
        pred = np.polyval(pc, np.log(P[ok]))
        R2 = 1 - np.sum((np.log(Ebars[ok]) - pred) ** 2) / np.sum((np.log(Ebars[ok]) - np.log(Ebars[ok]).mean()) ** 2)
        print(f"  {name:>8}: 指数p={pc[0]:+.3f} R²={R2:.4f}")
    nneg = sum(1 for v in results.values() if v["E_bar"] < 0)
    print(f"\nĒの符号: 負（引力的）{nneg}/10ペア")
    (HERE / "result_pre_p21c_phase_averaged_v1.json").write_text(
        json.dumps({"results": results, "components": comp}, indent=1, ensure_ascii=False))
    print(f"完了 {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()
