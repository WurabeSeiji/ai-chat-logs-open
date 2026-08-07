#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""P21: どの質量が重力するか——(t,R,Q)成分分解の重力判別

背景（事前記録）: v2 §11.4 の成分仮説（⟨ω⟩=t軸/μ_Gram=R軸/海関係=Q軸）は
重力に等価原理問題を生む。二体時計シフト E(d)∝M₁M₂（P15）の M が M_t=⟨ω⟩
だけなら、質量の担い軸が異なる種の間で重力結合が種依存になる。
閉塞恒等式 x²+y²+z²=t²+R²+Q²（空間側ノルム=質量不変量）が等価原理の機構
候補: 時間ゲージ(τ_t)は M_t に、空間ゲージ(τ_x)は不変量に結合するか。

P21a: 源の組成 f_src と振幅を独立に振った5源・全10ペアで二体結合
  E = Σ[τ12−τ1−τ2+τ0] を測り、log E を log(M_t1·M_t2) / log(Podd1·Podd2) /
  log(A1·A2) に回帰——どの積が最良に潰すか（R²比較）。
P21b: 各源単独で時間ゲージ歪み δτ_t と空間ゲージ歪み δτ_x の帯平均比
  γ_model = ⟨|δτ_x|⟩/⟨|δτ_t|⟩ を測り、組成 f_src への依存を見る。
  閉塞強制の等価原理なら γ は組成に依存しないはず（判定: 変動<30%）。
使い方: python3 run_pre_p21_which_mass_gravitates_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from itertools import combinations
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_p21", UIM / "run_ignition_fate_exact_v3.py")
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

    def run_fields(lump_prof, Tburn=500, Tavg=200):
        a2 = ((0.2 * carrier + lump_prof)[:, None] * np.ones((1, ne))).astype(complex)
        b2 = -1j * a2
        acc_t = np.zeros(n); acc_x = np.zeros(n)
        for j in range(Tburn + Tavg):
            ap = a2.copy(); a2, b2 = step(a2, b2)
            if j >= Tburn:
                acc_t += np.angle(np.einsum("xe,xe->x", np.conj(ap), a2))
                acc_x += np.angle(np.einsum("xe,xe->x", np.conj(a2), np.roll(a2, -1, axis=0)))
        return acc_t / Tavg, acc_x / Tavg

    x_L = n // 2; sep = 32
    cA, cB = x_L - sep // 2, x_L + sep // 2
    dA = np.minimum(np.abs(x - cA), n - np.abs(x - cA)) <= 8
    dxc = np.minimum(np.abs(x - x_L), n - np.abs(x - x_L))
    band = (dxc > 20) & (dxc < 60)

    # 5源: 組成と振幅を独立に振る
    SRCS = {"S1": (0.3, 0.05), "S2": (0.9, 0.05), "S3": (0.3, 0.10),
            "S4": (0.9, 0.10), "S5": (0.6, 0.07)}
    tau0_t, tau0_x = run_fields(np.zeros(n))
    solo = {}
    print("== P21b: 単独源の成分測定と γ_model ==")
    print(f"{'源':>3} {'f':>4} {'amp':>5} | {'M_t=⟨ω⟩':>10} {'σ_ω':>10} {'P奇':>8} | {'γ=|δτx|/|δτt|':>14}")
    for name, (f, a) in SRCS.items():
        prof = ladder_at(f, a, cA)
        tt, tx = run_fields(prof)
        dtt = tt - tau0_t; dtx = tx - tau0_x
        Mt = float(np.mean(np.abs(dtt[dA])))
        sw = float(np.std(dtt[dA]))
        Podd = f * (a ** 2)
        g_t = float(np.mean(np.abs(dtt[band]))); g_x = float(np.mean(np.abs(dtx[band])))
        gamma = g_x / g_t if g_t > 1e-15 else float("nan")
        solo[name] = dict(f=f, amp=a, Mt=Mt, sigma=sw, Podd=Podd, gamma=gamma,
                          tt=tt, prof=prof)
        print(f"{name:>3} {f:>4.1f} {a:>5.2f} | {Mt:>10.3e} {sw:>10.3e} {Podd:>8.4f} | {gamma:>14.3f}")
    gammas = np.array([s["gamma"] for s in solo.values()])
    gcv = float(np.nanstd(gammas) / np.nanmean(gammas))
    print(f"γ_model: 平均={np.nanmean(gammas):.3f} 変動={gcv:.1%} → "
          + ("組成非依存（閉塞強制の等価原理と整合）" if gcv < 0.3 else "組成依存（等価原理は成分選択的）"))

    # P21a: 全10ペアの二体結合
    print("\n== P21a: 二体結合 E と質量積の対応 ==")
    rows = []
    for (nA, nB) in combinations(SRCS.keys(), 2):
        prof12 = ladder_at(*SRCS[nA], cA) + ladder_at(*SRCS[nB], cB)
        t12, _ = run_fields(prof12)
        tB, _ = run_fields(ladder_at(*SRCS[nB], cB))
        E = float(np.sum(t12 - solo[nA]["tt"] - tB + tau0_t))
        rows.append((nA, nB, E))
        print(f"  {nA}+{nB}: E={E:+.4e}")
    # 回帰: log|E| vs 各積
    Es = np.array([abs(r[2]) for r in rows])
    prods = {
        "M_t積": np.array([solo[a]["Mt"] * solo[b]["Mt"] for a, b, _ in rows]),
        "P奇積": np.array([solo[a]["Podd"] * solo[b]["Podd"] for a, b, _ in rows]),
        "振幅積": np.array([solo[a]["amp"] * solo[b]["amp"] for a, b, _ in rows]),
        "不変量積(Mt²+σ²)": np.array([
            np.sqrt(solo[a]["Mt"]**2 + solo[a]["sigma"]**2) *
            np.sqrt(solo[b]["Mt"]**2 + solo[b]["sigma"]**2) for a, b, _ in rows]),
    }
    print(f"\n  回帰 log|E| = c + p·log(積):")
    best = None
    for name, P in prods.items():
        ok = (Es > 1e-15) & (P > 1e-15)
        pc = np.polyfit(np.log(P[ok]), np.log(Es[ok]), 1)
        pred = np.polyval(pc, np.log(P[ok]))
        R2 = 1 - np.sum((np.log(Es[ok]) - pred) ** 2) / np.sum((np.log(Es[ok]) - np.log(Es[ok]).mean()) ** 2)
        print(f"    {name:>16}: 指数p={pc[0]:+.3f} R²={R2:.4f}")
        if best is None or R2 > best[1]:
            best = (name, R2)
    print(f"\n  最良の潰し: {best[0]}（R²={best[1]:.4f}）")
    out = {"solo": {k: {kk: vv for kk, vv in v.items() if kk not in ("tt", "prof")}
                    for k, v in solo.items()},
           "pairs": [{"A": a, "B": b, "E": e} for a, b, e in rows],
           "gamma_cv": gcv, "best_product": best[0], "best_R2": best[1],
           "runtime_sec": time.time() - t0}
    (HERE / "result_pre_p21_which_mass_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {time.time()-t0:.0f}s → result_pre_p21_which_mass_v1.json")

if __name__ == "__main__":
    main()
