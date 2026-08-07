#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""S2: σ軸独立ノブの不存在実験（§18-1 の記録・反証2設計）

目的: 静的lump族で σ_ω（R成分候補）を、足場（空間サイズ）と Mt を固定した
まま独立に振れるかを検定する。振れるなら三者検定のσ軸が閉じる。

設計A（無補償）: 偶奇帯域分割 (σ_e,σ_o)=(32−δ,32+δ)。
設計B（二重較正）: 足場をスケール s の静的較正で、Mt を f のセカント法で固定。

判定（事前固定）:
 (S2) 独立ノブ存在 ⟺ 家族内で σ_ω の変動幅>20% を達成しつつ、
      足場と Mt の変動がともに<5%。

結果（実行済・本文§18-1）: 反証——設計Aは足場が+134%まで膨張、設計Bは
f 補償だけで足場が約2倍変動（3ノブ (f,σ_e,σ_o) が (Mt,σ_ω,足場) に本質的に
絡む）。静的lump族に独立ノブは存在しない。σ軸は時間二線源計器（§18-6）待ち。
使い方: python3 run_s2_sigma_knob_refutation_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_s2", UIM / "run_ignition_fate_exact_v3.py")
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


    def prof2(center, fsrc, sig_e, sig_o):
        prof = np.zeros(n); AMP = 0.05 * np.sqrt(n) * 0.1
        for parity, wgt, sg in (("even", np.sqrt(fsrc), sig_e),
                                ("odd", np.sqrt(1 - fsrc), sig_o)):
            ks = np.arange(4 if parity == "even" else 5, n // 2, 2)
            Wk = np.exp(-0.5 * (ks / sg) ** 2)
            sub = np.zeros(n)
            for kk, w in zip(ks, Wk):
                sub += w * np.cos(2 * np.pi * kk * (x - center) / n)
            prof += wgt * sub / np.sqrt(np.sum(sub ** 2))
        return AMP * prof

    def lump2(m, fsrc, sig_e, sig_o, center):
        return (prof2(center, fsrc, sig_e, sig_o)[:, None]
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

    cA = 100
    dx = (x - cA + n // 2) % n - n // 2

    def footprint(fsrc, se, so):
        p = prof2(cA, fsrc, se, so); w = p ** 2 / np.sum(p ** 2)
        return float(np.sqrt(np.sum(w * dx ** 2)))

    tau0 = run_tau(np.zeros((n, ne), complex))

    def dyn_read(fsrc, se, so):
        tau_s = run_tau(lump2(0, fsrc, se, so, cA))
        p = prof2(cA, fsrc, se, so); w = p ** 2 / np.sum(p ** 2)
        dω = tau_s - tau0
        Mt = float(np.sum(w * dω))
        sg = float(np.sqrt(max(np.sum(w * dω ** 2) - Mt ** 2, 0.0)))
        return Mt, sg

    FOOT0 = footprint(0.6, 32.0, 32.0)
    MT0, SG0 = dyn_read(0.6, 32.0, 32.0)
    print(f"基準: foot={FOOT0:.3f} Mt={MT0:+.4e} σ_ω={SG0:.4e}")
    rows = []
    # 設計A: 無補償
    for d in (0, 4, 8, 12):
        Mt, sg = dyn_read(0.6, 32.0 - d, 32.0 + d)
        fp = footprint(0.6, 32.0 - d, 32.0 + d)
        rows.append({"design": "A", "delta": d, "f": 0.6, "Mt": Mt,
                     "sigma": sg, "foot": fp})
        print(f"A δ={d:2d}: Mt={Mt:+.4e}({100*(Mt/MT0-1):+.1f}%) "
              f"σ_ω={sg:.4e}({100*(sg/SG0-1):+.1f}%) foot={fp:.1f}({100*(fp/FOOT0-1):+.1f}%)")
    # 設計B: 足場静的較正 + Mt動的補償（f セカント）
    for dt in (0.0, 0.15, 0.30, 0.45):
        best = (1e9, None)
        for s_ in np.linspace(8.0, 96.0, 353):
            se, so = s_ * (1 - dt), s_ * (1 + dt)
            if se < 6 or so < 6:
                continue
            fp = footprint(0.6, se, so)
            if abs(fp - FOOT0) < best[0]:
                best = (abs(fp - FOOT0), s_)
        s_ = best[1]; se, so = s_ * (1 - dt), s_ * (1 + dt)
        f1, f2 = 0.5, 0.7
        m1, _ = dyn_read(f1, se, so); m2, _ = dyn_read(f2, se, so)
        for _ in range(3):
            if abs(m2 - m1) < 1e-6:
                break
            f3 = min(max(f2 + (MT0 - m2) * (f2 - f1) / (m2 - m1), 0.05), 0.95)
            f1, m1 = f2, m2; f2 = f3; m2, _ = dyn_read(f2, se, so)
        Mt, sg = dyn_read(f2, se, so)
        fp = footprint(f2, se, so)
        rows.append({"design": "B", "delta_t": dt, "f": f2, "s": s_, "Mt": Mt,
                     "sigma": sg, "foot": fp})
        print(f"B δ̃={dt:.2f}: f={f2:.3f} Mt={Mt:+.4e}({100*(Mt/MT0-1):+.1f}%) "
              f"σ_ω={sg:.4e}({100*(sg/SG0-1):+.1f}%) foot={fp:.1f}({100*(fp/FOOT0-1):+.1f}%)")
    # 判定: いずれかの家族で σ変動>20% かつ 足場・Mt変動<5% の点集合(基準含め3点以上)
    ok = False
    for des in ("A", "B"):
        sub = [r for r in rows if r["design"] == des
               and abs(r["foot"] / FOOT0 - 1) < 0.05 and abs(r["Mt"] / MT0 - 1) < 0.05]
        if len(sub) >= 2:
            sgs = [r["sigma"] for r in sub] + [SG0]
            if (max(sgs) - min(sgs)) / SG0 > 0.20:
                ok = True
    verdict = ("独立ノブ存在（σ軸閉鎖可）" if ok else
               "反証: 静的lump族に σ_ω の独立ノブは存在しない（足場・Mt束縛下で"
               "σ変動を作れない）——σ軸は時間二線源計器（§18-6）の開発が前提")
    print(verdict)
    out = {"FOOT0": FOOT0, "MT0": MT0, "SG0": SG0, "rows": rows,
           "S2_independent_knob": bool(ok), "verdict": verdict,
           "runtime_sec": time.time() - t0}
    (HERE / "result_s2_sigma_knob_refutation_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()
