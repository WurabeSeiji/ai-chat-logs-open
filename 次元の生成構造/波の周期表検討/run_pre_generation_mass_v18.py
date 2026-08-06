#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v18: 世代=χ帯の質量比（P4）＋離調-質量相関（P5統合）

方法（事前記録）: 同一巻き m=1 の種を χ帯 {(10,12,14),(30,32,34),(50,52,54)} で
構成（世代1,2,3候補）し、海25%中 settle=2000 後に分散補償Gram（v10）で質量²を、
各帯の自前時計レート（復調前の位相前進）で海時計からの離調を測る。
判定: R4: 質量²の帯階層（比がα⁻¹級かは記述）／R5: 質量²と|離調|の正相関
（W/Z予言: t巻き=自前時計離調⇒質量、引き込み実験と整合するか）。
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_v18", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
v1, toy, base = ex.v1, ex.toy, ex.base

S = 8.0; J_SETTLE = 2000; J_WIN = 40

def main():
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    shape = (n, ne)
    ms = np.arange(ne); mm = np.where(ms <= ne // 2, ms, ms - ne)

    def SW(v):
        f = np.fft.fft(v.reshape(shape), axis=0, norm="ortho")
        f[0, :] = 0.0; f[n // 2:, :] = 0.0
        return np.fft.ifft(f, axis=0, norm="ortho").reshape(v.shape)
    def PJ(v, mset):
        f = np.fft.fft(v.reshape(shape), axis=1, norm="ortho")
        keep = np.isin(mm, list(mset)); f[:, ~keep] = 0.0
        return np.fft.ifft(f, axis=1, norm="ortho").reshape(v.shape)

    sea_a = PJ(SW(v1.make_bundle(sp, (21,), "A", scale=1.0)) * S, {1})
    # 海: m=0 に正しく作る（v14種のシフト法・v17の教訓を反映）
    eta = 2 * np.pi * np.arange(ne) / ne
    def SH(v, dm):
        return (v.reshape(shape) * np.exp(1j * dm * eta)[None, :]).reshape(v.shape)
    sea_a = SH(PJ(SW(v1.make_bundle(sp, (37, 39, 41), "A", scale=1.0)) * S, {1}), -1)
    sea_b = SH(PJ(SW(v1.make_bundle(sp, (37, 39, 41), "B", scale=1.0)) * S, {1}), -1)
    pws = float(np.sum(np.abs(sea_a) ** 2) + np.sum(np.abs(sea_b) ** 2))
    scs = np.sqrt(16.0 / pws); sea_a *= scs; sea_b *= scs

    out = {"J_SETTLE": J_SETTLE, "J_WIN": J_WIN, "rows": []}
    print(f"{'帯(世代)':>14} {'質量²(補償)':>12} {'S':>8} {'離調|Δω|':>10}")
    for gname, bins in (("g1(10,12,14)", (10, 12, 14)),
                          ("g2(30,32,34)", (30, 32, 34)),
                          ("g3(50,52,54)", (50, 52, 54))):
        a = PJ(SW(v1.make_bundle(sp, bins, "A", scale=1.0)) * S, {1})
        b = PJ(SW(v1.make_bundle(sp, bins, "B", scale=1.0)) * S, {1})
        pw = float(np.sum(np.abs(a) ** 2) + np.sum(np.abs(b) ** 2))
        sc = np.sqrt(48.0 / pw); a *= sc; b *= sc
        a = a + sea_a; b = b + sea_b
        for _ in range(J_SETTLE):
            a, b, _ = ex.collision_step_exact(a, b, sp)
        A = np.zeros((J_WIN, n, ne), complex); B = np.zeros((J_WIN, n, ne), complex)
        for t in range(J_WIN):
            a, b, _ = ex.collision_step_exact(a, b, sp)
            fa = np.fft.fft(a.reshape(shape), axis=0, norm="ortho")
            fb = np.fft.fft(b.reshape(shape), axis=0, norm="ortho")
            A[t] = np.fft.fft(fa, axis=1, norm="ortho")
            B[t] = np.fft.fft(fb, axis=1, norm="ortho")
        P = np.mean(np.abs(A) ** 2 + np.abs(B) ** 2, axis=0) / 2
        band = (mm[None, :] == 1) & (P > P.max() * 1e-6)
        # 海帯（m=0）の平均時計
        band0 = (mm[None, :] == 0) & (P > P.max() * 1e-6)
        def clockrate(band_):
            idx = np.argwhere(band_)
            ws, wts = [], []
            for (ki, mi) in idx:
                At = A[:, ki, mi]
                ws.append(np.angle(np.sum(At[1:] * np.conj(At[:-1]))))
                wts.append(P[ki, mi])
            wts = np.array(wts) / max(sum(wts), 1e-300)
            return float(np.sum(wts * np.array(ws)))
        w_sp = clockrate(band); w_sea = clockrate(band0)
        det = abs(w_sp - w_sea)
        # 分散補償Gram（v10）
        idx = np.argwhere(band)
        m2s, sms, wts = [], [], []
        for (ki, mi) in idx:
            At = A[:, ki, mi]; Bt = B[:, ki, mi]
            wA = np.angle(np.sum(At[1:] * np.conj(At[:-1])))
            wB = np.angle(np.sum(Bt[1:] * np.conj(Bt[:-1])))
            t_ = np.arange(J_WIN)
            Ad = At * np.exp(-1j * wA * t_); Bd = Bt * np.exp(-1j * wB * t_)
            Gaa = np.mean(np.abs(Ad) ** 2); Gbb = np.mean(np.abs(Bd) ** 2)
            Gab = np.mean(Ad * np.conj(Bd))
            T = 0.5 * (Gaa + Gbb)
            if T <= 0: continue
            m2s.append((Gaa * Gbb - abs(Gab) ** 2) / T ** 2)
            X = Gab.real; Y = -Gab.imag; Z = 0.5 * (Gaa - Gbb)
            sms.append(np.sqrt(X**2+Y**2+Z**2) / T)
            wts.append(P[ki, mi])
        wts = np.array(wts) / max(sum(wts), 1e-300)
        m2 = float(np.sum(wts * np.array(m2s))); sm = float(np.sum(wts * np.array(sms)))
        print(f"{gname:>14} {m2:>12.5f} {sm:>8.4f} {det:>10.5f}")
        out["rows"].append({"gen": gname, "mass2": m2, "S": sm, "detuning": det})
    r = out["rows"]
    if len(r) == 3:
        print(f"\n質量²比 g2/g1={r[1]['mass2']/r[0]['mass2']:.3f}  g3/g2={r[2]['mass2']/r[1]['mass2']:.3f}"
              f"  （α⁻¹=137級か: 記述）")
        dets = [x["detuning"] for x in r]; m2s_ = [x["mass2"] for x in r]
        corr = float(np.corrcoef(dets, m2s_)[0, 1])
        print(f"離調-質量相関 r={corr:+.3f}（P5: 正なら整合）")
        out["ratios"] = [r[1]['mass2']/r[0]['mass2'], r[2]['mass2']/r[1]['mass2']]
        out["detuning_mass_corr"] = corr
    out["runtime_sec"] = time.time() - t0
    (HERE / "pre_generation_mass_result_v18.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s")

if __name__ == "__main__":
    main()
