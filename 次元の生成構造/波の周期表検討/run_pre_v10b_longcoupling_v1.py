#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v10b 長結合版: 海結合を τ 級まで蓄積させた m番地×質量・スピン列の実測

背景: v10初走行では settle=200 ≪ ウォークτ~3000 のため海結合の効果が数%しか
入らず、海中 m=1..7 の列が4桁一致のままだった。本版は settle∈{2000,4000} で
結合を蓄積させてから分散補償Gramで測る。

測定（事前固定）:
  各 m=1..7: S_m = 純粋m種＋海25%。settle J_S ∈ {2000, 4000}、窓40衝突。
  (a) 質量²(補償)・S・s_z（m帯・分散補償Gram＝v10 Part1と同一計器）
  (b) m帯パワー保持率 = 窓平均P(m帯) / 初期P(m帯) ——寿命列のm依存
      （±1 vs ±5,±7 の実測区別＝素電荷一意性の残る判定点のデータ）
  対照: 海単独（同settle）。

判定（事前固定・記述的）:
  H_m: m列の質量²・S・保持率のいずれかに、settle=4000 で 3桁目以上の
  m依存の広がりが出る → 表のm列が実測で埋まる。
  出なければ「海結合はm列を分化させない（等価性が海中でも近似成立）」を記録。

使い方: python3 run_pre_v10b_longcoupling_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_v10b", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
v1, toy, base = ex.v1, ex.toy, ex.base

S = 8.0; J_WIN = 40; SETTLES = (2000, 4000)

def main() -> None:
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    shape = (n, ne)
    ms = np.arange(ne); mm = np.where(ms <= ne // 2, ms, ms - ne)
    eta = 2 * np.pi * np.arange(ne) / ne

    def single_winding(v):
        f = np.fft.fft(v.reshape(shape), axis=0, norm="ortho")
        f[0, :] = 0.0; f[n // 2:, :] = 0.0
        return np.fft.ifft(f, axis=0, norm="ortho").reshape(v.shape)
    def project_eta(v, m_set):
        f = np.fft.fft(v.reshape(shape), axis=1, norm="ortho")
        keep = np.isin(mm, list(m_set)); f[:, ~keep] = 0.0
        return np.fft.ifft(f, axis=1, norm="ortho").reshape(v.shape)
    def shift_eta(v, dm):
        return (v.reshape(shape) * np.exp(1j * dm * eta)[None, :]).reshape(v.shape)

    a0 = single_winding(v1.make_bundle(sp, (30, 32, 34), "A", scale=1.0)) * S
    b0 = single_winding(v1.make_bundle(sp, (30, 32, 34), "B", scale=1.0)) * S
    a0 = a0 + single_winding(v1.make_bundle(sp, (21,), "A", scale=1.0)) * (0.2 * S)
    pow0 = float(np.sum(np.abs(a0) ** 2) + np.sum(np.abs(b0) ** 2))
    a1 = project_eta(a0, {1}); b1 = project_eta(b0, {1})
    pw = float(np.sum(np.abs(a1) ** 2) + np.sum(np.abs(b1) ** 2))
    sc = np.sqrt(pow0 / pw); a1 *= sc; b1 *= sc
    sea_a = project_eta(v1.make_bundle(sp, (29, 31, 33), "A", scale=1.0) * S, {0})
    sea_b = project_eta(v1.make_bundle(sp, (29, 31, 33), "B", scale=1.0) * S, {0})
    pws = float(np.sum(np.abs(sea_a) ** 2) + np.sum(np.abs(sea_b) ** 2))
    scs = np.sqrt(0.25 * pow0 / pws); sea_a *= scs; sea_b *= scs

    def band_power(a, b, m_t):
        fa = np.fft.fft(np.fft.fft(a.reshape(shape), axis=0, norm="ortho"), axis=1, norm="ortho")
        fb = np.fft.fft(np.fft.fft(b.reshape(shape), axis=0, norm="ortho"), axis=1, norm="ortho")
        P = (np.abs(fa) ** 2 + np.abs(fb) ** 2) / 2
        return float(P[:, mm == m_t].sum())

    def comp_gram_band(A, B, band):
        idx = np.argwhere(band)
        P = np.mean(np.abs(A) ** 2 + np.abs(B) ** 2, axis=0) / 2
        m2s, sms, szs, ws = [], [], [], []
        for (ki, mi) in idx:
            At = A[:, ki, mi]; Bt = B[:, ki, mi]
            wA = np.angle(np.sum(At[1:] * np.conj(At[:-1])))
            wB = np.angle(np.sum(Bt[1:] * np.conj(Bt[:-1])))
            t_ = np.arange(A.shape[0])
            Ad = At * np.exp(-1j * wA * t_); Bd = Bt * np.exp(-1j * wB * t_)
            Gaa = np.mean(np.abs(Ad) ** 2); Gbb = np.mean(np.abs(Bd) ** 2)
            Gab = np.mean(Ad * np.conj(Bd))
            T = 0.5 * (Gaa + Gbb)
            if T <= 0: continue
            det = Gaa * Gbb - abs(Gab) ** 2
            X = Gab.real; Y = -Gab.imag; Z = 0.5 * (Gaa - Gbb)
            m2s.append(det / T ** 2)
            sms.append(np.sqrt(X**2 + Y**2 + Z**2) / T)
            szs.append(Z / T)
            ws.append(P[ki, mi])
        ws = np.array(ws) / max(sum(ws), 1e-300)
        return (float(np.sum(ws * np.array(m2s))), float(np.sum(ws * np.array(sms))),
                float(np.sum(ws * np.array(szs))))

    out = {"J_WIN": J_WIN, "SETTLES": list(SETTLES), "rows": []}
    for js in SETTLES:
        print(f"==== settle={js} ====")
        print(f"{'m':>3} {'質量²(補償)':>12} {'S':>8} {'s_z':>8} {'保持率':>8}")
        for m_t in range(1, 8):
            a = shift_eta(a1, m_t - 1) + sea_a
            b = shift_eta(b1, m_t - 1) + sea_b
            p_init = band_power(a, b, m_t)
            for _ in range(js):
                a, b, _ = ex.collision_step_exact(a, b, sp)
            A = np.zeros((J_WIN, n, ne), complex); B = np.zeros((J_WIN, n, ne), complex)
            for t in range(J_WIN):
                a, b, _ = ex.collision_step_exact(a, b, sp)
                fa = np.fft.fft(a.reshape(shape), axis=0, norm="ortho")
                fb = np.fft.fft(b.reshape(shape), axis=0, norm="ortho")
                A[t] = np.fft.fft(fa, axis=1, norm="ortho")
                B[t] = np.fft.fft(fb, axis=1, norm="ortho")
            P = np.mean(np.abs(A) ** 2 + np.abs(B) ** 2, axis=0) / 2
            band = (mm[None, :] == m_t) & (P > P.max() * 1e-8)
            if band.sum() == 0:
                print(f"{m_t:>+3}  帯空"); continue
            m2, sm, sz = comp_gram_band(A, B, band)
            ret = float(P[:, mm == m_t].sum() / max(p_init, 1e-300))
            print(f"{m_t:>+3} {m2:>12.5f} {sm:>8.4f} {sz:>+8.4f} {ret:>8.4f}")
            out["rows"].append({"settle": js, "m": m_t, "mass2": m2, "S": sm,
                                 "sz": sz, "retention": ret})
        # 海単独対照
        a, b = sea_a.copy(), sea_b.copy()
        for _ in range(js):
            a, b, _ = ex.collision_step_exact(a, b, sp)
        A = np.zeros((J_WIN, n, ne), complex); B = np.zeros((J_WIN, n, ne), complex)
        for t in range(J_WIN):
            a, b, _ = ex.collision_step_exact(a, b, sp)
            fa = np.fft.fft(a.reshape(shape), axis=0, norm="ortho")
            fb = np.fft.fft(b.reshape(shape), axis=0, norm="ortho")
            A[t] = np.fft.fft(fa, axis=1, norm="ortho")
            B[t] = np.fft.fft(fb, axis=1, norm="ortho")
        P = np.mean(np.abs(A) ** 2 + np.abs(B) ** 2, axis=0) / 2
        band = P > P.max() * 1e-8
        m2, sm, sz = comp_gram_band(A, B, band)
        print(f"海単独: 質量²={m2:.5f} S={sm:.4f}")
        out["rows"].append({"settle": js, "m": 0, "mass2": m2, "S": sm,
                             "sz": sz, "retention": None})
    # H_m 判定: settle=4000 の m=1..7 での広がり
    r4 = [r for r in out["rows"] if r["settle"] == 4000 and r["m"] >= 1]
    if len(r4) >= 5:
        m2s = np.array([r["mass2"] for r in r4])
        rets = np.array([r["retention"] for r in r4])
        spread_m2 = float(m2s.max() - m2s.min())
        spread_ret = float(rets.max() - rets.min())
        h_m = bool(spread_m2 > 1e-3 or spread_ret > 1e-2)
        print(f"\nH_m: 質量²広がり={spread_m2:.2e} 保持率広がり={spread_ret:.2e} → {h_m}")
        out["H_m"] = {"spread_mass2": spread_m2, "spread_retention": spread_ret,
                       "pass": h_m}
    out["runtime_sec"] = time.time() - t0
    (HERE / "pre_v10b_longcoupling_result_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s → pre_v10b_longcoupling_result_v1.json")

if __name__ == "__main__":
    main()
