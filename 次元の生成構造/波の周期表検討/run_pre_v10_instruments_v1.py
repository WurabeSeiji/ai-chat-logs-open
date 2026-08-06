#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v10計器: (1)分散補償Gram=番地×質量・スピン列の海中実測
          (2)スピン量子数読出し=SU(2)回帰角判別（θ/2 vs θ 調和）

Part1 規約（事前固定）:
  各モード(k,m)の各チャネルの自前時計 ω_X = arg Σ_t X_{t+1}X̄_t を推定し、
  X̃_t = X_t e^{-iω_X t} と復調してから 2×2 Gram。質量²=detΓ/T²、
  Bloch (X,Y,Z)/T。復調により決定論的分散（逆行波の位相前進）は除去され、
  残る非コヒーレンス＝質量が測れる（空間時間論文の復調Gramと同処方）。
  状態: S_m = 純粋m種＋海25%（m=1..7）。等価性は海が破る→m依存が表の中身。
  対照: 海なし純粋m（等価性→全m同値のはず）・海単独。

Part2 規約（事前固定）:
  チャネル二重項 ψ=(a,b) に U(θ)=exp(-iθσ_y/2) を作用させ、
  F(θ)=Re⟨ψ(0)|ψ(θ)⟩/||ψ||² を θ∈[0,4π]（33点）で測る。
  調和分解: c_half（cos θ/2 成分・4π回帰=スピノル的）と
  c_one（cos θ 成分・2π回帰=ベクトル的）。スピノル重み w=|c_half|/(|c_half|+|c_one|)。
  フェルミオン帯（χ偶|k|≥4）とボゾン帯（χ奇）で別々に測り比較。

使い方: python3 run_pre_v10_instruments_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_v10", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
v1, toy, base = ex.v1, ex.toy, ex.base

S = 8.0; J_SETTLE = 160; J_WIN = 40

def main() -> None:
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    shape = (n, ne)
    ms = np.arange(ne); mm = np.where(ms <= ne // 2, ms, ms - ne)
    ks = np.arange(n); kk = np.where(ks <= n // 2, ks, ks - n)
    ferm_k = (np.abs(kk) % 2 == 0) & (np.abs(kk) >= 4)
    bos_k = (np.abs(kk) % 2 == 1)
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

    def capture(a, b):
        for _ in range(J_SETTLE):
            a, b, _ = ex.collision_step_exact(a, b, sp)
        A = np.zeros((J_WIN, n, ne), complex); B = np.zeros((J_WIN, n, ne), complex)
        for t in range(J_WIN):
            a, b, _ = ex.collision_step_exact(a, b, sp)
            fa = np.fft.fft(a.reshape(shape), axis=0, norm="ortho")
            fb = np.fft.fft(b.reshape(shape), axis=0, norm="ortho")
            A[t] = np.fft.fft(fa, axis=1, norm="ortho")
            B[t] = np.fft.fft(fb, axis=1, norm="ortho")
        return A, B, a, b

    def comp_gram(A, B, band):
        """分散補償Gram: 各モード・各チャネルを自前時計で復調してからGram。"""
        idx = np.argwhere(band)
        P = np.mean(np.abs(A) ** 2 + np.abs(B) ** 2, axis=0) / 2
        w_tot = float(P[band].sum())
        m2s, sms, szs, ws = [], [], [], []
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
            det = Gaa * Gbb - abs(Gab) ** 2
            X = Gab.real; Y = -Gab.imag; Z = 0.5 * (Gaa - Gbb)
            m2s.append(det / T ** 2)
            sms.append(np.sqrt(X**2 + Y**2 + Z**2) / T)
            szs.append(Z / T)
            ws.append(P[ki, mi])
        ws = np.array(ws) / max(sum(ws), 1e-300)
        return (float(np.sum(ws * np.array(m2s))), float(np.sum(ws * np.array(sms))),
                float(np.sum(ws * np.array(szs))), w_tot)

    out = {"J_SETTLE": J_SETTLE, "J_WIN": J_WIN, "part1": [], "part2": []}
    print("== Part1: 分散補償Gram（海中の番地×質量・スピン） ==")
    print(f"{'状態':>14} {'帯':>4} {'質量²(補償)':>11} {'S(補償)':>9} {'s_z':>7}")
    # 海単独
    A, B, _, _ = capture(sea_a.copy(), sea_b.copy())
    P = np.mean(np.abs(A) ** 2 + np.abs(B) ** 2, axis=0) / 2
    band = P > P.max() * 1e-6
    m2, sm, sz, _ = comp_gram(A, B, band)
    print(f"{'海単独':>14} {'全':>4} {m2:>11.3e} {sm:>9.3f} {sz:>+7.3f}")
    out["part1"].append({"state": "sea", "band": "all", "mass2": m2, "S": sm, "sz": sz})
    # 海なし対照（等価性チェック: m=1と m=4 が同値のはず）
    for m_t in (1, 4):
        a = shift_eta(a1, m_t - 1); b = shift_eta(b1, m_t - 1)
        A, B, _, _ = capture(a, b)
        P = np.mean(np.abs(A) ** 2 + np.abs(B) ** 2, axis=0) / 2
        band = (mm[None, :] == m_t) & (P > P.max() * 1e-6)
        m2, sm, sz, _ = comp_gram(A, B, band)
        print(f"{'孤立m=%+d' % m_t:>13} {'m':>4} {m2:>11.3e} {sm:>9.3f} {sz:>+7.3f}")
        out["part1"].append({"state": f"iso_m{m_t}", "band": f"m={m_t}",
                              "mass2": m2, "S": sm, "sz": sz})
    # 海入り m=1..7
    for m_t in range(1, 8):
        a = shift_eta(a1, m_t - 1) + sea_a; b = shift_eta(b1, m_t - 1) + sea_b
        A, B, _, _ = capture(a, b)
        P = np.mean(np.abs(A) ** 2 + np.abs(B) ** 2, axis=0) / 2
        band = (mm[None, :] == m_t) & (P > P.max() * 1e-6)
        if band.sum() == 0:
            print(f"{'海中m=%+d' % m_t:>13}: 帯空"); continue
        m2, sm, sz, wt = comp_gram(A, B, band)
        band0 = (mm[None, :] == 0) & (P > P.max() * 1e-6)
        m2s, sms, szs, _ = comp_gram(A, B, band0)
        print(f"{'海中m=%+d' % m_t:>13} {'m':>4} {m2:>11.3e} {sm:>9.3f} {sz:>+7.3f}"
              f"   (同状態の海帯: 質量²={m2s:.2e})")
        out["part1"].append({"state": f"sea_m{m_t}", "band": f"m={m_t}",
                              "mass2": m2, "S": sm, "sz": sz,
                              "sea_band_mass2": m2s, "band_power": wt})

    print("\n== Part2: スピン量子数読出し（SU(2)回帰角判別） ==")
    thetas = np.linspace(0, 4 * np.pi, 33)
    results2 = {}
    for label, (aa, bb) in (("帯電census(D)", (a0.copy(), b0.copy())),
                             ("海単独", (sea_a.copy(), sea_b.copy()))):
        for _ in range(J_SETTLE):
            aa, bb, _ = ex.collision_step_exact(aa, bb, sp)
        fa0 = np.fft.fft(np.fft.fft(aa.reshape(shape), axis=0, norm="ortho"), axis=1, norm="ortho")
        fb0 = np.fft.fft(np.fft.fft(bb.reshape(shape), axis=0, norm="ortho"), axis=1, norm="ortho")
        for bname, bmask in (("F帯(χ偶≥4)", ferm_k), ("B帯(χ奇)", bos_k)):
            sel = bmask[:, None] & np.ones((1, ne), bool)
            A0 = fa0[sel]; B0 = fb0[sel]
            norm = float(np.sum(np.abs(A0) ** 2 + np.abs(B0) ** 2))
            if norm <= 0: continue
            F = []
            for th in thetas:
                c, s_ = np.cos(th / 2), np.sin(th / 2)
                Ath = c * A0 - s_ * B0
                Bth = s_ * A0 + c * B0
                ov = np.sum(np.conj(A0) * Ath + np.conj(B0) * Bth)
                F.append(float(ov.real) / norm)
            F = np.array(F)
            c_half = 2 * np.mean(F * np.cos(thetas / 2))
            c_one = 2 * np.mean(F * np.cos(thetas))
            w_spinor = abs(c_half) / (abs(c_half) + abs(c_one) + 1e-300)
            print(f"{label:>14} {bname}: F(2π)={F[16]:+.4f} F(4π)={F[32]:+.4f} "
                  f"c_θ/2={c_half:+.4f} c_θ={c_one:+.4f} スピノル重み={w_spinor:.3f}")
            results2[f"{label}|{bname}"] = {"F_2pi": float(F[16]), "F_4pi": float(F[32]),
                "c_half": float(c_half), "c_one": float(c_one),
                "w_spinor": float(w_spinor), "F_curve": F.tolist()}
    out["part2"] = results2
    out["thetas_n"] = 33
    out["runtime_sec"] = time.time() - t0
    (HERE / "pre_v10_instruments_result_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s → pre_v10_instruments_result_v1.json")

if __name__ == "__main__":
    main()
