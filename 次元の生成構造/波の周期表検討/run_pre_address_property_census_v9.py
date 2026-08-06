#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""予備実験 v9: 番地×性質census——各巻き数番地のスピン・質量・分類の実測

目的: 周期表第0版の空欄（スピン・質量・粒子割当）を実測で埋める。
計器（空間時間論文で確立済みの読出しをそのまま適用）:
  質量² = detΓ/T²（二チャネルGramの非コヒーレンス）
  スピン: |S| = √(X²+Y²+Z²)/T（Blochベクトル強度）、s_z = Z/T
  分類: χパリティ（|k|偶≥4=フェルミオン的 / 奇=ボゾン的）の帯パワー比
状態: 純粋+1帯電構成（v4bで厳密安定を確認済み）を e^{i(m−1)η} でシフトし
純粋m状態（m=0..8）を作る。J=160整定後、窓40衝突でGram測定。
帯: パワー上位の (k,m) モードを対象、パワー加重平均で番地の代表値とする。

割当規則（事前記録・候補判定）:
  無質量(質量²<1e-3)・中性(fold3=0)・コヒーレント → 光子的
  有質量・読電荷±1・フェルミオン帯優勢 → 電子/陽電子的
  読電荷0(mod3中性)・|m|=3(3構成子複合) → 中性子的候補
  その他 → 性質を記録し未割当（無理な割当をしない）

使い方: python3 run_pre_address_property_census_v9.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_v9", UIM / "run_ignition_fate_exact_v3.py")
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

    def fold3(m): return ((m + 1) % 3) - 1

    out = {"J_SETTLE": J_SETTLE, "J_WIN": J_WIN, "ne": int(ne), "addresses": []}
    print(f"{'m':>3} {'読電荷':>5} {'質量²':>10} {'|S|':>7} {'s_z':>7} "
          f"{'F帯比':>6} {'割当候補':>12}")
    for m_t in range(0, 9):
        a = shift_eta(a1, m_t - 1); b = shift_eta(b1, m_t - 1)
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
        Gaa = np.mean(np.abs(A) ** 2, axis=0); Gbb = np.mean(np.abs(B) ** 2, axis=0)
        Gab = np.mean(A * np.conj(B), axis=0)
        T = (Gaa + Gbb) / 2
        X = np.real(Gab); Y = np.imag(np.conj(Gab)); Z = (Gaa - Gbb) / 2
        det = np.real(Gaa * Gbb - np.abs(Gab) ** 2)
        Tf = np.maximum(T, 1e-300)
        mass2 = det / Tf ** 2
        smag = np.sqrt(X ** 2 + Y ** 2 + Z ** 2) / Tf
        sz = Z / Tf
        # 対象帯: 巻き数 m_t のモードでパワー上位
        band = (mm[None, :] == m_t) & (P > P.max() * 1e-6)
        wsum = float(np.sum(P[band]))
        if wsum <= 0:
            print(f"{m_t:>+3}  帯が空——スキップ"); continue
        w = P[band] / wsum
        m2 = float(np.sum(w * mass2[band]))
        sm = float(np.sum(w * smag[band]))
        szv = float(np.sum(w * sz[band]))
        fb_ = float(np.sum(P[band & ferm_k[:, None]]) / wsum)
        q3 = fold3(m_t)
        # 割当候補（事前規則）
        if m2 < 1e-3 and q3 == 0:
            assign = "光子的"
        elif abs(q3) == 1 and fb_ > 0.5 and m2 > 1e-3:
            assign = "電子族" if q3 == +1 else "陽電子族"
        elif q3 == 0 and abs(m_t) == 3:
            assign = "中性子的候補"
        else:
            assign = "未割当"
        print(f"{m_t:>+3} {q3:>+4} {m2:>10.2e} {sm:>7.3f} {szv:>+7.3f} "
              f"{fb_:>6.3f} {assign:>12}")
        out["addresses"].append({"m": m_t, "q_read": q3, "mass2": m2,
                                  "spin_mag": sm, "sz": szv, "ferm_frac": fb_,
                                  "assign": assign, "band_power": wsum})

    # ---- 対照測定（v9b・同一計器）----
    def measure_state(a, b, settle, label):
        for _ in range(settle):
            a, b, _ = ex.collision_step_exact(a, b, sp)
        A = np.zeros((J_WIN, n, ne), complex); B = np.zeros((J_WIN, n, ne), complex)
        for t in range(J_WIN):
            a, b, _ = ex.collision_step_exact(a, b, sp)
            fa = np.fft.fft(a.reshape(shape), axis=0, norm="ortho")
            fb = np.fft.fft(b.reshape(shape), axis=0, norm="ortho")
            A[t] = np.fft.fft(fa, axis=1, norm="ortho")
            B[t] = np.fft.fft(fb, axis=1, norm="ortho")
        P = np.mean(np.abs(A) ** 2 + np.abs(B) ** 2, axis=0) / 2
        Gaa = np.mean(np.abs(A) ** 2, axis=0); Gbb = np.mean(np.abs(B) ** 2, axis=0)
        Gab = np.mean(A * np.conj(B), axis=0)
        T = (Gaa + Gbb) / 2
        Xv = np.real(Gab); Yv = np.imag(np.conj(Gab)); Zv = (Gaa - Gbb) / 2
        det = np.real(Gaa * Gbb - np.abs(Gab) ** 2)
        Tf = np.maximum(T, 1e-300)
        band = P > P.max() * 1e-6
        w = P[band] / P[band].sum()
        m2 = float(np.sum(w * (det / Tf ** 2)[band]))
        sm = float(np.sum(w * (np.sqrt(Xv**2+Yv**2+Zv**2) / Tf)[band]))
        fb_ = float(np.sum(P[band & ferm_k[:, None]]) / P[band].sum())
        print(f"[対照 {label}] 質量²={m2:.2e} |S|={sm:.3f} F帯比={fb_:.3f}")
        out["controls"] = out.get("controls", [])
        out["controls"].append({"label": label, "mass2": m2, "spin_mag": sm,
                                 "ferm_frac": fb_, "settle": settle})
    measure_state(a1.copy(), b1.copy(), 0, "コヒーレント+1（整定なし=光子的対照）")
    measure_state(a1.copy(), b1.copy(), 2000, "長期熱化+1（海化対照）")
    measure_state(a1.copy(), a1.copy(), 0, "A-Bロック+1（位相ロック複製=光子的対照）")
    measure_state(a1.copy(), a1.copy(), 200, "A-Bロック+1（整定200後）")
    out["runtime_sec"] = time.time() - t0
    (HERE / "pre_address_property_census_result_v9.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s → pre_address_property_census_result_v9.json")

if __name__ == "__main__":
    main()
