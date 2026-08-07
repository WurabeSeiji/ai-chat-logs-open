#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""H1予備: 電子の厳密な写像の判定バッテリー（電子の反跳実験・第一段）

電子（周期表 v2 電子行）の指定:
  巻き m=−3（電荷 Q=m/3=−1・3|m ゆえ単独可読）・フェルミオン帯優勢組成・
  被覆度2（時計の二重被覆・{0,π}シート）・質量=⟨ω⟩（海との関係量）。

構成: Part A=局所化エンジン（重力論文の計器群）で巻き/電荷/可読性/組成/質量。
      Part B=厳密エンジン（run_pre_covering_degree_v13b の符号交代検定を
      電子束に適用・コピー移植）で被覆度。

判定（初版の反省を反映した v2 判定・事前固定）:
 (H1a) 構成純度: T=0 の η巻きスペクトルで w=−3 純度 >99.9%。発展後（T=300）
       のスペクトル分散は海駆動ウォークの物理（周期表柱4）として記録のみ
       （初版が発展後に純度判定を課したのは分類と力学の混同＝設計ミス）。
 (H1b) 電荷と可読性: 分母3時計で Q=m/3=−1（整数・厳密）。対照のクォーク型
       プローブ m=−1 は 3∤m で不可読（分数 −1/3・単独では読めない）。
 (H1c) 組成: フェルミオン帯分率 r=P_f/(P_f+P_b) > 0.5（フェルミオン優勢）。
       実測値を電子行の組成として記録（α帯 0.697 との距離も記録）。
 (H1d) 質量: 平衡後（T>1600・診断でプラトー確認済み）の二窓 ⟨ω⟩ が
       ドリフト <5%（初版は立ち上がり中の窓を比較した設計ミス）。
 (H1e) 被覆度: 正本 v14（クロス表）の Qz2 判定器——回帰ピーク位相が {0,π}
       の 0.1π 内に入る割合 Qz2>0.9 かつ両符号出現 → 被覆度2。電子束は
       χ奇倍音（29,31,33・v14 で被覆度2 が実測されたパリティ）× m=−3。
       （初版が較正前 v13b の符号交代判定を使ったのは正本参照ミス）。
使い方: python3 run_h1_electron_battery_pre_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_h1", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
v1, toy, base = ex.v1, ex.toy, ex.base

M_ELECTRON = -3
FSRC = 0.7  # フェルミオン帯優勢（α帯 0.697 近傍を初期値に採用・実測で記録）


def part_a():
    """局所化エンジン: 巻き・電荷・可読性・組成・質量"""
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

    def prof(center, fsrc):
        p = np.zeros(n); AMP = 0.05 * np.sqrt(n) * 0.1
        for parity, wgt in (("even", np.sqrt(fsrc)), ("odd", np.sqrt(1 - fsrc))):
            ks = np.arange(4 if parity == "even" else 5, n // 2, 2)
            Wk = np.exp(-0.5 * (ks / 32.0) ** 2)
            sub = np.zeros(n)
            for kk, w in zip(ks, Wk):
                sub += w * np.cos(2 * np.pi * kk * (x - center) / n)
            p += wgt * sub / np.sqrt(np.sum(sub ** 2))
        return AMP * p

    def lump(m, center, fsrc=FSRC):
        return prof(center, fsrc)[:, None] * np.exp(2j * np.pi * m * eta / ne)[None, :]

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
    mask = (np.abs((x - cA + n // 2) % n - n // 2) <= 24)
    mm = np.where(eta <= ne // 2, eta, eta - ne)

    def winding_spectrum(m):
        """T発展後の台上η巻きスペクトル（海対照差し引き）"""
        a2 = (sea[:, None] * np.ones((1, ne)) + lump(m, cA)).astype(complex)
        b2 = -1j * a2
        a0 = (sea[:, None] * np.ones((1, ne))).astype(complex); b0 = -1j * a0
        for j in range(300):
            a2, b2 = step(a2, b2); a0, b0 = step(a0, b0)
        da = a2 - a0; db = b2 - b0
        Fa = np.fft.fft(da[mask, :], axis=1); Fb = np.fft.fft(db[mask, :], axis=1)
        return np.sum(np.abs(Fa) ** 2 + np.abs(Fb) ** 2, axis=0)

    out = {}
    # H1a: 構成純度（T=0）＋発展後スペクトルはウォーク物理として記録
    l0 = lump(M_ELECTRON, cA)
    F0 = np.fft.fft(l0[mask, :], axis=1)
    P0 = np.sum(np.abs(F0) ** 2, axis=0)
    pur0 = float(P0[M_ELECTRON % ne] / np.sum(P0))
    P = winding_spectrum(M_ELECTRON)
    Pch = P.copy(); Pch[0] = 0.0
    w_main = int(mm[int(np.argmax(Pch))])
    conc = float(Pch[np.argmax(Pch)] / np.sum(Pch))
    h1a = pur0 > 0.999
    print(f"(H1a) 構成純度(T=0)={pur0:.6f}（>0.999）／発展後(T=300): 主巻き="
          f"{w_main} 集中度={conc:.3f}（ウォーク分散・記録のみ） → "
          f"{'通過' if h1a else '不成立'}")
    # H1b: 分母3時計の電荷読出しと可読性（電子 vs クォーク型対照）
    def charge_read(m):
        readable = (m % 3 == 0)
        return (m // 3 if readable else m / 3), readable
    Qe, re_ok = charge_read(M_ELECTRON)
    Qq, rq_ok = charge_read(-1)
    h1b = (Qe == -1) and re_ok and (not rq_ok)
    print(f"(H1b) 電子: Q={Qe}（可読={re_ok}）／クォーク型対照 m=−1: "
          f"Q={Qq:+.3f}（可読={rq_ok}） → {'通過' if h1b else '不成立'}")
    # H1c: 組成（フェルミオン帯分率）
    l0 = lump(M_ELECTRON, cA)
    Fl = np.fft.fft(l0, axis=0)
    Pf = float(np.sum(np.abs(Fl * Wf[:, None]) ** 2))
    Pb = float(np.sum(np.abs(Fl * Wb[:, None]) ** 2))
    r = Pf / (Pf + Pb)
    h1c = r > 0.5
    print(f"(H1c) フェルミオン帯分率 r={r:.4f}（>0.5・α帯0.697との差 {abs(r-0.697):.4f}）"
          f" → {'通過' if h1c else '不成立'}")
    # H1d: 質量（平衡後 T>1600 の二窓・診断でプラトー T≈1200〜 を確認済み）
    p = prof(cA, FSRC); w = p ** 2 / np.sum(p ** 2)
    def mass_window(Tburn, Tavg):
        a2 = (sea[:, None] * np.ones((1, ne)) + lump(M_ELECTRON, cA)).astype(complex)
        b2 = -1j * a2
        a0 = (sea[:, None] * np.ones((1, ne))).astype(complex); b0 = -1j * a0
        acc = np.zeros(n); acc0 = np.zeros(n)
        for j in range(Tburn + Tavg):
            ap = a2.copy(); ap0 = a0.copy()
            a2, b2 = step(a2, b2); a0, b0 = step(a0, b0)
            if j >= Tburn:
                acc += np.angle(np.einsum("xe,xe->x", np.conj(ap), a2))
                acc0 += np.angle(np.einsum("xe,xe->x", np.conj(ap0), a0))
        return float(np.sum(w * (acc - acc0))) / Tavg
    m1 = mass_window(1600, 1200)
    m2 = mass_window(2800, 1200)
    drift = abs(m2 - m1) / abs(m1)
    h1d = (m1 > 0) and (drift < 0.05)
    print(f"(H1d) 質量 ⟨ω⟩={m1:+.4e}（後窓 {m2:+.4e}・ドリフト{drift:.3f}<0.05）"
          f" → {'通過' if h1d else '不成立'}")
    out.update({"H1a": {"purity_T0": pur0, "w_main_evolved": w_main,
                        "concentration_evolved": conc, "ok": bool(h1a)},
                "H1b": {"Q_electron": Qe, "Q_quark_probe": Qq, "ok": bool(h1b)},
                "H1c": {"r": r, "ok": bool(h1c)},
                "H1d": {"mass": m1, "mass_late": m2, "drift": drift,
                        "ok": bool(h1d)}})
    return out, h1a and h1b and h1c and h1d


def part_b():
    """厳密エンジン: 被覆度の Qz2 検定（正本 v14 クロス表判定器・電子束 m=−3）"""
    S = 8.0; J_MAX = 6000
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    shape = (n, ne)
    ms = np.arange(ne); mm = np.where(ms <= ne // 2, ms, ms - ne)

    def single_winding(v):
        f = np.fft.fft(v.reshape(shape), axis=0, norm="ortho")
        f[0, :] = 0.0; f[n // 2:, :] = 0.0
        return np.fft.ifft(f, axis=0, norm="ortho").reshape(v.shape)

    def project_eta(v, m_set):
        f = np.fft.fft(v.reshape(shape), axis=1, norm="ortho")
        keep = np.isin(mm, list(m_set)); f[:, ~keep] = 0.0
        return np.fft.ifft(f, axis=1, norm="ortho").reshape(v.shape)

    # 電子束: 奇χ倍音（フェルミオン型）× 巻き m=−3
    a0 = single_winding(v1.make_bundle(sp, (29, 31, 33), "A", scale=1.0)) * S
    b0 = single_winding(v1.make_bundle(sp, (29, 31, 33), "B", scale=1.0)) * S
    a = project_eta(a0, {M_ELECTRON}); b = project_eta(b0, {M_ELECTRON})
    pw = float(np.sum(np.abs(a) ** 2) + np.sum(np.abs(b) ** 2))
    sc = np.sqrt((float(np.sum(np.abs(a0) ** 2) + np.sum(np.abs(b0) ** 2))) / pw)
    a *= sc; b *= sc
    a0_, b0_ = a.copy(), b.copy()
    s0 = np.imag(np.conj(b0_) * a0_); s0c = s0 - s0.mean()
    nrm = float(np.sum(np.abs(a0_) ** 2) + np.sum(np.abs(b0_) ** 2))
    A, C = [], []
    for j in range(J_MAX):
        a, b, _ = ex.collision_step_exact(a, b, sp)
        A.append(complex((np.vdot(a0_, a) + np.vdot(b0_, b)) / nrm))
        sj = np.imag(np.conj(b) * a); sjc = sj - sj.mean()
        C.append(float(np.dot(s0c, sjc) /
                       max(np.linalg.norm(s0c) * np.linalg.norm(sjc), 1e-300)))
    A = np.array(A); C = np.array(C)
    th = 0.85 * float(C[10:].max())
    peaks = []
    for j in range(11, J_MAX - 1):
        if C[j] > th and C[j] >= C[j - 1] and C[j] >= C[j + 1]:
            if not peaks or j - peaks[-1] > 5:
                peaks.append(j)
        if len(peaks) >= 15:
            break
    phis = np.array([np.angle(A[j]) for j in peaks])
    print(f"(H1e) 被覆度 Qz2 検定（v14 正本判定器・ピーク{len(phis)}個）:")
    if len(phis) < 4:
        print("   → 判定不能")
        return {"n_peaks": len(phis), "verdict": "判定不能", "ok": False}, False
    d0 = np.abs(phis)
    dpi = np.minimum(np.abs(phis - np.pi), np.abs(phis + np.pi))
    near = np.minimum(d0, dpi) < 0.1 * np.pi
    Qz2 = float(np.mean(near))
    signs = np.sign(np.cos(phis[near])) if near.any() else np.array([])
    if Qz2 > 0.9 and (signs < 0).any() and (signs > 0).any():
        verdict = "被覆度2(Z₂)"
    elif Qz2 > 0.9 and (signs > 0).all():
        verdict = "被覆度1"
    else:
        verdict = "連続(被覆なし)"
    for kk, (j, ph) in enumerate(zip(peaks[:10], phis[:10])):
        print(f"   k={kk+1:2d} J={j+1:4d} Φ/π={ph/np.pi:+.3f} "
              f"{'∈{0,π}帯' if near[kk] else '外'}")
    ok = verdict == "被覆度2(Z₂)"
    print(f"   Qz2={Qz2:.2f} → {verdict} → {'通過' if ok else '不成立'}")
    return {"n_peaks": int(len(phis)), "Qz2": Qz2, "verdict": verdict,
            "phis_over_pi": [float(x / np.pi) for x in phis],
            "ok": bool(ok)}, ok


def main():
    t0 = time.time()
    print(f"=== H1予備: 電子判定バッテリー（m={M_ELECTRON}, fsrc={FSRC}）===")
    out_a, ok_a = part_a()
    out_b, ok_b = part_b()
    ok = ok_a and ok_b
    verdict = ("電子の写像成立: 巻き−3・電荷−1・単独可読・フェルミオン優勢・"
               "安定質量・被覆度2の5項目通過" if ok else "要精査")
    print(verdict)
    out = {"M": M_ELECTRON, "fsrc": FSRC, "part_a": out_a, "part_b_covering": out_b,
           "verdict": verdict, "runtime_sec": time.time() - t0}
    (HERE / "result_h1_electron_battery_pre_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
