#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""G9f: 超選択セクターの長時間検定（§18-10 の実施）

M3（G9d）: 倍加軌道が海(0)に到達 ⟺ q|m（q=ne の奇数部分）。
ne=16: 全電荷可溶（1→2→4→8→0）。
ne=12: m=1 は永続（1→2→4→8→4 周期軌道・海に届かない）、m=3 は溶解
（3→6→12≡0）。

設計の要（厳密な内部対照）: ne=16 では 3∈Z₁₆^*（単元）ゆえ m=1 と m=3 は
単元自己同型により厳密縮退する（G9e で実証済）。ne=12 では 3 は単元でない
（gcd(3,12)=3）。ゆえに M3 が動力学的に実在するなら:

判定（初版L1は浮動小数点カオスを無視した設計ミス→較正形に修正・反証記録）:
 (L1') セクター分岐がカオス床を超える: ne=16 の単元対 (m=1, m=3) は厳密算術
      では自己同型により縮退するが、η和の加算順序差がカオス増幅され数値床
      deg16 を作る（T=8000 実測 ~0.9%）。これを対照床として、
      ne=12 の分岐 split12 > 10×deg16 を要求する。
 (L2) 分岐の向き: ne=12 で溶解種 m=3 の荷電パワーは永続種 m=1 より
      速く減る（後半窓の傾き slope(m=3) < slope(m=1)）。
 (L3) 周期軌道蓄積: 後半窓平均で C₄₈(m=1) > C₄₈(m=3)。
 (L4) 厳密部分群閉じ込め（M1の動力学版）: m=3@12 の {4,8} 内容は全時間で
      機械零（<1e-20）——内容は部分群 ⟨3⟩={0,3,6,9} を厳密に出ない。
使い方: python3 run_g9f_superselection_longtime_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_g9f", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
base = ex.base

T_TOTAL = 8000
SNAP = 400


def trace_for(ne, m, n):
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

    cA = 100
    prof = np.zeros(n); AMP = 0.05 * np.sqrt(n) * 0.1
    for parity, wgt in (("even", np.sqrt(0.6)), ("odd", np.sqrt(0.4))):
        ks = np.arange(4 if parity == "even" else 5, n // 2, 2)
        Wk = np.exp(-0.5 * (ks / 32.0) ** 2)
        sub = np.zeros(n)
        for kk, w in zip(ks, Wk):
            sub += w * np.cos(2 * np.pi * kk * (x - cA) / n)
        prof += wgt * sub / np.sqrt(np.sum(sub ** 2))
    prof *= AMP
    lump = prof[:, None] * np.exp(2j * np.pi * m * eta / ne)[None, :]
    mask = (np.abs((x - cA + n // 2) % n - n // 2) <= 24)
    a2 = (sea[:, None] * np.ones((1, ne)) + lump).astype(complex); b2 = -1j * a2
    a0 = (sea[:, None] * np.ones((1, ne))).astype(complex); b0 = -1j * a0
    Ts, As, C48s = [], [], []
    for j in range(T_TOTAL + 1):
        if j % SNAP == 0:
            da = a2 - a0; db = b2 - b0
            Fa = np.fft.fft(da[mask, :], axis=1); Fb = np.fft.fft(db[mask, :], axis=1)
            P = np.sum(np.abs(Fa) ** 2 + np.abs(Fb) ** 2, axis=0)
            Ts.append(j)
            As.append(float(np.sum(P) - P[0]))
            C48s.append(float(P[4 % ne] + P[8 % ne]))
        if j < T_TOTAL:
            a2, b2 = step(a2, b2); a0, b0 = step(a0, b0)
    return np.array(Ts), np.array(As), np.array(C48s)


def late_slope(Ts, ys, frac=0.5):
    i0 = int(len(Ts) * frac)
    A = np.vstack([np.ones_like(Ts[i0:], dtype=float), Ts[i0:]]).T
    coef, *_ = np.linalg.lstsq(A, ys[i0:], rcond=None)
    return float(coef[1])


def main():
    t0 = time.time()
    n = base.build_source_params(base.Params(high_n=63,
                                             recursive_collision_count=200)).chi_grid_n
    tr = {}
    for ne, m in [(16, 1), (16, 3), (12, 1), (12, 3)]:
        Ts, As, C48s = trace_for(ne, m, n)
        tr[(ne, m)] = (Ts, As, C48s)
        print(f"ne={ne} m={m}: A(0)={As[0]:.4e} A(末)={As[-1]:.4e} "
              f"A末/A初={As[-1]/As[0]:.4f}  C48末={C48s[-1]:.3e}")
    i0 = len(tr[(16, 1)][0]) // 2
    # L1: ne=16 厳密縮退 vs ne=12 分岐（後半窓のA平均で比較）
    m16 = [float(np.mean(tr[(16, mm)][1][i0:])) for mm in (1, 3)]
    m12 = [float(np.mean(tr[(12, mm)][1][i0:])) for mm in (1, 3)]
    deg16 = abs(m16[0] - m16[1]) / m16[0]
    split12 = abs(m12[0] - m12[1]) / m12[0]
    l1 = split12 > 10.0 * deg16
    # L2: 分岐の向き（後半窓の傾き）
    s12_1 = late_slope(*tr[(12, 1)][:2])
    s12_3 = late_slope(*tr[(12, 3)][:2])
    l2 = s12_3 < s12_1
    # L3: 周期軌道への蓄積
    c1 = float(np.mean(tr[(12, 1)][2][i0:]))
    c3 = float(np.mean(tr[(12, 3)][2][i0:]))
    l3 = c1 > c3
    c3max = float(np.max(np.abs(tr[(12, 3)][2])))
    l4 = c3max < 1e-20
    print(f"\n(L1') セクター分岐vs.カオス床: ne=16 床={deg16:.2e}・"
          f"ne=12 分岐={split12:.4f}（>10×床={10*deg16:.4f}） → {'通過' if l1 else '不成立'}")
    print(f"(L2) 分岐の向き: slope(m=3)={s12_3:+.3e} < slope(m=1)={s12_1:+.3e}"
          f" → {'通過' if l2 else '不成立'}")
    print(f"(L3) 周期軌道蓄積: C48(m=1)={c1:.3e} > C48(m=3)={c3:.3e}"
          f" → {'通過' if l3 else '不成立'}")
    print(f"(L4) 厳密部分群閉じ込め: max|C48(m=3@12)|={c3max:.2e}（<1e-20）"
          f" → {'通過' if l4 else '不成立'}")
    ok = l1 and l2 and l3 and l4
    verdict = ("超選択セクターの動力学的実在を確認: 部分群閉じ込めは厳密（機械零）・"
               "溶解/永続の分岐はカオス床の10倍超（M1/M3 の動力学版・ne=構造定数）"
               if ok else "要精査")
    print(verdict)
    out = {"traces": {f"ne{ne}_m{m}": {"T": tr[(ne, m)][0].tolist(),
                                       "A": tr[(ne, m)][1].tolist(),
                                       "C48": tr[(ne, m)][2].tolist()}
                      for (ne, m) in tr},
           "L1_deg16": deg16, "L1_split12": split12,
           "L2_slope12_m1": s12_1, "L2_slope12_m3": s12_3,
           "L3_C48_m1": c1, "L3_C48_m3": c3,
           "L4_c3max": c3max,
           "L1": bool(l1), "L2": bool(l2), "L3": bool(l3), "L4": bool(l4),
           "verdict": verdict, "runtime_sec": time.time() - t0}
    (HERE / "result_g9f_superselection_longtime_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
