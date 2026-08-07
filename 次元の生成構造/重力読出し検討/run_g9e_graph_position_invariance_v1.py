#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""G9e: 倍加グラフ位置の不変性——ne=12 超選択予言の初期実測

G9d の定理 M3: 倍加軌道が海(0)に到達 ⟺ q|m（q=ne の奇数部分）。
ne=16 は全電荷が海に溶け、ne=12 では 3∤m の電荷が永続セクターになる。

予言（事前固定・グラフ位置仮説）: solo源の η スペクトル署名は m や ne の
値でなく、倍加グラフ上の位置（海までの距離・周期軌道・共役との合流）で
決まる。
 (E1) 「海から2歩」の同値: spec(m=4, ne=16) = spec(m=3, ne=12)
      ——一次/共役/倍加の相対重み3成分が相対差<2%で一致。
 (E2) 合流位置の固有署名: ne=12 の m=4 は −4≡8≡2·4（共役と倍加が同一巻きに
      合流）ゆえ、一次保持が一般位置より>20%増強され、合流巻きの重みは
      一般位置の共役+倍加の和より小さい（干渉再編）。
 (E3) 一般位置の普遍性: spec(m=1, ne=12) = spec(m=1, ne=16)（3成分<2%）
      ——海まで≥3歩の位置は局所グラフが同型で署名が共通。
使い方: python3 run_g9e_graph_position_invariance_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_g9e", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
base = ex.base


def spectrum_for(ne, m, n, T=700, snaps=(500, 600, 700)):
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
    P = np.zeros(ne); cnt = 0
    for j in range(T + 1):
        if j in snaps:
            da = a2 - a0; db = b2 - b0
            Fa = np.fft.fft(da[mask, :], axis=1); Fb = np.fft.fft(db[mask, :], axis=1)
            P += np.sum(np.abs(Fa) ** 2 + np.abs(Fb) ** 2, axis=0); cnt += 1
        if j < T:
            a2, b2 = step(a2, b2); a0, b0 = step(a0, b0)
    return P / cnt


def sig3(P, ne, m):
    """署名3成分: (一次, 共役, 倍加) の相対重み（w=0 海変形を除く全パワーで規格化）"""
    tot = np.sum(P) - P[0]
    return (float(P[m % ne] / tot), float(P[(-m) % ne] / tot),
            float(P[(2 * m) % ne] / tot))


def main():
    t0 = time.time()
    n = base.build_source_params(base.Params(high_n=63,
                                             recursive_collision_count=200)).chi_grid_n
    S = {}
    for ne, m in [(16, 1), (16, 4), (12, 1), (12, 3), (12, 4)]:
        P = spectrum_for(ne, m, n)
        S[(ne, m)] = sig3(P, ne, m)
        print(f"ne={ne} m={m}: (一次,共役,倍加)=({S[(ne,m)][0]:.4f},"
              f"{S[(ne,m)][1]:.4f},{S[(ne,m)][2]:.4f})")

    def rel3(a, b):
        return max(abs(x - y) / max(x, 1e-12) for x, y in zip(a, b))

    e1 = rel3(S[(16, 4)], S[(12, 3)])
    e3 = rel3(S[(16, 1)], S[(12, 1)])
    prim_gen = S[(12, 1)][0]
    prim_conf = S[(12, 4)][0]
    conf_merge = S[(12, 4)][2]           # 合流巻き（8=−4=2·4）の重み
    gen_sum = S[(12, 1)][1] + S[(12, 1)][2]
    e2 = (prim_conf > 1.2 * prim_gen) and (conf_merge < gen_sum)
    print(f"\n(E1) 海から2歩の同値 spec(4@16)=spec(3@12): 相対差={e1:.4f}（判定<0.02）"
          f" → {'通過' if e1 < 0.02 else '不成立'}")
    print(f"(E2) 合流位置の固有署名: 一次 {prim_conf:.4f} vs 一般 {prim_gen:.4f}"
          f"（>1.2倍）・合流巻き {conf_merge:.4f} < 共役+倍加 {gen_sum:.4f}"
          f" → {'通過' if e2 else '不成立'}")
    print(f"(E3) 一般位置の普遍性 spec(1@12)=spec(1@16): 相対差={e3:.4f}（判定<0.02）"
          f" → {'通過' if e3 < 0.02 else '不成立'}")
    ok = (e1 < 0.02) and e2 and (e3 < 0.02)
    verdict = ("グラフ位置仮説成立: η スペクトル署名は倍加グラフ上の位置の不変量"
               "（M定理群の物理的実在の初期実証・超選択の長時間検定は§18登録）"
               if ok else "要精査")
    print(verdict)
    out = {"signatures": {f"ne{ne}_m{m}": list(v) for (ne, m), v in S.items()},
           "E1_rel": e1, "E2": bool(e2), "E3_rel": e3,
           "verdict": verdict, "runtime_sec": time.time() - t0}
    (HERE / "result_g9e_graph_position_invariance_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
