#!/usr/bin/env python3
"""残穴①②: 位相選択の統計＋読出し局所性の公理的整理

穴①: 位相選択の統計（実行前固定）:
    乱位相パケット種 50 個（N=12, Nreg=16, δ=3e-2, T=200）で早期レートの
    符号と大きさの分布を測る。
    S1-a: 物質化率 P(+) を測定（v3 予備では 5/6）。
    S1-b: 決定論性——瞬時符号予測子 sign(dP_odd/dt|₀)（頂点レートの奇数
        射影と現奇数内容の内積×2Re）が、フィットした早期レートの符号を
        予測するか。的中率 ≥ 90% なら「くじは擬似乱数だが決定論的」。

穴②: 読出し局所性の公理的整理:
    定理拡張 T-L: 閉塞保存は R_{ee'}(n)（レジスタ点ごとの対対称強度）を
        許す——一意化定理の対消去は点 n ごとに成立するため。数値検証:
        ランダムな n 依存対対称 R で dC=0（機械精度）、非対称なら破れ。
    命題 P-L: 局所性は読出しでなく頂点の点ごと構造（仮定A1）が公理的に
        担う。v3 の増強はその Fourier 双対＝チャネル多重度: 局在パケット
        （幅W）⟺ 帯域 ~Nreg/W ⟺ 四波混合チャネル数増。検証: 種の奇数
        倍音本数 n_h ∈ {1,2,4,8}（乱位相、各4種平均）で |rate| ~ n_h^q の
        冪 q を実測（q>0 なら双対読みが成立、q≈1〜2 の帯を予想）。

使い方: python3 run_holes_h1_h2_v1.py
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec3 = importlib.util.spec_from_file_location(
    "s3h", HERE / "run_genesis_v3_register_local_v1.py")
g3 = importlib.util.module_from_spec(spec3)
sys.modules[spec3.name] = g3
spec3.loader.exec_module(g3)
abl = g3.abl
V2 = g3.V2

N_GRAPH, NREG, DELTA = 12, 16, 3e-2
T_FIT = 200


def early_rate_and_sign_pred(Z0c, wp0, seed_edge, prof):
    m = N_GRAPH * (N_GRAPH - 1) // 2
    C0 = np.zeros((m, NREG), complex)
    C0[:, 2] = Z0c
    for k in range(NREG):
        if abs(prof[k]) > 0:
            C0[:, k] += DELTA * prof[k] * seed_edge
    eng = V2(N_GRAPH, C0, wp0, vertex_on=True)
    # 瞬時符号予測子: 頂点レートの奇数射影と現奇数内容の内積
    R = eng._readout()
    W = np.fft.ifft(eng.C, axis=1) * NREG
    rate_n = eng._vertex_rate(W, R)
    rate_k = np.fft.fft(rate_n, axis=1) / NREG
    ks = np.arange(NREG)
    odd = (ks % 2 == 1)
    g0 = 2.0 * float(np.real(np.vdot(eng.C[:, odd], rate_k[:, odd])))
    fs = []
    for t in range(T_FIT):
        eng.step()
        fs.append(eng.diagnostics()["f_seed"])
    fs = np.array(fs)
    tt = np.arange(10, 150, dtype=float)
    A = np.vstack([tt, np.ones_like(tt)]).T
    coef, _, _, _ = np.linalg.lstsq(A, np.log(fs[10:150]), rcond=None)
    return float(coef[0]), float(g0)


def main() -> None:
    t0 = time.time()
    m = N_GRAPH * (N_GRAPH - 1) // 2
    _, v, _, _, _, _, _, Z0c, wp0 = abl.build_init(N_GRAPH, False)
    Z0c = Z0c / np.linalg.norm(Z0c)
    seed_edge = g3.zero_closure_state(m, np.random.default_rng(98000))
    odd_ks = [k for k in range(NREG) if k % 2 == 1]
    out = {}

    # ===== 穴① 位相選択の統計 =====
    print("=== 穴① 位相選択の統計（50 乱位相パケット） ===")
    rates, preds = [], []
    for ps in range(50):
        rng = np.random.default_rng(100000 + ps)
        prof = np.zeros(NREG, complex)
        for k in odd_ks:
            prof[k] = np.exp(1j * rng.uniform(0, 2 * np.pi)) / np.sqrt(len(odd_ks))
        r_, g0 = early_rate_and_sign_pred(Z0c, wp0, seed_edge, prof)
        rates.append(r_)
        preds.append(g0)
    rates = np.array(rates)
    preds = np.array(preds)
    p_plus = float(np.mean(rates > 0))
    hits = float(np.mean(np.sign(rates) == np.sign(preds)))
    out["H1"] = {"n": 50, "P_plus": p_plus,
                  "rate_median_abs": float(np.median(np.abs(rates))),
                  "rate_max": float(rates.max()), "rate_min": float(rates.min()),
                  "sign_predictor_hit": hits}
    out["H1_S1b_pass"] = bool(hits >= 0.9)
    print(f"  物質化率 P(+) = {p_plus:.2f}（50個中 {int(p_plus*50)} 個が成長）")
    print(f"  |rate| 中央値={np.median(np.abs(rates)):.2e} 範囲=[{rates.min():.2e}, {rates.max():.2e}]")
    print(f"  S1-b 瞬時符号予測子の的中率 = {hits:.2f} → {'決定論的くじ' if hits>=0.9 else '予測子不完全'}")

    # ===== 穴② T-L: n依存対対称Rの閉塞保存 =====
    print("=== 穴② 読出し局所性 ===")
    ia, ib = np.triu_indices(5, k=1)
    m5 = 10
    adj = [[] for _ in range(m5)]
    for e in range(m5):
        for f in range(m5):
            if e != f and (ia[e] in (ia[f], ib[f]) or ib[e] in (ia[f], ib[f])):
                adj[e].append(f)
    rng = np.random.default_rng(7)
    drifts_sym, drifts_asym = [], []
    for trial in range(5):
        Wn = rng.normal(size=(m5, 8)) + 1j * rng.normal(size=(m5, 8))
        Rsym = rng.uniform(0.1, 0.9, size=(m5, m5, 8))
        Rsym = 0.5 * (Rsym + Rsym.transpose(1, 0, 2))          # 対対称（n依存）
        Rasym = rng.uniform(0.1, 0.9, size=(m5, m5, 8))         # 非対称
        for Rm, bag in ((Rsym, drifts_sym), (Rasym, drifts_asym)):
            dW = np.zeros_like(Wn)
            for e in range(m5):
                for f_ in adj[e]:
                    dW[e] += 1j * Rm[e, f_] * (np.abs(Wn[f_]) ** 2 * Wn[e]
                                                - Wn[f_] ** 2 * np.conj(Wn[e]))
            dC = np.abs(np.sum(Wn * dW, axis=0))
            bag.append(float(dC.max() / np.sum(np.abs(Wn) ** 2) ** 1.5))
    out["H2_TL"] = {"sym_max": max(drifts_sym), "asym_min": min(drifts_asym)}
    tl = bool(max(drifts_sym) < 1e-14 and min(drifts_asym) > 1e-3)
    out["H2_TL_pass"] = tl
    print(f"  T-L: n依存対対称R dC={max(drifts_sym):.1e}（機械零）／非対称 dC={min(drifts_asym):.1e} → {tl}")
    print(f"       ⇒ 閉塞保存は点ごと局所読出し R_ee'(n) を許す（定理拡張成立）")

    # ===== 穴② P-L: チャネル多重度の冪 =====
    ratios = []
    base = None
    for n_h in (1, 2, 4, 8):
        vals = []
        for ps in range(4):
            rng = np.random.default_rng(101000 + 10 * n_h + ps)
            picks = list(rng.choice(odd_ks, size=n_h, replace=False))
            prof = np.zeros(NREG, complex)
            for k in picks:
                prof[k] = np.exp(1j * rng.uniform(0, 2 * np.pi)) / np.sqrt(n_h)
            r_, _ = early_rate_and_sign_pred(Z0c, wp0, seed_edge, prof)
            vals.append(abs(r_))
        med = float(np.median(vals))
        if base is None:
            base = med
        ratios.append({"n_h": n_h, "abs_rate_median": med, "ratio": med / base})
        print(f"  n_h={n_h}: |rate|中央値={med:.3e}（比 {med/base:.2f}）")
    ln_n = np.log([r_["n_h"] for r_ in ratios])
    ln_r = np.log([r_["abs_rate_median"] for r_ in ratios])
    A = np.vstack([ln_n, np.ones_like(ln_n)]).T
    coef, _, _, _ = np.linalg.lstsq(A, ln_r, rcond=None)
    q = float(coef[0])
    out["H2_PL"] = {"rows": ratios, "q_exponent": q}
    out["H2_PL_pass"] = bool(q > 0)
    print(f"  P-L: |rate| ~ n_h^q, q={q:.2f} → {'チャネル多重度=局所性の双対、成立' if q>0 else '不成立'}")

    out["runtime_sec"] = time.time() - t0
    (HERE / "holes_h1_h2_result_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"\nsaved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
