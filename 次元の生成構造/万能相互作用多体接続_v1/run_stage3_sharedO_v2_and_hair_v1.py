#!/usr/bin/env python3
"""段階3: 共有O再設計 v2エンジン＋受入＋v1物理頑健性＋毛(η)拡張census P2完全版

v2設計（媒介頂点md §7 の未決事項の決定）:
    線形部 = 全スライスに同一の Cayley O(θ̄) を適用（共有O）。
    θ̄_e = arg(Σ_k c_e[k])。単一スライス極限で Σ が当該列そのもの（零加算は
    bitwise 恒等）となり対応原理が厳密に成立する。この規約はレジスタ点 n=0
    の読みに等しく（Σ_k c_k = w[0]）、並進アンカーを持つ——無名性上の
    代替案（パワー重み位相平均）は bitwise 対応を壊すため v2 では採らず、
    トレードオフとして記録する。
    共有Oは辺空間の直交行列を全レジスタ点に一様適用するため、
    **点ごと場閉塞 Σ_e w_e[n]² が線形部でも厳密保存**（v1の設計宿題の解消）。

受入基準 v2（実行前固定）:
    B1 対応原理: 単一スライスで abl 軌道と bitwise 一致（T=500）。
    B2 保存則（結合力学 T=200・δ=0.1）: 点ごと場閉塞の max 変化 ≤1e-10、
       スライス毎閉塞 ≤1e-10、全ノルム相対 ≤1e-10。
    B3 物理頑健性（v1 との比較）: E1a の p ∈ [1.9,2.1]・C が v1 値 6.38 の
       [×0.5,×2]。D3 四ポンプの枝順序（+枝2種 > −枝2種）が保存されるか
       ——保存されればカイラリティ効果はアーキテクチャ非依存、
       反転/消失すれば v1 の産物と判定（どちらでも記録）。
    B4 ロック持続: D1-P2 完全予測子コヒーレンスの10step減衰が v1
       （0.9999→0.9936）より浅い（共有Oは差動回転による脱整合を持たない）。

毛(η)拡張 census P2 完全版（実行前固定）:
    レジスタを2次元 (n×η)（Nn=5, Nη=8）に拡張。頂点は点ごと（2D格子の
    各点の積）なので毛の和則 m_partner = 2m_B − m_s が自動成立するはず。
    ポンプ (k=2, m=+2)・種 (k=1, m=+1) → 相棒予言 (k=3, m=+3=2·2−1 mod 8)。
    H1 毛の排他性: T=30 で P_{k=3,m=3} / P_{k=3,m≠3} ≥ 100。
    H2 毛込み対相関: 2D完全予測子（パラメトリック＋ビートXPM、毛込み）
       との一step増分コヒーレンス ≥ 0.99。
    ——二体census P2（毛の帳簿を持つ相棒）の多体完全版。

使い方: python3 run_stage3_sharedO_v2_and_hair_v1.py
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec_e = importlib.util.spec_from_file_location(
    "s2base", HERE / "run_stage2_vertex_engine_v1.py")
s2 = importlib.util.module_from_spec(spec_e)
sys.modules[spec_e.name] = s2
spec_e.loader.exec_module(s2)
abl = s2.abl
gen3 = s2.gen3


class VertexEngineV2(s2.VertexEngine):
    """共有O線形部（単一 wp・単一 LowRankSystem）。頂点部は v1 を継承。"""

    def __init__(self, n, C0, wp, vertex_on=True, strength_scale=1.0):
        super().__init__(n, C0, {}, vertex_on=vertex_on, strength_scale=strength_scale)
        self.wp = wp.copy()
        self.sys_shared = abl.LowRankSystem(n)

    def _linear(self):
        zsum = np.sum(self.C, axis=1)
        self.sys_shared.set_theta(np.angle(zsum))
        se, self.wp = self.sys_shared.sigma_max_power(self.wp)
        for k in range(self.nreg):
            col = self.C[:, k]
            if np.linalg.norm(col) > 0.0:
                self.C[:, k] = self.sys_shared.cayley_step(col, se)


def fval_factory(v0):
    p = v0.real / np.linalg.norm(v0.real)
    q = v0.imag - (v0.imag @ p) * p
    q = q / np.linalg.norm(q)

    def f(Z):
        Zp = Z - p * (p @ Z) - q * (q @ Z)
        return float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Z) @ Z))

    return f


def main() -> None:
    t0 = time.time()
    n, m, nreg, delta = 5, 10, 5, 1e-2
    out = {}
    _, v, _, _, _, _, _, Z0c, wp0 = abl.build_init(n, False)
    r2 = gen3.make_parent(n, seed=2)
    Csec = np.fft.fft(r2.relation_waves, axis=1) / n
    seed_state = Csec[:, 1] / np.linalg.norm(Csec[:, 1])
    vplus = r2.parent_vector / np.linalg.norm(r2.parent_vector)
    ia, ib = np.triu_indices(n, k=1)

    def C0_of(pump, d=delta, reg=nreg):
        C0 = np.zeros((m, reg), complex)
        C0[:, 2] = pump
        C0[:, 1] = d * seed_state
        return C0

    # ---- B1: 対応原理 ----
    T = 500
    sys_ref = abl.LowRankSystem(n)
    sys_ref.set_theta(np.angle(Z0c))
    Z = Z0c.copy()
    w = wp0.copy()
    C0 = np.zeros((m, nreg), complex)
    C0[:, 1] = Z0c
    eng = VertexEngineV2(n, C0, wp0, vertex_on=False)
    ok_b1 = True
    for t in range(T):
        Z, w = abl.evolve(sys_ref, Z, w)
        eng.step()
        if not np.array_equal(eng.C[:, 1], Z):
            ok_b1 = False
            break
    out["B1_pass"] = bool(ok_b1)
    print(f"  B1 対応原理 bitwise（T={T}）: {ok_b1}")

    # ---- B2: 保存則（結合力学） ----
    eng = VertexEngineV2(n, C0_of(Z0c / np.linalg.norm(Z0c), d=0.1), wp0, vertex_on=True)
    W0 = np.fft.ifft(eng.C, axis=1) * nreg
    cl0 = np.abs(np.sum(W0 ** 2, axis=0))
    sl0 = {k: complex(eng.C[:, k] @ eng.C[:, k]) for k in (1, 2)}
    n0 = float(np.sum(np.abs(W0) ** 2))
    for t in range(200):
        eng.step()
    W1 = np.fft.ifft(eng.C, axis=1) * nreg
    cl1 = np.abs(np.sum(W1 ** 2, axis=0))
    point_drift = float(np.max(np.abs(cl1 - cl0)))
    slice_drift = max(abs(complex(eng.C[:, k] @ eng.C[:, k]) - sl0[k]) for k in (1, 2))
    norm_drift = abs(float(np.sum(np.abs(W1) ** 2)) - n0) / n0
    out["B2"] = {"pointwise_closure_drift": point_drift,
                  "slice_closure_drift": float(slice_drift),
                  "norm_drift_rel": float(norm_drift)}
    out["B2_pass"] = bool(point_drift <= 1e-10 and slice_drift <= 1e-10
                           and norm_drift <= 1e-10)
    print(f"  B2 保存則: 点ごと閉塞={point_drift:.2e} スライス閉塞={slice_drift:.2e} "
          f"ノルム={norm_drift:.2e} → {out['B2_pass']}")

    # ---- B3: 物理頑健性 ----
    rows = []
    for d in (1e-3, 3e-3, 1e-2, 3e-2):
        eng = VertexEngineV2(n, C0_of(Z0c / np.linalg.norm(Z0c), d=d), wp0)
        fs = []
        for t in range(60):
            eng.step()
            fs.append(eng.diagnostics()["f_seed"])
        fs = np.array(fs)
        tt = np.arange(5, 40, dtype=float)
        A = np.vstack([tt, np.ones_like(tt)]).T
        coef, _, _, _ = np.linalg.lstsq(A, np.log(fs[5:40]), rcond=None)
        rows.append({"delta": d, "f0": float(fs[0]), "rate": float(coef[0])})
    ln_f = np.log([q["f0"] for q in rows])
    ln_r = np.log([q["rate"] for q in rows])
    A = np.vstack([ln_f, np.ones_like(ln_f)]).T
    coef, _, _, _ = np.linalg.lstsq(A, ln_r, rcond=None)
    p_exp, C_v2 = float(coef[0]), float(np.exp(coef[1]))
    d3 = {}
    for tag, pump in (("minus", Z0c / np.linalg.norm(Z0c)),
                       ("minus_conj", np.conj(Z0c) / np.linalg.norm(Z0c)),
                       ("plus", vplus), ("plus_conj", np.conj(vplus))):
        eng = VertexEngineV2(n, C0_of(pump), np.random.default_rng(96000).normal(size=m))
        fs = []
        for t in range(60):
            eng.step()
            fs.append(eng.diagnostics()["f_seed"])
        tt = np.arange(5, 40, dtype=float)
        A = np.vstack([tt, np.ones_like(tt)]).T
        coef, _, _, _ = np.linalg.lstsq(A, np.log(np.array(fs)[5:40]), rcond=None)
        d3[tag] = float(coef[0])
    branch_order = (min(d3["plus"], d3["minus_conj"])
                    > max(d3["minus"], d3["plus_conj"]))
    out["B3"] = {"p": p_exp, "C": C_v2, "d3_rates": d3,
                  "branch_order_preserved": bool(branch_order)}
    out["B3_pass"] = bool(1.9 <= p_exp <= 2.1 and 0.5 <= C_v2 / 6.378 <= 2.0)
    print(f"  B3 点火則: p={p_exp:.3f} C={C_v2:.3f}（v1: 2.001/6.378） → {out['B3_pass']}")
    print(f"     D3再検: {['%s=%.2e' % (k, r_) for k, r_ in d3.items()]} 枝順序保存={branch_order}")

    # ---- B4: ロック持続 ----
    def agg(x):
        S = np.zeros(n, complex)
        np.add.at(S, ia, x)
        np.add.at(S, ib, x)
        return S[ia] + S[ib] - 2 * x

    def full_pred(C, R):
        c1, c2 = C[:, 1], C[:, 2]
        beat = c2 * np.conj(c1)
        par = c2 ** 2
        return (R * (agg(beat) * c2 - agg(par) * np.conj(c1))
                + (agg(R * beat) * c2 - agg(R * par) * np.conj(c1)))

    cohs = []
    eng = VertexEngineV2(n, C0_of(Z0c / np.linalg.norm(Z0c)), wp0)
    for t in range(10):
        R = eng._readout()
        c3b = eng.C[:, 3].copy()
        pred = full_pred(eng.C, R)
        eng.step()
        dc3 = eng.C[:, 3] - c3b
        cohs.append(float(abs(np.vdot(pred, dc3))
                          / max(np.linalg.norm(pred) * np.linalg.norm(dc3), 1e-300)))
    out["B4"] = {"coherences": cohs, "final": cohs[-1]}
    out["B4_pass"] = bool(cohs[-1] > 0.9936)
    print(f"  B4 ロック持続: t=1..10 = {[round(c,4) for c in cohs]}（v1末尾0.9936） → {out['B4_pass']}")

    # ---- 毛(η)拡張 census P2 完全版 ----
    print("  === 毛(η)拡張（Nn=5, Nη=8） ===")
    Nn, Neta = 5, 8
    regflat = Nn * Neta

    def put(C, k, mm, vec):
        C[:, k * Neta + ((mm) % Neta)] = vec   # 平坦化: idx = k*Nη + m

    # 平坦化レジスタでは点ごと頂点が (n,η) 2D格子の点ごと積になるよう
    # 2D FFT を自前で行う: C2[k,m] per edge。
    class HairEngine(VertexEngineV2):
        def __init__(self, n_, C2_0, wp, **kw):
            # C2_0: M×Nn×Nη → 平坦 M×(Nn·Nη)
            C0f = C2_0.reshape(C2_0.shape[0], -1)
            super().__init__(n_, C0f, wp, **kw)
            self.Nn, self.Neta = Nn, Neta
            ks = np.arange(Nn)
            self.odd_k = (ks % 2 == 1)
            self.even_k = (ks % 2 == 0) & (ks != 0)

        def C2(self):
            return self.C.reshape(self.m, self.Nn, self.Neta)

        def _readout(self):
            P2 = np.abs(self.C2()) ** 2
            Pk = P2.sum(axis=2)                       # M×Nn（毛を周辺化）
            Av = np.zeros((self.n, self.Nn))
            np.add.at(Av, self.ia, Pk)
            np.add.at(Av, self.ib, Pk)
            Sagg = Av[self.ia] + Av[self.ib] - 2 * Pk
            comb = Pk + Sagg
            Pf = comb[:, self.odd_k].sum(axis=1)
            Pb = comb[:, self.even_k].sum(axis=1)
            th = np.arctan2(np.sqrt(np.maximum(Pf, 0)), np.sqrt(np.maximum(Pb, 0)))
            return self.scale * np.sin(th) ** 2

        def _nonlinear(self):
            R = self._readout()
            if not np.any(R > 0):
                return
            C2 = self.C2()
            W = np.fft.ifft2(C2, axes=(1, 2)) * (self.Nn * self.Neta)
            Wf = W.reshape(self.m, -1)
            rate0 = self._vertex_rate(Wf, R)
            Lmax = float(np.max(np.abs(rate0))) / max(float(np.max(np.abs(Wf))), 1e-300)
            nsub = max(1, int(np.ceil(Lmax / s2.H_MAX)))
            h = 1.0 / nsub
            for _ in range(nsub):
                k1 = self._vertex_rate(Wf, R)
                k2 = self._vertex_rate(Wf + 0.5 * h * k1, R)
                k3 = self._vertex_rate(Wf + 0.5 * h * k2, R)
                k4 = self._vertex_rate(Wf + h * k3, R)
                Wf = Wf + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
            W = Wf.reshape(self.m, self.Nn, self.Neta)
            self.C = (np.fft.fft2(W, axes=(1, 2)) / (self.Nn * self.Neta)
                      ).reshape(self.m, -1)

    C2_0 = np.zeros((m, Nn, Neta), complex)
    C2_0[:, 2, 2] = Z0c / np.linalg.norm(Z0c)      # ポンプ (k=2, m=+2)
    C2_0[:, 1, 1] = delta * seed_state              # 種 (k=1, m=+1)
    heng = HairEngine(n, C2_0, wp0)
    for t in range(30):
        heng.step()
    C2 = heng.C2()
    P_k3 = np.sum(np.abs(C2[:, 3, :]) ** 2, axis=0)    # 毛分解した k=3 パワー
    m_star = 3                                          # 予言 2·2−1=3
    excl = float(P_k3[m_star] / max(P_k3.sum() - P_k3[m_star], 1e-300))
    h1 = bool(excl >= 100)
    # H2: 2D完全予測子（毛込み）
    heng2 = HairEngine(n, C2_0, wp0)
    for t in range(1):
        R = heng2._readout()
        c1v = heng2.C2()[:, 1, 1]
        c2v = heng2.C2()[:, 2, 2]
        beat = c2v * np.conj(c1v)
        par = c2v ** 2
        pred = (R * (agg(beat) * c2v - agg(par) * np.conj(c1v))
                + (agg(R * beat) * c2v - agg(R * par) * np.conj(c1v)))
        before = heng2.C2()[:, 3, m_star].copy()
        heng2.step()
        dc = heng2.C2()[:, 3, m_star] - before
        coh_h = float(abs(np.vdot(pred, dc))
                      / max(np.linalg.norm(pred) * np.linalg.norm(dc), 1e-300))
    h2 = bool(coh_h >= 0.99)
    out["HAIR"] = {"P_k3_by_m": [float(x) for x in P_k3], "m_star": m_star,
                    "exclusivity": excl, "lock_coherence": coh_h,
                    "H1": h1, "H2": h2}
    out["HAIR_pass"] = bool(h1 and h2)
    print(f"  H1 毛の排他性: P(m=3)/P(m≠3)={excl:.2e} → {h1}")
    print(f"  H2 毛込みロック: コヒーレンス={coh_h:.4f} → {h2}")

    ok = all(out[k] for k in ("B1_pass", "B2_pass", "B3_pass", "B4_pass", "HAIR_pass"))
    out["all_pass"] = bool(ok)
    out["runtime_sec"] = time.time() - t0
    print(f"\n総合: {'ALL PASS' if ok else '不成立あり'}")
    (HERE / "stage3_sharedO_hair_result_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
