#!/usr/bin/env python3
"""段階2: 媒介非弾性頂点エンジン v1——単体検証＋多体点火実験（S2-E1）

アーキテクチャ（媒介頂点の一意化と多体アーキテクチャ_v1.md §3-4）:
    状態 = レジスタ行列 C ∈ C^{M×Nreg}（辺 × 倍音）。1stepは
    (a) 線形部: 占有スライスごとに独立 Cayley（段階1エンジン、bitwise検証済み）
    (b) 非弾性部: n空間で媒介頂点の流れを RK4 部分刻みで τ∈[0,1] 積分
        δw_e[n] = (i/2)[R_e(𝒜_e w_e − ℬ_e w̄_e) + (𝒜ᴿ_e w_e − ℬᴿ_e w̄_e)]
        （定理§3の一意形・算術平均強度・O(M)/点・恒等因数分解）
    強度 R_e = sin²θ_e、θ_e = atan2(√P_f, √P_b)。P_f/P_b は辺 e の
    レジスタパワースペクトルと頂点集約パワースペクトル（隣接辺の
    パワー和——GATE-0のγ教訓により振幅和でなくパワー和）の合成に
    奇数倍音（フェルミオン）/偶数倍音マスクを適用（step毎に読出し・
    substep中は固定——二体toyの衝突毎読出しと同文法）。
    DC/Nyquistは初期データから排除（定理）。頂点が生成するDC内容は
    投影せず観測量として記録する（存在資格なし内容の生成率＝診断）。

単体検証（実行前固定）:
    U1 対応原理: 頂点OFF（strength倍率0）で段階1エンジンと bitwise 一致（T=200）。
    U2 保存則: 頂点ONで、各レジスタ点の場の閉塞 |Σ_e w_e[n]²| と
       全ノルム Σ|w|² のドリフトが T=200 で ≤1e-10（RK4・h_max=0.02。
       二体v3の 2.9e-3/3000衝突 からの改善目標）。
    U3 生成の存在と唯一性: 種因子 f_seed=P_odd/P_tot が、頂点OFFでは
       一定（≤1e-14、段階1のA4）、頂点ONでは有意に変化（≥1e-6）。

S2-E1 多体点火実験（実行前固定の予言）:
    ポンプ: 偶スライス k=2 に control 親（コヒーレント凝縮体）。
    種: 奇スライス k=1 に白色セクター状態 ×δ、δ ∈ {1e-3,3e-3,1e-2,3e-2}。
    P-E1a（点火則）: 初期成長率 rate = d(ln P_odd)/d(step) は f_seed に
       冪 p でスケール。二体は p=2（rate=C·f²）。多体の p を実測し、
       C_multi の一定性（δ掃引で±30%以内なら「一定」）を判定。
    P-E1b（ゲート因子）: 同一ノルム・同一δで、ポンプを非凝縮状態
       （一般零閉塞ランダム状態）に替えると初期成長率が低下する
       （凝縮/非凝縮 比 >1）。比の大きさは探索的記録。

使い方: python3 run_stage2_vertex_engine_v1.py
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SERIES = HERE.parent
ABL = SERIES / "第8論文_二段階seed除去による準安定相の因果分離" / "code" / "run_preliminary_seed_ablation_v1.py"
GEN3 = SERIES / "make_parent_white_managed_v1" / "make_parent_white_harmonics_n_only_v3.py"

spec = importlib.util.spec_from_file_location("abl_s2", ABL)
abl = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = abl
spec.loader.exec_module(abl)
spec2 = importlib.util.spec_from_file_location("gen3_s2", GEN3)
gen3 = importlib.util.module_from_spec(spec2)
sys.modules[spec2.name] = gen3
spec2.loader.exec_module(gen3)

H_MAX = 0.02


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


class VertexEngine:
    """線形スライス部＋媒介非弾性頂点。"""

    def __init__(self, n, C0, wps, vertex_on=True, strength_scale=1.0):
        self.n = n
        self.nreg = C0.shape[1]
        self.m = C0.shape[0]
        self.C = C0.copy()
        self.wps = {k: w.copy() for k, w in wps.items()}
        self.vertex_on = vertex_on
        self.scale = strength_scale
        self.ia, self.ib = np.triu_indices(n, k=1)
        ks = np.arange(self.nreg)
        self.odd_mask = (ks % 2 == 1)
        self.even_mask = (ks % 2 == 0) & (ks != 0)
        if self.nreg % 2 == 0:
            self.even_mask[self.nreg // 2] = False
        self.sys = {}

    def _occupied(self):
        return [k for k in range(self.nreg) if np.linalg.norm(self.C[:, k]) > 0.0]

    def _linear(self):
        for k in self._occupied():
            if k not in self.sys:
                self.sys[k] = abl.LowRankSystem(self.n)
                self.sys[k].set_theta(np.angle(self.C[:, k]))
            if k not in self.wps:
                self.wps[k] = np.random.default_rng(92000 + k).normal(size=self.m)
            Znew, self.wps[k] = abl.evolve(self.sys[k], self.C[:, k], self.wps[k])
            self.C[:, k] = Znew

    def _readout(self):
        P = np.abs(self.C) ** 2                      # M×Nreg 辺別パワースペクトル
        Av = np.zeros((self.n, self.nreg))
        np.add.at(Av, self.ia, P)
        np.add.at(Av, self.ib, P)
        Sagg = Av[self.ia] + Av[self.ib] - 2 * P      # 隣接辺パワー集約（自己控除）
        comb = P + Sagg
        Pf = comb[:, self.odd_mask].sum(axis=1)
        Pb = comb[:, self.even_mask].sum(axis=1)
        theta = np.arctan2(np.sqrt(np.maximum(Pf, 0.0)), np.sqrt(np.maximum(Pb, 0.0)))
        return self.scale * np.sin(theta) ** 2

    def _vertex_rate(self, W, R):
        a2 = np.abs(W) ** 2
        z2 = W ** 2
        A = np.zeros((self.n, self.nreg))
        B = np.zeros((self.n, self.nreg), complex)
        AR = np.zeros((self.n, self.nreg))
        BR = np.zeros((self.n, self.nreg), complex)
        Ra2 = R[:, None] * a2
        Rz2 = R[:, None] * z2
        np.add.at(A, self.ia, a2); np.add.at(A, self.ib, a2)
        np.add.at(B, self.ia, z2); np.add.at(B, self.ib, z2)
        np.add.at(AR, self.ia, Ra2); np.add.at(AR, self.ib, Ra2)
        np.add.at(BR, self.ia, Rz2); np.add.at(BR, self.ib, Rz2)
        cA = A[self.ia] + A[self.ib] - 2 * a2
        cB = B[self.ia] + B[self.ib] - 2 * z2
        cAR = AR[self.ia] + AR[self.ib] - 2 * Ra2
        cBR = BR[self.ia] + BR[self.ib] - 2 * Rz2
        Wc = np.conj(W)
        return 0.5j * (R[:, None] * (cA * W - cB * Wc) + (cAR * W - cBR * Wc))

    def _nonlinear(self):
        R = self._readout()
        if not np.any(R > 0):
            return
        W = np.fft.ifft(self.C, axis=1) * self.nreg
        rate0 = self._vertex_rate(W, R)
        Lmax = float(np.max(np.abs(rate0))) / max(float(np.max(np.abs(W))), 1e-300)
        nsub = max(1, int(np.ceil(Lmax / H_MAX)))
        h = 1.0 / nsub
        for _ in range(nsub):
            k1 = self._vertex_rate(W, R)
            k2 = self._vertex_rate(W + 0.5 * h * k1, R)
            k3 = self._vertex_rate(W + 0.5 * h * k2, R)
            k4 = self._vertex_rate(W + h * k3, R)
            W = W + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        self.C = np.fft.fft(W, axis=1) / self.nreg

    def step(self):
        self._linear()
        if self.vertex_on:
            self._nonlinear()

    def diagnostics(self):
        W = np.fft.ifft(self.C, axis=1) * self.nreg
        closure = np.abs(np.sum(W ** 2, axis=0))          # レジスタ点ごとの場の閉塞
        norm = float(np.sum(np.abs(W) ** 2))
        P = np.abs(self.C) ** 2
        P_odd = float(P[:, self.odd_mask].sum())
        P_tot = float(P.sum())
        P_dc = float(P[:, 0].sum())
        return {"closure_max": float(closure.max()), "norm": norm,
                "f_seed": P_odd / P_tot if P_tot > 0 else 0.0, "P_dc": P_dc}


def main() -> None:
    t0 = time.time()
    n = 5
    m = 10
    nreg = 5
    print(f"段階2 頂点エンジン N={n} Nreg={nreg}  ABL {sha256(ABL)[:16]}…")
    out = {"imports": {"abl": sha256(ABL), "generator_v3": sha256(GEN3)},
           "criteria": {"H_MAX": H_MAX,
                         "U2_drift": "<=1e-10 (T=200)",
                         "P_E1a": "rate ~ C*f_seed^p; C constant within ±30% over delta sweep",
                         "P_E1b": "condensed pump rate > incoherent pump rate"}}

    _, v, _, _, _, _, _, Z0c, wp0 = abl.build_init(n, False)
    r = gen3.make_parent(n, seed=2)
    Csec = np.fft.fft(r.relation_waves, axis=1) / n

    # ---- U1: 頂点OFF＝段階1とbitwise ----
    C0 = np.zeros((m, nreg), complex)
    C0[:, 2] = Z0c
    C0[:, 1] = 0.01 * Csec[:, 1] / np.linalg.norm(Csec[:, 1])
    wps = {1: np.random.default_rng(92001).normal(size=m),
           2: wp0.copy()}
    e_off = VertexEngine(n, C0, wps, vertex_on=False)
    e_ref = VertexEngine(n, C0, wps, vertex_on=True, strength_scale=0.0)
    ok_u1 = True
    for t in range(200):
        e_off.step()
        e_ref.step()
        if not np.array_equal(e_off.C, e_ref.C):
            ok_u1 = False
            break
    out["U1_pass"] = bool(ok_u1)
    print(f"  U1 頂点OFF bitwise（scale=0 vs vertex_on=False, T=200）: {ok_u1}")

    # ---- U2+U3: 保存則と生成の存在 ----
    eng = VertexEngine(n, C0, wps, vertex_on=True)
    d0 = eng.diagnostics()
    fs = [d0["f_seed"]]
    for t in range(200):
        eng.step()
    d1 = eng.diagnostics()
    closure_drift = abs(d1["closure_max"] - d0["closure_max"])
    norm_drift = abs(d1["norm"] - d0["norm"]) / d0["norm"]
    dfseed_on = abs(d1["f_seed"] - d0["f_seed"])
    lin = VertexEngine(n, C0, wps, vertex_on=False)
    l0 = lin.diagnostics()
    for t in range(200):
        lin.step()
    l1 = lin.diagnostics()
    dfseed_off = abs(l1["f_seed"] - l0["f_seed"])
    out["U2"] = {"closure_drift": closure_drift, "norm_drift_rel": norm_drift,
                  "P_dc_final": d1["P_dc"]}
    out["U2_pass"] = bool(closure_drift <= 1e-10 and norm_drift <= 1e-10)
    out["U3"] = {"dfseed_vertex_on": dfseed_on, "dfseed_vertex_off": dfseed_off}
    out["U3_pass"] = bool(dfseed_on >= 1e-6 and dfseed_off <= 1e-14)
    print(f"  U2 保存則（T=200）: 閉塞ドリフト={closure_drift:.2e} ノルム相対={norm_drift:.2e} "
          f"DC生成={d1['P_dc']:.2e} → {out['U2_pass']}")
    print(f"  U3 生成: Δf_seed ON={dfseed_on:.2e} OFF={dfseed_off:.2e} → {out['U3_pass']}")

    # ---- S2-E1a: 点火則（δ掃引） ----
    print("  S2-E1a 点火則（凝縮ポンプ k=2 = control親、種 k=1）")
    rows = []
    seed_state = Csec[:, 1] / np.linalg.norm(Csec[:, 1])
    for delta in (1e-3, 3e-3, 1e-2, 3e-2):
        C0 = np.zeros((m, nreg), complex)
        C0[:, 2] = Z0c
        C0[:, 1] = delta * seed_state
        eng = VertexEngine(n, C0, {1: np.random.default_rng(92001).normal(size=m),
                                    2: wp0.copy()})
        po = []
        for t in range(60):
            eng.step()
            d = eng.diagnostics()
            po.append(d["f_seed"])
        po = np.array(po)
        f0 = po[0]
        win = np.log(po[5:40])
        tt = np.arange(5, 40, dtype=float)
        A = np.vstack([tt, np.ones_like(tt)]).T
        coef, _, _, _ = np.linalg.lstsq(A, win, rcond=None)
        rate = float(coef[0])
        rows.append({"delta": delta, "f_seed0": float(f0), "rate": rate})
        print(f"    δ={delta:.0e}: f_seed(0)={f0:.2e} rate={rate:.3e}/step")
    ln_f = np.log([r_["f_seed0"] for r_ in rows])
    ln_rate = np.log([max(r_["rate"], 1e-300) for r_ in rows])
    A = np.vstack([ln_f, np.ones_like(ln_f)]).T
    coef, _, _, _ = np.linalg.lstsq(A, ln_rate, rcond=None)
    p_exp = float(coef[0])
    Cs = [r_["rate"] / (r_["f_seed0"] ** p_exp) for r_ in rows]
    c_spread = (max(Cs) - min(Cs)) / np.mean(Cs)
    out["E1a"] = {"rows": rows, "power_p": p_exp, "C_values": Cs,
                   "C_spread": float(c_spread)}
    out["E1a_pass"] = bool(c_spread < 0.3)
    print(f"    冪 p={p_exp:.3f}（二体はp=2相当）C一定性: 広がり{c_spread*100:.1f}% → {out['E1a_pass']}")

    # ---- S2-E1b: ゲート因子（凝縮 vs 非凝縮ポンプ） ----
    delta = 1e-2
    rng = np.random.default_rng(94000)
    A2 = rng.normal(size=(m, 2))
    Q, _ = np.linalg.qr(A2)
    Zinc = (Q[:, 0] + 1j * Q[:, 1]) / np.sqrt(2)   # 一般零閉塞（非凝縮）
    rates = {}
    for tag, pump in (("condensed", Z0c), ("incoherent", Zinc)):
        C0 = np.zeros((m, nreg), complex)
        C0[:, 2] = pump / np.linalg.norm(pump)
        C0[:, 1] = delta * seed_state
        eng = VertexEngine(n, C0, {1: np.random.default_rng(92001).normal(size=m),
                                    2: np.random.default_rng(96000).normal(size=m)})
        po = []
        for t in range(60):
            eng.step()
            po.append(eng.diagnostics()["f_seed"])
        po = np.array(po)
        tt = np.arange(5, 40, dtype=float)
        A = np.vstack([tt, np.ones_like(tt)]).T
        coef, _, _, _ = np.linalg.lstsq(A, np.log(po[5:40]), rcond=None)
        rates[tag] = float(coef[0])
        print(f"  S2-E1b {tag}: rate={rates[tag]:.3e}/step")
    ratio = rates["condensed"] / rates["incoherent"] if rates["incoherent"] > 0 else np.inf
    out["E1b"] = {"rates": rates, "ratio": float(ratio) if np.isfinite(ratio) else None}
    out["E1b_pass"] = bool(rates["condensed"] > rates["incoherent"])
    print(f"    凝縮/非凝縮 比 = {ratio if np.isfinite(ratio) else 'inf'} → {out['E1b_pass']}")

    ok = all(out[k] for k in ("U1_pass", "U2_pass", "U3_pass", "E1a_pass", "E1b_pass"))
    out["all_pass"] = bool(ok)
    out["runtime_sec"] = time.time() - t0
    print(f"\n総合: {'ALL PASS' if ok else '不成立あり——反証記録として保存'}")
    (HERE / "stage2_vertex_engine_result_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
