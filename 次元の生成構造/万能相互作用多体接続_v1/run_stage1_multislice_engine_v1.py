#!/usr/bin/env python3
"""段階1: マルチスライス・レジスタ行列エンジン v1——受入基準4項の検証

アーキテクチャ（媒介頂点の一意化と多体アーキテクチャ_v1.md §4）:
    状態 = レジスタ行列 C ∈ C^{M×N_slice}（辺 e × スライス k）。
    線形部 = スライスごとに独立に既存 abl.evolve（Cayley一歩）を適用。
    零スライスは厳密に不動（線形性）なのでスキップは恒等（数値ガード）。
    スライスを混ぜられる演算は存在しない——線形無生成が設計定理。

受入基準（実行前固定・媒介頂点md §5）:
    A1 対応原理: 単一スライス初期条件で既存 abl 軌道と bitwise 一致
       （N=5: T=1000 / N=40: T=300）。
    A2 非結合: 多スライス同時走行の各スライスが、単独走行と bitwise 一致
       （N=5、白色セクター4本同時、T=2000）。
    A3 二分再現: 混載レジスタ（スライス0=control親、1..3=白色セクター）で
       control=潜伏バースト（crossing≈1166）／セクター=即時（11〜16）。
    A4 種因子計装: レジスタ奇偶比 f_seed=P_odd/(P_odd+P_even) が全軌道で
       厳密一定（|Δf_seed| ≤ 1e-14）——線形無生成の直接実測。

使い方: python3 run_stage1_multislice_engine_v1.py
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

spec = importlib.util.spec_from_file_location("abl_s1", ABL)
abl = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = abl
spec.loader.exec_module(abl)
spec2 = importlib.util.spec_from_file_location("gen3_s1", GEN3)
gen3 = importlib.util.module_from_spec(spec2)
sys.modules[spec2.name] = gen3
spec2.loader.exec_module(gen3)


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


class MultiSliceEngine:
    """レジスタ行列の線形部。スライスごとに独立な abl.evolve。"""

    def __init__(self, n, C0, wps):
        self.n = n
        self.C = C0.copy()
        self.wps = {k: w.copy() for k, w in wps.items()}
        self.occupied = [k for k in range(C0.shape[1])
                         if np.linalg.norm(C0[:, k]) > 0.0]
        self.sys = {k: abl.LowRankSystem(n) for k in self.occupied}
        for k in self.occupied:
            self.sys[k].set_theta(np.angle(self.C[:, k]))

    def step(self):
        for k in self.occupied:
            Znew, self.wps[k] = abl.evolve(self.sys[k], self.C[:, k], self.wps[k])
            self.C[:, k] = Znew


def reference_run(n, Z0, wp, T):
    sys_lr = abl.LowRankSystem(n)
    sys_lr.set_theta(np.angle(Z0))
    Z = Z0.copy()
    w = wp.copy()
    traj = [Z.copy()]
    for _ in range(T):
        Z, w = abl.evolve(sys_lr, Z, w)
        traj.append(Z.copy())
    return traj


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
    print(f"段階1 マルチスライスエンジン  ABL {sha256(ABL)[:16]}…  GEN3 {sha256(GEN3)[:16]}…")
    out = {"imports": {"abl": sha256(ABL), "generator_v3": sha256(GEN3)}}

    # ---- A1: 対応原理（bitwise） ----
    a1 = {}
    for n, T in ((5, 1000), (40, 300)):
        m = n * (n - 1) // 2
        _, v, _, _, _, _, _, Z0c, wp0 = abl.build_init(n, False)
        ref = reference_run(n, Z0c, wp0, T)
        C0 = np.zeros((m, 8), complex)
        C0[:, 1] = Z0c
        eng = MultiSliceEngine(n, C0, {1: wp0})
        exact = True
        for t in range(1, T + 1):
            eng.step()
            if not np.array_equal(eng.C[:, 1], ref[t]):
                exact = False
                break
        a1[f"N{n}"] = {"T": T, "bitwise": bool(exact)}
        print(f"  A1 N={n} T={T}: bitwise一致 = {exact}")
    out["A1"] = a1
    out["A1_pass"] = all(d["bitwise"] for d in a1.values())

    # ---- A2: 非結合（bitwise・N=5 セクター4本同時） ----
    n, T = 5, 2000
    m = 10
    r = gen3.make_parent(n, seed=2)
    Csec = np.fft.fft(r.relation_waves, axis=1) / n
    ks = [1, 2, 3, 4]
    refs = {}
    for k in ks:
        Zk = Csec[:, k] / np.linalg.norm(Csec[:, k])
        refs[k] = reference_run(n, Zk, np.random.default_rng(92000 + k).normal(size=m), T)
    C0 = np.zeros((m, 8), complex)
    wps = {}
    for k in ks:
        C0[:, k] = Csec[:, k] / np.linalg.norm(Csec[:, k])
        wps[k] = np.random.default_rng(92000 + k).normal(size=m)
    eng = MultiSliceEngine(n, C0, wps)
    exact = {k: True for k in ks}
    for t in range(1, T + 1):
        eng.step()
        for k in ks:
            if exact[k] and not np.array_equal(eng.C[:, k], refs[k][t]):
                exact[k] = False
    out["A2"] = {f"k{k}": bool(v) for k, v in exact.items()}
    out["A2_pass"] = all(exact.values())
    print(f"  A2 N=5 セクター4本同時 T={T}: 各スライスbitwise一致 = {out['A2']}")

    # ---- A3: 二分再現＋ A4: 種因子計装（混載レジスタ） ----
    n, T = 5, 2000
    _, v, _, _, _, _, _, Z0c, wp0 = abl.build_init(n, False)
    C0 = np.zeros((m, 8), complex)
    C0[:, 0] = Z0c                       # スライス0（偶アドレス）= control親
    wps = {0: wp0}
    for k in (1, 2, 3):
        C0[:, k] = Csec[:, k] / np.linalg.norm(Csec[:, k])
        wps[k] = np.random.default_rng(92000 + k).normal(size=m)
    eng = MultiSliceEngine(n, C0, wps)
    fvals = {k: fval_factory(C0[:, k]) for k in (0, 1, 2, 3)}
    crossings = {k: None for k in (0, 1, 2, 3)}
    P_odd0 = sum(np.linalg.norm(C0[:, k]) ** 2 for k in range(8) if k % 2 == 1)
    P_tot0 = sum(np.linalg.norm(C0[:, k]) ** 2 for k in range(8))
    fseed0 = P_odd0 / P_tot0
    max_dev = 0.0
    for t in range(1, T + 1):
        eng.step()
        for k in crossings:
            if crossings[k] is None and fvals[k](eng.C[:, k]) > 0.05:
                crossings[k] = t
        P_odd = sum(np.linalg.norm(eng.C[:, k]) ** 2 for k in range(8) if k % 2 == 1)
        P_tot = sum(np.linalg.norm(eng.C[:, k]) ** 2 for k in range(8))
        max_dev = max(max_dev, abs(P_odd / P_tot - fseed0))
    out["A3"] = {"crossings": {str(k): c for k, c in crossings.items()}}
    a3_pass = (crossings[0] is not None and crossings[0] > 1000
               and all(crossings[k] is not None and crossings[k] < 30 for k in (1, 2, 3)))
    out["A3_pass"] = bool(a3_pass)
    print(f"  A3 混載: crossings = {crossings}（予言: スライス0≈1166潜伏バースト、1-3は11〜16即時）→ {a3_pass}")
    out["A4"] = {"fseed0": float(fseed0), "max_deviation": float(max_dev)}
    out["A4_pass"] = bool(max_dev <= 1e-14)
    print(f"  A4 種因子: f_seed(0)={fseed0:.6f} 最大偏差={max_dev:.2e}（≤1e-14で線形無生成の実測）→ {out['A4_pass']}")

    ok = out["A1_pass"] and out["A2_pass"] and out["A3_pass"] and out["A4_pass"]
    out["all_pass"] = bool(ok)
    out["runtime_sec"] = time.time() - t0
    print(f"\n受入判定: {'ALL PASS——段階1完了' if ok else 'FAIL あり'}")
    (HERE / "stage1_acceptance_result_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
