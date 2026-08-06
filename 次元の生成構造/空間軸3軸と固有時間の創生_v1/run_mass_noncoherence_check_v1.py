#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""質量読出し(i)＝非コヒーレンスの機械検証（本稿初出導出の検収）

導出（本文§5に書き下ろす内容の要約）:
    二成分の複素読出し u=(c₁,c₂) の窓平均 Gram 行列 Γ=⟨u u†⟩ は 2×2 エルミート。
    Pauli 分解 Γ = T·I + X σx + Y σy + Z σz により
        det Γ = T² − X² − Y² − Z².
    Γ は半正定値ゆえ det Γ ≥ 0 が自動成立（光錐束縛）。
    det Γ = 0 ⟺ rank 1 ⟺ 窓内で u が単一複素振幅に比例（完全コヒーレント＝光的）。
    det Γ > 0 は読出し対の非コヒーレンスであり、これを質量²の読出しとする。

判定（実行前固定）:
    MC1 合成 rank-1 状態: detΓ/T² < 1e-12（機械零）。
    MC2 合成 非コヒーレント状態（独立2成分）: detΓ/T² が O(1)。
    MC3 模型実測: 親固有モード軌道（凝縮前の光的状態）の detΓ/T² が
        準安定凝縮体のそれより小さい（コヒーレント→質量小、凝縮→質量あり）。
        値は記述的に記録（探索）。

使い方: python3 run_mass_noncoherence_check_v1.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec1 = importlib.util.spec_from_file_location("pre1mc", HERE / "run_pre_2plus1_structure_v1.py")
pre1 = importlib.util.module_from_spec(spec1)
sys.modules[spec1.name] = pre1
spec1.loader.exec_module(pre1)
abl = pre1.abl


def gram_mass2(c1, c2):
    """窓平均Gram Γ=⟨u u†⟩ の Pauli 分解と detΓ（=T²−X²−Y²−Z²）。"""
    G = np.array([[np.mean(c1 * np.conj(c1)), np.mean(c1 * np.conj(c2))],
                  [np.mean(c2 * np.conj(c1)), np.mean(c2 * np.conj(c2))]])
    T = float(np.real(G[0, 0] + G[1, 1]) / 2)
    X = float(np.real(G[0, 1]))
    Y = float(np.imag(G[1, 0]))
    Z = float(np.real(G[0, 0] - G[1, 1]) / 2)
    det = float(np.real(np.linalg.det(G)))
    ident = det - (T ** 2 - X ** 2 - Y ** 2 - Z ** 2)
    return {"T": T, "X": X, "Y": Y, "Z": Z, "det": det,
            "mass2_over_T2": det / max(T ** 2, 1e-300),
            "pauli_identity_residual": ident}


def main() -> None:
    t0 = time.time()
    rng = np.random.default_rng(7)
    ns = 2000
    out = {}

    # ---- MC1 rank-1（光的）: u(t) = a(t)·(v1,v2) ----
    a = rng.normal(size=ns) + 1j * rng.normal(size=ns)
    v1c, v2c = 0.8 + 0.3j, -0.2 + 0.9j
    r1 = gram_mass2(a * v1c, a * v2c)
    print(f"[MC1] rank-1合成: 質量²/T² = {r1['mass2_over_T2']:.2e}"
          f"（Pauli恒等式残差 {r1['pauli_identity_residual']:.1e}）")
    mc1 = abs(r1["mass2_over_T2"]) < 1e-12
    out["MC1"] = {**r1, "pass": bool(mc1)}

    # ---- MC2 非コヒーレント合成: 独立2成分 ----
    b = rng.normal(size=ns) + 1j * rng.normal(size=ns)
    r2 = gram_mass2(a, b)
    print(f"[MC2] 非コヒーレント合成: 質量²/T² = {r2['mass2_over_T2']:.3f}（O(1)期待）")
    mc2 = r2["mass2_over_T2"] > 0.1
    out["MC2"] = {**r2, "pass": bool(mc2)}

    # ---- MC3 模型実測: 親固有モード（光的） vs 準安定凝縮体 ----
    N = 8
    sys_lr, v, _, _, _, p, q, Z, wp = abl.build_init(N, True)
    M = sys_lr.m
    # (a) 凝縮前の親固有モード軌道（seed OFF で純回転: 潜伏期の最初 400 step）
    sys2, v2_, _, _, _, p2, q2, Z2, wp2 = abl.build_init(N, False)
    S_par = []
    for t in range(400):
        Z2, wp2 = abl.evolve(sys2, Z2, wp2)
        S_par.append(Z2.copy())
    S_par = np.array(S_par)
    # (b) 準安定凝縮体
    S_meta = []
    for t in range(12000):
        Z, wp = abl.evolve(sys_lr, Z, wp)
        if t >= 6000 and t % 5 == 0:
            S_meta.append(Z.copy())
    S_meta = np.array(S_meta)

    def two_channel(S, p_, q_):
        """読出し2成分: 親平面複素座標 c1 と、親平面外の主方向複素座標 c2。"""
        c1 = (S @ p_) + 1j * (S @ q_)
        Sp = np.array([z - p_ * (p_ @ z) - q_ * (q_ @ z) for z in S])
        X2 = np.hstack([Sp.real, Sp.imag])
        Xc = X2 - X2.mean(axis=0)
        _, _, Vt = np.linalg.svd(Xc, full_matrices=False)
        u = Vt[0]
        Mloc = S.shape[1]
        d = u[:Mloc] + 1j * u[Mloc:]
        d /= np.linalg.norm(d)
        c2 = Sp @ np.conj(d)
        return c1, c2

    c1p, c2p = two_channel(S_par, p2, q2)
    c1m, c2m = two_channel(S_meta, p, q)
    rp = gram_mass2(c1p, c2p)
    rm = gram_mass2(c1m, c2m)
    print(f"[MC3] 親固有モード（光的候補）: 質量²/T² = {rp['mass2_over_T2']:.3e}")
    print(f"      準安定凝縮体:            質量²/T² = {rm['mass2_over_T2']:.3e}")
    mc3 = rp["mass2_over_T2"] < rm["mass2_over_T2"]
    print(f"      光的 < 凝縮体 = {mc3}（比 {rm['mass2_over_T2']/max(rp['mass2_over_T2'],1e-300):.1e}）")
    out["MC3"] = {"parent_lightlike": rp, "metastable": rm, "pass": bool(mc3)}

    out["all_pass"] = bool(mc1 and mc2 and mc3)
    out["runtime_sec"] = time.time() - t0
    print(f"判定: {'ALL PASS' if out['all_pass'] else '不成立あり'}")
    (HERE / "mass_noncoherence_check_result_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
