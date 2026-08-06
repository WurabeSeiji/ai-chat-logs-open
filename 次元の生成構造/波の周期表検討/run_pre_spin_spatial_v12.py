#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v12: スピン量子数の空間側計器——凝縮体軸まわりの回転荷 ℓ の実測

原理（事前記録）: v11でチャネル側にスピン統計が不在（厳密対称）と確定。
空間側: 凝縮体の回転軸生成子 K で状態を物理回転 R(β)=exp(βK̂) させ、
各SVD平面読出し c_j(β) の位相前進から回転荷 ℓ_j = dφ_j/dβ を測る。
ℓ=±1: ベクトル（2π回帰）／ℓ=±2: テンソル・スピン2的（π回帰）／
半整数: スピノル・4π回帰（124/248二重被覆の空間実現）。

判定（事前固定）:
  各平面の ℓ_j を線形フィット（β∈[0,4π], 33点）し、最近整数/半整数との偏差と
  回帰忠実度 F(2π), F(4π) を記録。
  H_half: いずれかの平面で ℓ が半整数（|ℓ−(k+1/2)|<0.05）かつ F(2π)<0.5<F(4π)
  → 空間スピノル発見。H_int: 全て整数 → この凝縮体の観測モードは
  ベクトル/テンソル世界（スピノルは読出し枠に現れない）を記録。
  N=5,6,8 で実施。

使い方: python3 run_pre_spin_spatial_v12.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
SPACE = HERE.parent / "空間軸3軸と固有時間の創生_v1"
spec1 = importlib.util.spec_from_file_location("pre1_v12", SPACE / "run_pre_2plus1_structure_v1.py")
pre1 = importlib.util.module_from_spec(spec1)
sys.modules[spec1.name] = pre1
spec1.loader.exec_module(pre1)
abl = pre1.abl
edge_adjacency, build_K = pre1.edge_adjacency, pre1.build_K

T_END = 4000; WIN = (2000, 4000); SAMPLE_EVERY = 5
N_PLANES = 5

def analyze(n):
    sys_lr, v, _, _, _, p, q, Z, wp = abl.build_init(n, True)
    M = sys_lr.m
    adj = edge_adjacency(n)
    samples = []
    for t in range(T_END):
        Z, wp = abl.evolve(sys_lr, Z, wp)
        if WIN[0] <= t < WIN[1] and (t % SAMPLE_EVERY == 0):
            samples.append(Z.copy())
    S = np.array(samples)
    Sp = S - np.outer(S @ p, p) - np.outer(S @ q, q)
    X = np.hstack([Sp.real, Sp.imag])
    Xc = X - X.mean(axis=0)
    _, sv, Vt = np.linalg.svd(Xc, full_matrices=False)

    Zstar = S[-1].copy()
    theta = np.angle(Zstar)
    K = build_K(theta, adj)
    K20 = np.zeros((2 * M, 2 * M))
    K20[:M, :M] = K; K20[M:, M:] = K
    A = 0.5 * (K20 - K20.T)
    # 生成子の規格化: 上位平面の回転レートで割り、β=平面1の物理回転角にする
    U3 = Vt[:3].T
    A3 = U3.T @ A @ U3; A3 = 0.5 * (A3 - A3.T)
    w_vec = np.array([A3[2, 1], A3[0, 2], A3[1, 0]])
    omega = float(np.linalg.norm(w_vec))
    Ahat = A / max(omega, 1e-300)

    x0 = np.concatenate([Zstar.real, Zstar.imag])
    betas = np.linspace(0, 4 * np.pi, 33)
    R_step = expm(Ahat * (betas[1] - betas[0]))
    rows = []
    phases = {j: [] for j in range(N_PLANES)}
    overlaps = {j: [] for j in range(N_PLANES)}
    x = x0.copy()
    dcs = []
    for j in range(N_PLANES):
        d1, d2 = Vt[2 * j], Vt[2 * j + 1]
        dcs.append(d1 + 1j * d2)
    c0 = [np.dot(np.conj(dc), x0.astype(complex)) for dc in dcs]
    for bi, b in enumerate(betas):
        if bi > 0:
            x = R_step @ x
        for j in range(N_PLANES):
            c = np.dot(np.conj(dcs[j]), x.astype(complex))
            phases[j].append(np.angle(c))
            overlaps[j].append(abs(c) / max(abs(c0[j]), 1e-300))
    print(f"N={n} (ω={omega:.3e}):")
    for j in range(N_PLANES):
        ph = np.unwrap(np.array(phases[j]))
        ell = float(np.polyfit(betas, ph, 1)[0])
        resid = float(np.max(np.abs(ph - (ph[0] + ell * betas))))
        F2 = float(overlaps[j][16] * np.cos(ph[16] - ph[0] - 0))  # 2πでの複素回帰実部
        # 回帰忠実度: |⟨c(β), c(0)⟩|/|c0|² 相当 → 位相込み
        f2pi = float(np.cos(ph[16] - ph[0]) * overlaps[j][16])
        f4pi = float(np.cos(ph[32] - ph[0]) * overlaps[j][32])
        near_int = round(ell)
        near_half = round(ell * 2) / 2
        is_half = bool(abs(ell - near_half) < 0.05 and abs(near_half % 1) == 0.5)
        print(f"  平面{j+1}: ℓ={ell:+.4f}（最近整数{near_int:+d} 偏差{ell-near_int:+.4f}・"
              f"直線性残差{resid:.3f}） F(2π)={f2pi:+.3f} F(4π)={f4pi:+.3f}"
              f"{'  ←半整数!' if is_half else ''}")
        rows.append({"plane": j + 1, "ell": ell, "nearest_int": near_int,
                      "dev_int": ell - near_int, "linearity_resid": resid,
                      "F_2pi": f2pi, "F_4pi": f4pi, "is_half": is_half})
    return {"N": n, "omega": omega, "planes": rows}

def main():
    t0 = time.time()
    out = {"T_END": T_END, "WIN": list(WIN), "betas_n": 33, "scan": []}
    for n in (5, 6, 8):
        out["scan"].append(analyze(n))
    halves = [pl for r in out["scan"] for pl in r["planes"] if pl["is_half"]]
    ints = all(abs(pl["dev_int"]) < 0.1 for r in out["scan"] for pl in r["planes"]
               if pl["linearity_resid"] < 0.5)
    print(f"\nH_half（空間スピノル）= {bool(halves)}   "
          f"H_int（直線性良好な平面は全て整数ℓ）= {ints}")
    out["H_half"] = bool(halves); out["H_int"] = bool(ints)
    out["runtime_sec"] = time.time() - t0
    (HERE / "pre_spin_spatial_result_v12.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s → pre_spin_spatial_result_v12.json")

if __name__ == "__main__":
    main()
