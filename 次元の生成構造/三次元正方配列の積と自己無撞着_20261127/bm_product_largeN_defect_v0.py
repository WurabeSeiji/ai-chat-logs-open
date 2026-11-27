#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
大 N (N = 12, 16, 24 の結晶解 T* および N = 11, 13 の反対称波周期結晶) に対する
純粋 BM 積の 1 マス欠陥発展: N を大きくしても全滅 (0 に消滅) するかの確認。

乱数なし。3 次斉次スケーリング分離を用い、||P|| < 1e-12 を厳密ゼロとして停止する。
"""
import math
import os
import time
import numpy as np


def bm_prod(T):
    return np.einsum("itk,ijt,tjk->ijk", T, T, T)


def proj_dist(A, B):
    na = np.linalg.norm(A)
    nb = np.linalg.norm(B)
    if na == 0 or nb == 0:
        return 1.0
    ov = abs(np.vdot(A / na, B / nb))
    return float(np.sqrt(max(0.0, 1.0 - ov * ov)))


def sum_z2_ratio(T):
    denom = float(np.sum(np.abs(T) ** 2))
    return 0.0 if denom == 0 else float(abs(np.sum(T ** 2)) / denom)


def build_sn_crystal_sz2_zero(N):
    """Step 2 の厳密結晶解 T* (a=0, b1=0, b2*b3=c^2, Σ(T*)^2=0, Prod(T*)=(N-2)T*)"""
    # u^2 + (N-2)u + 1 = 0 の負実根
    u = 0.5 * (-(N - 2.0) + math.sqrt((N - 2.0) ** 2 - 4.0))
    b2 = 1j * math.sqrt(-u)
    b3 = 1.0 / b2  # = -i / sqrt(-u)
    c = 1.0 + 0j
    T = np.full((N, N, N), c, dtype=np.complex128)
    idx = np.arange(N)
    # R2: i == k != j -> b2
    for i in range(N):
        T[i, :, i] = b2
    # R3: i == j != k -> b3
    for i in range(N):
        T[i, i, :] = b3
    # R1: j == k != i -> 0, R0: i == j == k -> 0
    for j in range(N):
        T[:, j, j] = 0.0
    return T


def build_antisym_tensor(N, a=1, b=2):
    w = np.exp(2j * np.pi / N)
    idx = np.arange(N)
    u = (idx[None, :, None] - idx[:, None, None]) % N
    v = (idx[None, None, :] - idx[:, None, None]) % N
    return w ** ((a * u + b * v) % N) - w ** ((a * v + b * u) % N)


def run_defect(T0, pos, eps, max_steps=100):
    U = T0.copy()
    U[pos] += eps
    U = U / np.linalg.norm(U)
    hist = [U]
    for step in range(1, max_steps + 1):
        P = bm_prod(U)
        nrm = float(np.linalg.norm(P))
        if nrm < 1e-12:
            return f"zero at step {step}", step, 0, 0.0
        U_next = P / nrm
        for prev_i in range(max(0, len(hist) - 14), len(hist)):
            if proj_dist(hist[prev_i], U_next) < 1e-10 and step >= 15:
                p = step - prev_i
                nz = int(np.sum(np.abs(U_next) > 1e-5))
                return f"converged period={p} at step {step}", step, nz, sum_z2_ratio(U_next)
        hist.append(U_next)
        U = U_next
    nz = int(np.sum(np.abs(U) > 1e-5))
    min_rec = min(proj_dist(h, U) for h in hist[-31:-1])
    return f"running at step {max_steps} (min_rec={min_rec:.2e})", max_steps, nz, sum_z2_ratio(U)


def run():
    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    t0 = time.time()
    out("=" * 72)
    out("大 N での純粋 BM 積 1 マス欠陥検証 (結晶解 N=12,16,24 / 反対称波 N=11,13)")
    out("=" * 72)

    out("\n--- 1. S_N 結晶解 T* (ΣT^2 = 0) の大 N 検証 ---")
    defects_sn = [
        ("対角 (0,0,0) +iε", (0, 0, 0), 1e-4j),
        ("対角 (0,0,0) +ε", (0, 0, 0), 1e-4),
        ("b1 (1,0,0) +ε", (1, 0, 0), 1e-4),
        ("b2 (0,1,0) +ε", (0, 1, 0), 1e-4),
        ("c (0,1,2) +iε", (0, 1, 2), 1e-4j),
    ]
    for N in (12, 16, 24):
        T0 = build_sn_crystal_sz2_zero(N)
        P0 = bm_prod(T0)
        err_p = float(np.max(np.abs(P0 - (N - 2) * T0)))
        err_z2 = sum_z2_ratio(T0)
        assert err_p < 1e-10 and err_z2 < 1e-12
        out(f"\n[N = {N}] (成分数 {N**3}, 検算誤差 {err_p:.1e}, ΣT^2={err_z2:.1e})")
        for label, pos, eps in defects_sn:
            res, st, nz, sz2 = run_defect(T0, pos, eps, max_steps=80)
            out(f"  {label:20s} → {res} (nz={nz}/{N**3}, Σz^2={sz2:.3e})")

    out("\n--- 2. 反対称波周期結晶 F_{1,2} (ΣT^2 = 0) の大 N 検証 ---")
    defects_as = [
        ("対角 (0,0,0) +ε", (0, 0, 0), 1e-4),
        ("対角 (0,0,0) +iε", (0, 0, 0), 1e-4j),
        ("非対角 (0,1,0) +ε", (0, 1, 0), 1e-4),
        ("R4 (0,1,2) +iε", (0, 1, 2), 1e-4j),
    ]
    for N in (11, 13):
        T0 = build_antisym_tensor(N, 1, 2)
        out(f"\n[N = {N}] (成分数 {N**3}, 理論周期 ord_N(2)={10 if N==11 else 12})")
        for label, pos, eps in defects_as:
            res, st, nz, sz2 = run_defect(T0, pos, eps, max_steps=80)
            out(f"  {label:20s} → {res} (nz={nz}/{N**3}, Σz^2={sz2:.3e})")

    out(f"\n所要時間: {time.time() - t0:.2f} 秒")
    os.makedirs("results", exist_ok=True)
    out_path = os.path.join("results", "bm_product_largeN_defect_v0.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"結果を保存しました: {out_path}")


if __name__ == "__main__":
    run()
