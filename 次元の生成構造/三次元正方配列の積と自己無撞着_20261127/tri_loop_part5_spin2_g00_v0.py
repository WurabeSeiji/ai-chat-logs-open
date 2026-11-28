#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Part V 追加検証 [MS-6]：
多殻大格子 K=2 (N=512) において、中心 d=0 にスピン 2 (σ_x) モード局在質量
  P_{spin2} = δ_{000}^{⊗2} ⊗ M_{sigma_x}
を投入したとき、弱結合共存相で殻深さ d=0, 1, 2 に応じた局所固有時間遅延 g_00(d)
（時計速度 φ_bg(d)/φ_0 または殻間位相差 Δθ(d) の距離勾配）が解禁されるかを厳密閉包で検証する。
乱数なし。
"""
import os
import time
import numpy as np
from bm_product_algebra_tables_v0 import make_mult_table
from tri_loop_large_crystal_v0 import tri_prod_fast, tensor_prod_3d
from tri_loop_part5_multishell_gravity_v0 import (
    build_closure,
    build_structure_tensor,
    tri_out_triple_recursive,
)


def step_tri_K_gen(C, G_out, G_in, K):
    r_in = G_in.shape[0]
    nz_m = [m for m in range(r_in) if np.max(np.abs(C[..., m])) > 1e-15]
    P = np.zeros_like(C)
    cache = {}
    for m in nz_m:
        for n in nz_m:
            for p in nz_m:
                gvec = G_in[m, n, p, :]
                if np.max(np.abs(gvec)) < 1e-15:
                    continue
                key = (m, n, p)
                if key not in cache:
                    cache[key] = tri_out_triple_recursive(C[..., m], C[..., n], C[..., p], G_out, K)
                P += np.einsum("...,q->...q", cache[key], gvec)
    nrm = float(np.linalg.norm(P))
    return (P / nrm, nrm) if nrm > 1e-15 else (np.zeros_like(C), 0.0)


def extract_bg_vecs(C, b_out, b_in, K):
    """各殻深さ d=0..K の背景チャネル (b in {0,1}) ベクトル (正規化前) と複素時計位相 θ(d) を返す"""
    v0 = np.array([b[0, 0, 0] for b in b_out], dtype=np.float64)
    v1 = np.array([b[1, 0, 1] for b in b_out], dtype=np.float64)
    bg_tensors = np.array([b[:, 0:2, :] for b in b_in], dtype=np.float64)  # (r_in, 8, 2, 8)
    c0_vec = np.array([0.5 * (b[4, 0, 4] + b[5, 0, 5]) for b in b_in], dtype=np.float64)
    cJ_vec = np.array([0.5 * (b[5, 0, 4] - b[4, 0, 5]) for b in b_in], dtype=np.float64)
    vecs = []
    Es = np.zeros(K + 1)
    ths = np.zeros(K + 1)
    for d in range(K + 1):
        ten = C
        for s in range(K):
            vec = v0 if s < (K - d) else v1
            ten = np.tensordot(vec, ten, axes=([0], [0]))
        bg_loc = np.tensordot(ten, bg_tensors, axes=([0], [0])).ravel()
        vecs.append(bg_loc)
        Es[d] = float(np.dot(bg_loc, bg_loc))
        ths[d] = float(np.arctan2(float(np.dot(ten, cJ_vec)), float(np.dot(ten, c0_vec))))
    return vecs, Es, ths


def run():
    t0 = time.time()
    os.makedirs("results", exist_ok=True)
    lines = []

    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    O8 = np.array(make_mult_table(8), dtype=np.float64)
    D0 = np.zeros((8, 8, 8), dtype=np.float64)
    D0[0, 0, 0] = 1.0

    T0 = tensor_prod_3d(O8, O8)
    U = T0.copy()
    U[0, 1, 0] += 1e-4
    U /= np.linalg.norm(U)
    for _ in range(200):
        P = tri_prod_fast(U)
        U = P / np.linalg.norm(P)
    H_mat = U.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5).reshape(512, 512)
    _, _, vt = np.linalg.svd(H_mat, full_matrices=False)
    M8 = vt[0].reshape(8, 8, 8)

    # スピン 2 (σ_x) 内側プローブ (Part I §6 と完全同一)
    Msx = np.zeros((8, 8, 8), dtype=np.float64)
    for a in [4, 5]:
        for b in [2, 3]:
            for c in [4, 5]:
                if a != c:
                    Msx[a, b, c] = 1.0
    Msx /= np.linalg.norm(Msx)

    b_out = build_closure([O8, D0], max_dim=16)
    b_in = build_closure([M8, Msx], max_dim=20)
    G_out = build_structure_tensor(b_out)
    G_in = build_structure_tensor(b_in)

    c_O8 = np.array([float(np.sum(O8 * b)) for b in b_out]); c_O8 /= np.linalg.norm(c_O8)
    c_D0 = np.array([float(np.sum(D0 * b)) for b in b_out]); c_D0 /= np.linalg.norm(c_D0)
    c_M8 = np.array([float(np.sum(M8 * b)) for b in b_in]); c_M8 /= np.linalg.norm(c_M8)
    c_Msx = np.array([float(np.sum(Msx * b)) for b in b_in]); c_Msx /= np.linalg.norm(c_Msx)

    log("=" * 86)
    log("Part V 追加検証 [MS-6]：スピン 2 (σ_x) 局在質量による多殻 g_00(d) 解禁テスト (K=2, N=512)")
    log("=" * 86)
    log(f"  厳密閉包次元：外側 1 殻 = {len(b_out)}D, 内側 {{M8, M_sigma_x}} = {len(b_in)}D")

    # 基準時計 φ_0 (真空 K=2, τ=20..40 平均)
    K = 2
    C_vac = np.tensordot(np.tensordot(c_O8, c_O8, axes=0), c_M8, axes=0)
    C_vac /= np.linalg.norm(C_vac)
    phi0_l = []
    for st in range(1, 41):
        Cn, _ = step_tri_K_gen(C_vac, G_out, G_in, K)
        if st >= 20:
            vp, _, _ = extract_bg_vecs(C_vac, b_out, b_in, K)
            vc, _, _ = extract_bg_vecs(Cn, b_out, b_in, K)
            cval = float(np.clip(np.dot(vp[0], vc[0]) / (np.linalg.norm(vp[0]) * np.linalg.norm(vc[0])), -1.0, 1.0))
            phi0_l.append(float(np.arccos(abs(cval))))
        C_vac = Cn
    phi0 = float(np.mean(phi0_l))
    log(f"  [Ref] 真空 K=2 基準時計 φ_0 = {phi0:.6f} rad/step\n")

    log(f"  {'ε_sx':>7s} | {'||P||(40)':>9s} | {'E(0)/E(2)':>10s} | {'φ(0)/φ_0':>9s} {'φ(1)/φ_0':>9s} {'φ(2)/φ_0':>9s} | {'Δθ(0) [rad]':>12s} {'Δθ(1) [rad]':>12s}")
    for eps in [0.002, 0.005, 0.008, 0.012, 0.015, 0.025]:
        bg = np.tensordot(c_O8, c_O8, axes=0)
        pr = np.tensordot(c_D0, c_D0, axes=0)
        C = np.tensordot(bg, c_M8, axes=0) + eps * np.tensordot(pr, c_Msx, axes=0)
        C /= np.linalg.norm(C)
        phis_acc = [[] for _ in range(K + 1)]
        for st in range(1, 41):
            Cn, nP = step_tri_K_gen(C, G_out, G_in, K)
            if st >= 20:
                vp, _, _ = extract_bg_vecs(C, b_out, b_in, K)
                vc, _, _ = extract_bg_vecs(Cn, b_out, b_in, K)
                for d in range(K + 1):
                    np_v, nc_v = float(np.linalg.norm(vp[d])), float(np.linalg.norm(vc[d]))
                    if np_v > 1e-15 and nc_v > 1e-15:
                        cval = float(np.clip(np.dot(vp[d], vc[d]) / (np_v * nc_v), -1.0, 1.0))
                        phis_acc[d].append(float(np.arccos(abs(cval))))
            C = Cn
        _, Ed, thd = extract_bg_vecs(C, b_out, b_in, K)
        r_phi = [float(np.mean(phis_acc[d])) / phi0 if phis_acc[d] else 0.0 for d in range(K + 1)]
        dth0 = (thd[0] - thd[2] + np.pi) % (2 * np.pi) - np.pi
        dth1 = (thd[1] - thd[2] + np.pi) % (2 * np.pi) - np.pi
        rat02 = Ed[0] / Ed[2] if Ed[2] > 1e-300 else np.inf
        log(f"  {eps:7.3f} | {nP:9.6f} | {rat02:10.4f} | {r_phi[0]:9.6f} {r_phi[1]:9.6f} {r_phi[2]:9.6f} | {dth0:+12.3e} {dth1:+12.3e}")

    log(f"\nDone in {time.time() - t0:.2f}s.")
    with open("results/tri_loop_part5_spin2_g00_v0.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    run()
