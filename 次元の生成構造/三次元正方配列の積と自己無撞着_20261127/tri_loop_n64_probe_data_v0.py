#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
図 2 用データ生成: 512 マス時間結晶背景 H = O8 ⊗ M8 と 5 種類の近傍プローブ
  - P3a: 大域複素回転波 (完全共存)
  - P3b: 大域ヌル/カイラル波 (吸収 & ||P|| 2倍相転移)
  - P3c: 大域スピン2(σ_x)波 (エネルギー10倍吸収 & 周期2ロック)
  - P4 : 局在単一サイト凝縮体 (衣をまとって安定束縛 & 減速)
  - P5 : 局在四元数凝縮体 (本体:衣 = 1:1 安定束縛)
の τ = 0 -> 50 全ステップ時系列を results/tri_loop_n64_probe_timeseries.npz に保存する。
乱数なし。
"""
import os
import time
import numpy as np
from bm_product_algebra_tables_v0 import make_mult_table
from tri_loop_large_crystal_v0 import tri_prod_fast, tensor_prod_3d, proj_dist


def step_tri(U):
    P = tri_prod_fast(U)
    nrm = float(np.linalg.norm(P))
    if nrm < 1e-15:
        return np.zeros_like(U), 0.0
    return P / nrm, nrm


def run():
    t0 = time.time()
    os.makedirs("results", exist_ok=True)

    O8 = np.array(make_mult_table(8), dtype=np.float64)
    H4 = np.array(make_mult_table(4), dtype=np.float64)
    T0 = tensor_prod_3d(O8, O8)

    U = T0.copy()
    U[0, 1, 0] += 1e-4
    U /= np.linalg.norm(U)
    for _ in range(200):
        U, _ = step_tri(U)
    H = U.copy()

    H6 = H.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5).reshape(512, 512)
    Uu, _, Vtv = np.linalg.svd(H6, full_matrices=False)
    outer = Uu[:, 0].reshape(8, 8, 8)
    inner = Vtv[0, :].reshape(8, 8, 8)
    if np.sum(outer * O8) < 0:
        outer = -outer
        inner = -inner

    H_c = H.astype(np.complex128)
    mask_H = np.abs(H_c) > 1e-8

    D_orth_b = np.zeros_like(inner)
    for ib, b in enumerate([2, 3]):
        M = inner[np.ix_([4, 5], [ib], [4, 5])][:, 0, :]
        for ia, a in enumerate([4, 5]):
            for ic, c in enumerate([4, 5]):
                D_orth_b[a, b, c] = M[ia, ic]

    proj_null = np.array([[1.0, -1.0j], [1.0j, 1.0]], dtype=np.complex128)
    P3b_in = np.zeros((8, 8, 8), dtype=np.complex128)
    for b in [2, 3]:
        for ia, a in enumerate([4, 5]):
            for ic, c in enumerate([4, 5]):
                P3b_in[a, b, c] = proj_null[ia, ic]

    sx = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
    P3c_in = np.zeros((8, 8, 8), dtype=np.complex128)
    for b in [2, 3]:
        for ia, a in enumerate([4, 5]):
            for ic, c in enumerate([4, 5]):
                P3c_in[a, b, c] = sx[ia, ic]

    out_loc1 = np.zeros((8, 8, 8), dtype=np.complex128)
    out_loc1[0, 0, 0] = 1.0
    out_H4 = np.zeros((8, 8, 8), dtype=np.complex128)
    out_H4[:4, :4, :4] = H4

    probes = {
        "P3a": tensor_prod_3d(outer.astype(np.complex128), D_orth_b.astype(np.complex128)),
        "P3b": tensor_prod_3d(outer.astype(np.complex128), P3b_in),
        "P3c": tensor_prod_3d(outer.astype(np.complex128), P3c_in),
        "P4": tensor_prod_3d(out_loc1, D_orth_b.astype(np.complex128)),
        "P5": tensor_prod_3d(out_H4, D_orth_b.astype(np.complex128)),
    }

    eps = 0.1
    max_tau = 50
    save = {"tau": np.arange(0, max_tau + 1, dtype=np.int32)}

    for key, Pr in probes.items():
        Pr = Pr / np.linalg.norm(Pr)
        mask_P = np.abs(Pr) > 1e-8
        Um = H_c + eps * Pr
        Um /= np.linalg.norm(Um)

        eH_arr = np.zeros(max_tau + 1)
        eP_arr = np.zeros(max_tau + 1)
        eO_arr = np.zeros(max_tau + 1)
        nrm_arr = np.zeros(max_tau + 1)
        d1P_arr = np.zeros(max_tau + 1)

        P0 = tri_prod_fast(Um)
        eH_arr[0] = float(np.sum(np.abs(Um[mask_H]) ** 2))
        eP_arr[0] = float(np.sum(np.abs(Um[mask_P]) ** 2))
        eO_arr[0] = max(0.0, 1.0 - eH_arr[0] - eP_arr[0])
        nrm_arr[0] = float(np.linalg.norm(P0))

        for tau in range(1, max_tau + 1):
            Um_next, nrm = step_tri(Um)
            d1p = proj_dist(Um[mask_P], Um_next[mask_P]) if np.linalg.norm(Um_next[mask_P]) > 1e-12 else 0.0
            Um = Um_next
            eH_arr[tau] = float(np.sum(np.abs(Um[mask_H]) ** 2))
            eP_arr[tau] = float(np.sum(np.abs(Um[mask_P]) ** 2))
            eO_arr[tau] = max(0.0, 1.0 - eH_arr[tau] - eP_arr[tau])
            nrm_arr[tau] = nrm
            d1P_arr[tau] = d1p

        save[f"{key}_eH"] = eH_arr
        save[f"{key}_eP"] = eP_arr
        save[f"{key}_eO"] = eO_arr
        save[f"{key}_nrm"] = nrm_arr
        save[f"{key}_d1P"] = d1P_arr
        print(f"done {key}: t=50 eH={eH_arr[-1]:.5f}, eP={eP_arr[-1]:.5f}, eO={eO_arr[-1]:.5f}, ||P||={nrm_arr[-1]:.5f}", flush=True)

    out_npz = os.path.join("results", "tri_loop_n64_probe_timeseries.npz")
    np.savez_compressed(out_npz, **save)
    print(f"Saved {out_npz} in {time.time()-t0:.2f}s")


if __name__ == "__main__":
    run()
