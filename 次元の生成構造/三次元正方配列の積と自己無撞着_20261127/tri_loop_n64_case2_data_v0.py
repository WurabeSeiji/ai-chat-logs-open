#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
図 1 用データ生成: Case 2 (大域丸め誤差種 ~1e-16) の τ = 0 -> 7,200 時系列と
最終 8 軸エネルギー分布を results/tri_loop_n64_case2_layer4.npz に保存する。
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
    T0 = tensor_prod_3d(O8, O8)
    U = T0.copy()
    U[0, 1, 0] += 1e-4
    U /= np.linalg.norm(U)
    for _ in range(200):
        U, _ = step_tri(U)

    H6 = U.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5).reshape(512, 512)
    Uu, _, Vtv = np.linalg.svd(H6, full_matrices=False)
    outer = Uu[:, 0].reshape(8, 8, 8)
    inner = Vtv[0, :].reshape(8, 8, 8)

    U = tensor_prod_3d(outer, inner)
    U /= np.linalg.norm(U)

    max_tau = 7200
    taus, nzs, neffs, nrms, d1s, sv4s = [], [], [], [], [], []
    eA_snaps = {}

    def outer_axis_energy(X):
        a2 = (np.abs(X) ** 2).reshape(8, 8, 8, 8, 8, 8)
        ei = np.sum(a2, axis=(1, 2, 3, 4, 5))
        ej = np.sum(a2, axis=(0, 1, 3, 4, 5))
        ek = np.sum(a2, axis=(0, 1, 2, 3, 5))
        return (ei + ej + ek) / 3.0

    for tau in range(0, max_tau + 1):
        if tau == 0:
            P0 = tri_prod_fast(U)
            nrm = float(np.linalg.norm(P0))
            d1 = 0.0
        else:
            U_next, nrm = step_tri(U)
            d1 = proj_dist(U, U_next)
            U = U_next

        keep = (
            (tau <= 160 and tau % 2 == 0)
            or (6200 <= tau <= 6600 and tau % 5 == 0)
            or (tau % 50 == 0)
        )
        if keep:
            absU = np.abs(U)
            nz = int(np.sum(absU > 1e-5))
            ipr = float(np.sum(absU ** 4))
            n_eff = (1.0 / ipr) if ipr > 0 else 0.0
            U6 = U.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5).reshape(512, 512)
            sv = np.linalg.svd(U6, compute_uv=False)
            taus.append(tau)
            nzs.append(nz)
            neffs.append(n_eff)
            nrms.append(nrm)
            d1s.append(d1)
            sv4s.append(sv[:4].copy())

        if tau in (0, 1500, 7000):
            eA_snaps[f"eA_{tau}"] = outer_axis_energy(U)

    out_npz = os.path.join("results", "tri_loop_n64_case2_layer4.npz")
    save = {
        "tau": np.array(taus, dtype=np.int32),
        "nz": np.array(nzs, dtype=np.int32),
        "n_eff": np.array(neffs, dtype=np.float64),
        "norm_P": np.array(nrms, dtype=np.float64),
        "d1": np.array(d1s, dtype=np.float64),
        "sv4": np.array(sv4s, dtype=np.float64),
    }
    save.update(eA_snaps)
    np.savez_compressed(out_npz, **save)
    print(f"Saved {out_npz} in {time.time()-t0:.2f}s (eA_7000={np.round(eA_snaps['eA_7000'], 5).tolist()})", flush=True)


if __name__ == "__main__":
    run()
