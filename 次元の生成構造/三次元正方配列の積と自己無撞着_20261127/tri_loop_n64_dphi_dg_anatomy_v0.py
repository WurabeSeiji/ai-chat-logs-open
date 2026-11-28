#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
δphi / δg ≈ 0.71 の厳密代数解剖スクリプト
  - P4 (A=0, eps=0.08, 0.16) において、サイト A=0 の 8 本の非ゼロエッジ (B, C) ごとに
    エネルギー変化 ΔE(B,C)、複素位相角 θ(0,B,C) vs θ(1,B',C')、および異方性漏れを全数分解。
  - さらに τ=20..80 の各ステップごとに、dphi(τ)/dg(τ) と 基準回転角 phi_0(τ) の関係を追跡し、
    δphi / δg ≈ 0.71 が定数（1/sqrt(2) や 5/7 等）なのか、phi_0(τ) の関数なのかを厳密判別する。
"""
import os
import numpy as np
from bm_product_algebra_tables_v0 import make_mult_table
from tri_loop_large_crystal_v0 import tri_prod_fast, tensor_prod_3d


def step_tri(U):
    P = tri_prod_fast(U)
    nrm = float(np.linalg.norm(P))
    return P / nrm, nrm


def wrap_pi(x):
    return (x + np.pi) % (2 * np.pi) - np.pi


def block_decomp(B2):
    c0 = 0.5 * (B2[0, 0] + B2[1, 1])
    cz = 0.5 * (B2[0, 0] - B2[1, 1])
    cx = 0.5 * (B2[0, 1] + B2[1, 0])
    cJ = 0.5 * (B2[1, 0] - B2[0, 1])
    r = float(np.hypot(c0, cJ))
    th = float(np.arctan2(cJ, c0))
    aniso = float(np.hypot(cx, cz))
    return c0, cJ, r, th, aniso


def run():
    O8 = np.array(make_mult_table(8), dtype=np.float64)
    T0 = tensor_prod_3d(O8, O8)
    U = T0.copy()
    U[0, 1, 0] += 1e-4
    U /= np.linalg.norm(U)
    for _ in range(200):
        U, _ = step_tri(U)
    H = U.copy()
    H_mat = H.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5).reshape(512, 512)
    _, _, vt_svd = np.linalg.svd(H_mat, full_matrices=False)
    M8 = vt_svd[0].reshape(8, 8, 8)

    P4 = np.zeros((64, 64, 64), dtype=np.float64)
    for a in [4, 5]:
        for b in [2, 3]:
            for c in [4, 5]:
                P4[0 * 8 + a, 0 * 8 + b, 0 * 8 + c] = M8[a, b - 2, c]
    P4 /= np.linalg.norm(P4)

    for eps in [0.08, 0.16]:
        print("=" * 78)
        print(f"厳密解剖: P4 at A=0, eps = {eps}")
        print("=" * 78)
        Ut = (H + eps * P4)
        Ut /= np.linalg.norm(Ut)
        for st in range(1, 51):
            Un, _ = step_tri(Ut)
            if st in [30, 31, 32, 40, 50]:
                T_p = Ut.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5)[:, :, :, :, 0:2, :]
                T_c = Un.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5)[:, :, :, :, 0:2, :]
                # サイト A=0 vs A=1 の全体 E, phi
                v0p, v0c = T_p[0].ravel(), T_c[0].ravel()
                v1p, v1c = T_p[1].ravel(), T_c[1].ravel()
                E0, E1 = float(np.dot(v0c, v0c)), float(np.dot(v1c, v1c))
                c0 = float(np.dot(v0p, v0c) / (np.linalg.norm(v0p) * np.linalg.norm(v0c)))
                c1 = float(np.dot(v1p, v1c) / (np.linalg.norm(v1p) * np.linalg.norm(v1c)))
                phi0 = float(np.arccos(abs(np.clip(c0, -1, 1))))
                phi1 = float(np.arccos(abs(np.clip(c1, -1, 1))))
                dg = 1.0 - E0 / E1
                dphi = 1.0 - phi0 / phi1
                # cos の差と tan(phi1) による厳密微分関係も確認
                # d(arccos(c)) = -dc / sin(phi)
                dc_abs = abs(c0) - abs(c1)
                print(
                    f"  tau={st:2d}: phi1={phi1:.6f} rad | dg={dg:.7f}, dphi={dphi:.7f}, "
                    f"dphi/dg={dphi/dg:.5f} | abs(c0)-abs(c1)={dc_abs:+.6e}"
                )
            if st == 40:
                # サイト A=0 の 8 本の非ゼロエッジ (B, C) を全数解剖 (b=0)
                print("  --- Site A=0 の非ゼロ (B,C) ブロック (b=0, tau=40) vs Site A=1 ---")
                T_c = Un.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5)
                # 基準として A=1, B=0, C=1 の (r_ref, th_ref)
                _, _, r_ref, th_ref, _ = block_decomp(T_c[1, 0, 1, 4:6, 0, 4:6])
                for B in range(8):
                    for C in range(8):
                        if abs(O8[0, B, C]) > 0.5:
                            _, _, r, th, an = block_decomp(T_c[0, B, C, 4:6, 0, 4:6])
                            # O8[0, B, C] の符号を補正して位相比較
                            sgn = O8[0, B, C]
                            th_adj = wrap_pi(th - (0.0 if sgn > 0 else np.pi) - th_ref)
                            print(
                                f"    (A=0, B={B}, C={C}) [O8={sgn:+.0f}]: r/r_ref={r/r_ref:.7f} "
                                f"(E/E_ref={(r/r_ref)**2:.7f}), dtheta={th_adj:+.6e} rad, aniso={an:.2e}"
                            )
                # さらに A=0 において O8[0,B,C]==0 なのに漏れている非ゼロ成分があるかチェック！
                leak_E = 0.0
                for B in range(8):
                    for C in range(8):
                        if abs(O8[0, B, C]) < 0.5:
                            leak_E += float(np.sum(T_c[0, B, C, :, 0:2, :] ** 2))
                print(f"    (A=0) O8=0 領域へのエネルギー漏れ (b=0,1): {leak_E:.6e}")
                # 同様に A=1 の 8 本の非ゼロエッジもチェック
                print("  --- Site A=1 の非ゼロ (B,C) ブロック (b=0, tau=40) ---")
                for B in range(8):
                    for C in range(8):
                        if abs(O8[1, B, C]) > 0.5:
                            _, _, r, th, an = block_decomp(T_c[1, B, C, 4:6, 0, 4:6])
                            sgn = O8[1, B, C]
                            th_adj = wrap_pi(th - (0.0 if sgn > 0 else np.pi) - th_ref)
                            print(
                                f"    (A=1, B={B}, C={C}) [O8={sgn:+.0f}]: r/r_ref={r/r_ref:.7f} "
                                f"(E/E_ref={(r/r_ref)**2:.7f}), dtheta={th_adj:+.6e} rad, aniso={an:.2e}"
                            )
            Ut = Un


if __name__ == "__main__":
    run()
