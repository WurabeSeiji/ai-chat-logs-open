#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
1: N = 4 厳密結晶解 T* (a=0, b1=0, b2=+i, b3=-i, c=1) に 1 マス欠陥を入れた
200 ステップ長期追跡と、最終構造（「溶け残り」や収束先・再帰）の確認。

注意:
  - 乱数なし。
  - 3 次斉次スケーリング分離 U_n = T_n / ||T_n|| を用いるが、
    浮動小数丸めノイズの増幅を防ぐため ||P|| < 1e-12 を「厳密ゼロ候補」として明示停止する。
"""
import os
import numpy as np


def build_n4_crystal():
    N = 4
    T = np.zeros((N, N, N), dtype=np.complex128)
    for i in range(N):
        for j in range(N):
            for k in range(N):
                if i == j == k or j == k:
                    T[i, j, k] = 0.0
                elif i == k:
                    T[i, j, k] = 1j
                elif i == j:
                    T[i, j, k] = -1j
                else:
                    T[i, j, k] = 1.0
    return T


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


def summarize_tensor(U, tol=1e-4):
    """非ゼロ成分の位置と値を要約"""
    N = U.shape[0]
    nz = []
    for i in range(N):
        for j in range(N):
            for k in range(N):
                if abs(U[i, j, k]) > tol:
                    nz.append((i, j, k, U[i, j, k]))
    return nz


def run():
    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=" * 72)
    out("1: N = 4 厳密結晶解 T* の 1 マス欠陥 200 ステップ長期追跡と最終構造の確認")
    out("=" * 72)

    T0 = build_n4_crystal()
    T0_unit = T0 / np.linalg.norm(T0)

    defects = [
        ("対角 (0,0,0) に +ε", (0, 0, 0), 1e-4),
        ("対角 (0,0,0) に +iε", (0, 0, 0), 1e-4j),
        ("b1 位置 (1,0,0) に +ε", (1, 0, 0), 1e-4),
        ("b2 位置 (0,1,0) に +ε", (0, 1, 0), 1e-4),
        ("c 位置 (0,1,2) に +ε", (0, 1, 2), 1e-4),
        ("c 位置 (0,1,2) に +iε", (0, 1, 2), 1e-4j),
    ]

    for label, pos, eps in defects:
        out(f"\n--- [{label}] ---")
        U = T0.copy()
        U[pos] += eps
        U = U / np.linalg.norm(U)
        hist = [U]
        status = "running"
        for step in range(1, 201):
            P = bm_prod(U)
            nrm = float(np.linalg.norm(P))
            if nrm < 1e-12:
                out(f"  step {step:3d}: ||P|| = {nrm:.3e} (< 1e-12) → 0 に消滅として停止")
                status = f"zero at step {step}"
                break
            U_next = P / nrm
            d_crys = proj_dist(T0_unit, U_next)
            d_prev = proj_dist(U, U_next)
            sz2 = sum_z2_ratio(U_next)
            if step in (10, 20, 30, 40, 50, 75, 100, 150, 200):
                nz_cnt = int(np.sum(np.abs(U_next) > 1e-6))
                out(f"  step {step:3d}: dist(T*)={d_crys:.4e}, dist(prev)={d_prev:.4e}, |Σz^2|/Σ|z|^2={sz2:.4e}, ||P||={nrm:.4f}, nz={nz_cnt}/64")
            # 固定点・周期収束チェック (step >= 5)
            found_p = None
            for prev_i in range(max(0, len(hist) - 12), len(hist)):
                if proj_dist(hist[prev_i], U_next) < 1e-10:
                    found_p = step - prev_i
                    break
            hist.append(U_next)
            U = U_next
            if found_p is not None and step >= 25:
                nz_cnt = int(np.sum(np.abs(U) > 1e-6))
                out(f"  step {step:3d}: 周期 {found_p} に収束 (誤差 < 1e-10), dist(T*)={d_crys:.4e}, |Σz^2|/Σ|z|^2={sz2:.4e}, nz={nz_cnt}/64")
                status = f"converged period={found_p} at step {step}"
                break

        if status == "running":
            # 後半 50 ステップでの最小再帰距離
            min_rec = min(proj_dist(h, U) for h in hist[-51:-1])
            status = f"non-converged at step 200 (min_rec_50={min_rec:.3e})"
        out(f"  判定: {status}")

        # 最終状態の非ゼロ成分を表示 (16 個以下なら全部)
        nz_list = summarize_tensor(U, tol=1e-4)
        out(f"  最終状態の非ゼロ成分数 (|u|>1e-4): {len(nz_list)}/64")
        if 0 < len(nz_list) <= 16:
            # 位相基準を最大絶対値成分の実数正に揃えて表示
            ref = max(nz_list, key=lambda x: abs(x[3]))[3]
            ph = ref / abs(ref)
            for (i, j, k, val) in nz_list:
                v = val / ph
                out(f"    ({i},{j},{k}): {v.real:+.6f} {v.imag:+.6f}i  (|v|={abs(v):.6f})")

    os.makedirs("results", exist_ok=True)
    out_path = os.path.join("results", "bm_product_n4_defect_long_v0.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\n結果を保存しました: {out_path}")


if __name__ == "__main__":
    run()
