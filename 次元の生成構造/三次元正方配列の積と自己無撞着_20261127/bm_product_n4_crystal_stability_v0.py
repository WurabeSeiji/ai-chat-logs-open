#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
B: Step 2 で見つかった N = 4 厳密結晶解 T* (a=0, b1=0, b2=+i, b3=-i, c=1) の
面外安定性解析 (全 64 方向のヤコビアン固有値) と、1 マス欠陥を入れたときの非線形発展。

注意:
  - 乱数は一切使わない。
  - BM 積 Prod(T, T, T) は 3 次斉次 (Prod(sT) = s^3 Prod(T)) である。
    非線形発展で成分が浮動小数の上限を超えるのを防ぐため、共通のスカラー係数 s_n を
    外へ括り出して U_n = T_n / s_n として追跡する (向きと比は厳密に不変)。
    最初の 4 ステップはスケーリングなしの直接計算と一致することも検算する。
"""
import os
import numpy as np


def build_n4_crystal():
    N = 4
    T = np.zeros((N, N, N), dtype=np.complex128)
    for i in range(N):
        for j in range(N):
            for k in range(N):
                if i == j == k:
                    T[i, j, k] = 0.0
                elif j == k:
                    T[i, j, k] = 0.0
                elif i == k:
                    T[i, j, k] = 1j
                elif i == j:
                    T[i, j, k] = -1j
                else:
                    T[i, j, k] = 1.0
    return T


def bm_prod(T):
    # Prod(T, T, T)_{ijk} = sum_t T_{itk} T_{ijt} T_{tjk}
    return np.einsum("itk,ijt,tjk->ijk", T, T, T)


def bm_jacobian(T):
    """T における BM 積の Fréchet 微分 (64x64 複素行列)"""
    N = T.shape[0]
    dim = N ** 3
    J = np.zeros((dim, dim), dtype=np.complex128)
    for idx in range(dim):
        E = np.zeros((N, N, N), dtype=np.complex128)
        E.flat[idx] = 1.0
        dP = (
            np.einsum("itk,ijt,tjk->ijk", E, T, T)
            + np.einsum("itk,ijt,tjk->ijk", T, E, T)
            + np.einsum("itk,ijt,tjk->ijk", T, T, E)
        )
        J[:, idx] = dP.reshape(-1)
    return J


def normalized_sum_z2(T):
    """|Σ T_{ijk}^2| / Σ |T_{ijk}|^2 (スケール不変な Σz^2=0 の検算指標)"""
    denom = np.sum(np.abs(T) ** 2)
    if denom == 0:
        return 0.0
    return abs(np.sum(T ** 2)) / denom


def run():
    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=" * 72)
    out("B: N = 4 厳密結晶解 T* (a=0, b1=0, b2=+i, b3=-i, c=1) の面外安定性と欠陥発展")
    out("=" * 72)

    T0 = build_n4_crystal()
    P0 = bm_prod(T0)
    assert np.allclose(P0, 2.0 * T0, atol=1e-15), "結晶解 T* の検算失敗"
    assert abs(np.sum(T0 ** 2)) < 1e-15, "ΣT^2 = 0 の検算失敗"
    out("[検算合格] Prod(T*) = 2 T*,  Σ(T*)^2 = 0.\n")

    # 1. 全 64 方向のヤコビアン固有値 (倍率 λ = 2 で割った射影空間の増幅率 μ = eig(J)/2)
    J = bm_jacobian(T0)
    evals = np.linalg.eigvals(J) / 2.0
    # 固有値の重複度を集計 (小数 6 桁で丸めてグループ化)
    groups = {}
    for ev in evals:
        r = 0.0 if abs(ev.real) < 1e-7 else round(ev.real, 5)
        im = 0.0 if abs(ev.imag) < 1e-7 else round(ev.imag, 5)
        key = (r, im)
        groups[key] = groups.get(key, 0) + 1

    out("--- 1. 全 64 方向の線形安定性 (J / 2 の固有値一覧) ---")
    out("  (注意: μ = 3 の 1 本は T* 自身のスカラー倍方向。残り 63 本が向き・形の摂動)")
    for (r, im), cnt in sorted(groups.items(), key=lambda x: (-abs(complex(*x[0])), -x[0][0])):
        mag = abs(complex(r, im))
        ang = np.degrees(np.angle(complex(r, im))) if mag > 1e-7 else 0.0
        out(f"  μ = {r:+.5f} {im:+.5f}i  (|μ| = {mag:.5f}, 角 = {ang:+7.2f}°)  × {cnt} 重")

    # 2. 1 マス欠陥を入れた非線形発展
    out("\n--- 2. 1 マス欠陥を入れた非線形発展 (ε = 1e-4、30 ステップ) ---")
    defects = [
        ("対角 (0,0,0) に +ε", (0, 0, 0), 1e-4),
        ("対角 (0,0,0) に +iε", (0, 0, 0), 1e-4j),
        ("b1 位置 (1,0,0) に +ε", (1, 0, 0), 1e-4),
        ("b2 位置 (0,1,0) に +ε", (0, 1, 0), 1e-4),
        ("c 位置 (0,1,2) に +ε", (0, 1, 2), 1e-4),
        ("c 位置 (0,1,2) に +iε", (0, 1, 2), 1e-4j),
    ]

    # まず最初の 3 ステップが「スケーリングなし直接計算」と一致することを 1 例で検算
    Td_test = T0.copy()
    Td_test[0, 1, 2] += 1e-4
    U_direct = Td_test.copy()
    U_scaled = Td_test / np.max(np.abs(Td_test))
    for _ in range(3):
        U_direct = bm_prod(U_direct)
        U_scaled = bm_prod(U_scaled)
        U_scaled = U_scaled / np.max(np.abs(U_scaled))
    U_direct_norm = U_direct / np.max(np.abs(U_direct))
    assert np.allclose(U_direct_norm, U_scaled, atol=1e-12), "斉次スケーリング検算失敗"
    out("[検算合格] 3 次斉次スケーリング分離は直接計算と 1e-12 以内で一致。\n")

    T0_unit = T0 / np.linalg.norm(T0)
    for label, pos, eps in defects:
        U = T0.copy()
        U[pos] += eps
        U = U / np.linalg.norm(U)
        out(f"[{label}]")
        prev_U = None
        for step in range(1, 31):
            P = bm_prod(U)
            nrm = np.linalg.norm(P)
            if nrm < 1e-30:
                out(f"  step {step:2d}: 0 に消滅")
                break
            U_next = P / nrm
            # 結晶 T* との射影距離 (位相因子を除いた最小距離)
            overlap = np.vdot(T0_unit, U_next)
            dist_crystal = np.sqrt(max(0.0, 1.0 - abs(overlap) ** 2))
            # 直前ステップとの射影距離
            dist_prev = np.sqrt(max(0.0, 1.0 - abs(np.vdot(U, U_next)) ** 2))
            sz2 = normalized_sum_z2(U_next)
            if step in (1, 2, 3, 5, 10, 15, 20, 30):
                out(f"  step {step:2d}: dist(T*)={dist_crystal:.4e}, dist(prev)={dist_prev:.4e}, |Σz^2|/Σ|z|^2={sz2:.4e}, norm_ratio={nrm:.4f}")
            U = U_next
        # 最終状態の非ゼロ個数・構造要約
        nz = int(np.sum(np.abs(U) > 1e-6))
        out(f"  最終 (step 30): 非ゼロ={nz}/64, max|U|={np.max(np.abs(U)):.4f}")

    os.makedirs("results", exist_ok=True)
    out_path = os.path.join("results", "bm_product_n4_crystal_stability_v0.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\n結果を保存しました: {out_path}")


if __name__ == "__main__":
    run()
