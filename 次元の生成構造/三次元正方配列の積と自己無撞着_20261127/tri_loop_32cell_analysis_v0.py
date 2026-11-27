#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
§18 で見つかった八元数積表 (8^3) の虚数欠陥 +iε による 2 種類の溶け残り:
  (A) 336 件 (3/4): 32 マス周期 2 溶け残り (||P|| ≈ 0.2500, 例: (0,1,2)+iε)
  (B) 112 件 (1/4): 8 マス固定点溶け残り (||P|| ≈ 0.8140, 例: (0,0,1)+iε)
の内部構造・スライス行列・周期変化・ΣT^2・安定性、および 336 : 112 (3:1) に
分かれる位置条件を厳密に確認する。
"""
import os
import numpy as np
from bm_product_algebra_tables_v0 import make_mult_table


def tri_prod_np(T):
    return np.einsum("ipq,qjr,rpk->ijk", T, T, T, optimize=True)


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


def evolve_to_limit(T0, pos, eps, steps=60):
    U = T0.astype(np.complex128).copy()
    U[pos] += eps
    U /= np.linalg.norm(U)
    for _ in range(steps):
        P = tri_prod_np(U)
        U = P / np.linalg.norm(P)
    return U


def inspect_state(name, U, out):
    N = U.shape[0]
    P = tri_prod_np(U)
    nrm = float(np.linalg.norm(P))
    U1 = P / nrm
    P2 = tri_prod_np(U1)
    nrm2 = float(np.linalg.norm(P2))
    U2 = P2 / nrm2

    mask = np.abs(U) > 1e-5
    nz = int(np.sum(mask))
    i_set = sorted({i for i in range(N) if np.any(mask[i, :, :])})
    j_set = sorted({j for j in range(N) if np.any(mask[:, j, :])})
    k_set = sorted({k for k in range(N) if np.any(mask[:, :, k])})

    out(f"\n=== {name} ===")
    out(f"  非ゼロ数: {nz}/512,  ||P(U)|| = {nrm:.8f},  ||P(U1)|| = {nrm2:.8f}")
    out(f"  dist(U, U1) = {proj_dist(U, U1):.3e},  dist(U, U2) = {proj_dist(U, U2):.3e}")
    out(f"  |ΣU^2|/Σ|U|^2 = {sum_z2_ratio(U):.6f}")
    out(f"  非ゼロ添字集合: I = {i_set} ({len(i_set)}個), J = {j_set} ({len(j_set)}個), K = {k_set} ({len(k_set)}個)")

    # 各 j in j_set のスライス M_j = U[I, j, K] を表示 (位相基準を最大成分に揃える)
    nz_vals = U[mask]
    ref = nz_vals[np.argmax(np.abs(nz_vals))]
    ph = ref / abs(ref)
    U_al = U / ph

    for j in j_set:
        M = U_al[np.ix_(i_set, [j], k_set)][:, 0, :]
        # M^H M (エルミート直交性チェック)
        MhM = M.conj().T @ M
        diag_MhM = np.diag(MhM)
        off_MhM = MhM - np.diag(diag_MhM)
        is_unitary_scaled = float(np.max(np.abs(off_MhM))) < 1e-6 and float(np.max(np.abs(diag_MhM - diag_MhM[0]))) < 1e-6
        # M^T M (実直交/複素直交チェック)
        MTM = M.T @ M
        out(f"  [スライス j = {j}] 大きさ {M.shape}, ユニタリ定数倍か={is_unitary_scaled} (M^H M 対角={diag_MhM[0].real:.6f}, 非対角最大={np.max(np.abs(off_MhM)):.2e})")
        for r_idx, i in enumerate(i_set):
            row_s = "  ".join(f"{M[r_idx, c_idx].real:+.4f}{M[r_idx, c_idx].imag:+.4f}i" for c_idx in range(len(k_set)))
            out(f"    i={i}: {row_s}")

    # 周期 2 の場合、U と U1 の関係 (共役か、転置か、符号反転か)
    if proj_dist(U, U1) > 1e-4:
        # 同じ位相基準で比較
        d_conj = proj_dist(np.conj(U), U1)
        d_T = proj_dist(np.transpose(U, (2, 1, 0)), U1)
        out(f"  U と U1 の関係: dist(conj(U), U1)={d_conj:.3e}, dist(U^T_ik, U1)={d_T:.3e}")

    # 内部・外部の安定性 (微小摂動 1e-6 で 10 ステップ)
    # 内部: 最初の非ゼロ位置に +1e-6 と +1e-6j
    first_nz = tuple(int(x[0]) for x in np.where(mask))
    # 外部: (I, J, K) に含まれない最初の位置 (あれば)
    ext_pos = None
    for i in range(N):
        for j in range(N):
            for k in range(N):
                if i not in i_set and j not in j_set and k not in k_set:
                    ext_pos = (i, j, k)
                    break
            if ext_pos:
                break
        if ext_pos:
            break
    for plabel, ppos, peps in [
        (f"内部 {first_nz}+1e-6", first_nz, 1e-6),
        (f"内部 {first_nz}+1e-6j", first_nz, 1e-6j),
    ] + ([(f"完全外部 {ext_pos}+1e-6", ext_pos, 1e-6)] if ext_pos else []):
        V = U.copy()
        V[ppos] += peps
        V /= np.linalg.norm(V)
        Ur = U.copy()
        ds = []
        for _ in range(10):
            Ur = tri_prod_np(Ur); Ur /= np.linalg.norm(Ur)
            V = tri_prod_np(V); V /= np.linalg.norm(V)
            ds.append(proj_dist(Ur, V))
        out(f"  摂動 {plabel:24s}: step 1,2,5,10 後距離 = {[f'{ds[m]:.2e}' for m in (0, 1, 4, 9)]}")


def analyze_3to1_split(T0_int, out):
    """448 個の元ゼロ位置が 336 (nz=32) と 112 (nz=8) に分かれる幾何学的条件を特定する"""
    out("\n=== 448 個の元ゼロ位置が 336 (nz=32) : 112 (nz=8) に分かれる条件 ===")
    T0 = np.array(T0_int, dtype=np.complex128)
    cat_counts = {}
    for i in range(8):
        for j in range(8):
            for k in range(8):
                if T0_int[i][j][k] != 0:
                    continue
                U = evolve_to_limit(T0, (i, j, k), 1e-4j, steps=45)
                nz = int(np.sum(np.abs(U) > 1e-5))
                # 特徴量: (i==0, j==0, k==0, i==j, j==k, i==k)
                feat = (
                    "i=0" if i == 0 else "i>0",
                    "j=0" if j == 0 else "j>0",
                    "k=0" if k == 0 else "k>0",
                    "i=k" if i == k else ("i=j" if i == j else ("j=k" if j == k else "all_diff")),
                )
                key = (feat, nz)
                cat_counts[key] = cat_counts.get(key, 0) + 1
    for (feat, nz), cnt in sorted(cat_counts.items()):
        out(f"  {feat} -> nz={nz:2d}: {cnt:3d} 件")


def run():
    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=" * 72)
    out("八元数積表 (8^3) 虚数欠陥 +iε による 32 マス・8 マス溶け残りの内部構造解析")
    out("=" * 72)

    T0_int = make_mult_table(8)
    T0 = np.array(T0_int, dtype=np.complex128)

    U32 = evolve_to_limit(T0, (0, 1, 2), 1e-4j, steps=60)
    inspect_state("A: 32 マス周期 2 溶け残り (例: (0,1,2)+iε, 336/448 件)", U32, out)

    U8 = evolve_to_limit(T0, (0, 0, 1), 1e-4j, steps=60)
    inspect_state("B: 8 マス固定点溶け残り (例: (0,0,1)+iε, 112/448 件)", U8, out)

    analyze_3to1_split(T0_int, out)

    os.makedirs("results", exist_ok=True)
    out_path = os.path.join("results", "tri_loop_32cell_analysis_v0.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\n結果を保存しました: {out_path}")


if __name__ == "__main__":
    run()
