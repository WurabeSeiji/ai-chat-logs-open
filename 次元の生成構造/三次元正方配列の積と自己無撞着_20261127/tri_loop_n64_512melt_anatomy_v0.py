#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N = 64 (O8 ⊗ O8) の実数 1 マス欠陥から生じた 512 マス複合溶け残り状態 (step 200) の
厳密な代数構造・テンソル積分解 (SVD)・2×2 複素共役構造・内部結合定数 (視野角 ρ) の動力学、
および 2 次摂動に対する復元性（ヒッグス的硬さ vs 位相モード vs 崩壊）を検証する。

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


def node_energies(U):
    a2 = np.abs(U) ** 2
    return (np.sum(a2, axis=(1, 2)) + np.sum(a2, axis=(0, 2)) + np.sum(a2, axis=(0, 1))) / 3.0


def analyze_2x2_complex(M):
    """
    2×2 実数行列 M = [[a, b], [c, d]] を
    複素数表現 (a=d, b=-c) と 共役表現 (a=-d, b=c) に直交分解し、それぞれのノルム比を返す。
    """
    a, b, c, d = M[0, 0], M[0, 1], M[1, 0], M[1, 1]
    # C-like (回転・スカラー倍): [[x, -y], [y, x]]
    x_rot = 0.5 * (a + d)
    y_rot = 0.5 * (c - b)
    n_rot = np.sqrt(2.0 * (x_rot ** 2 + y_rot ** 2))
    # Anti-C-like (鏡映・共役): [[u, v], [v, -u]]
    u_ref = 0.5 * (a - d)
    v_ref = 0.5 * (b + c)
    n_ref = np.sqrt(2.0 * (u_ref ** 2 + v_ref ** 2))
    return n_rot, n_ref, (x_rot, y_rot), (u_ref, v_ref)


def run():
    t0 = time.time()
    os.makedirs("results", exist_ok=True)
    lines = []

    def log(s=""):
        print(s)
        lines.append(s)

    O8 = np.array(make_mult_table(8), dtype=np.float64)
    T0 = tensor_prod_3d(O8, O8)

    U = T0.copy()
    U[0, 1, 0] += 1e-4
    U /= np.linalg.norm(U)

    for st in range(1, 201):
        U, nrm = step_tri(U)

    U200 = U.copy()
    log("=" * 78)
    log("1. N = 64 (O8 ⊗ O8) 512 マス状態 (step 200) の厳密代数構造")
    log("=" * 78)
    mask = np.abs(U200) > 1e-8
    idx_i = sorted(np.where(np.any(mask, axis=(1, 2)))[0].tolist())
    idx_j = sorted(np.where(np.any(mask, axis=(0, 2)))[0].tolist())
    idx_k = sorted(np.where(np.any(mask, axis=(0, 1)))[0].tolist())
    log(f"非ゼロ数 (>1e-8): {int(np.sum(mask))}, ||P||={nrm:.8f}")
    log(f"生存軸 i ({len(idx_i)}): {idx_i}")
    log(f"  -> (A=i//8, a=i%8): {[(x // 8, x % 8) for x in idx_i]}")
    log(f"生存軸 j ({len(idx_j)}): {idx_j}")
    log(f"  -> (B=j//8, b=j%8): {[(x // 8, x % 8) for x in idx_j]}")
    log(f"生存軸 k ({len(idx_k)}): {idx_k}")
    log(f"  -> (C=k//8, c=k%8): {[(x // 8, x % 8) for x in idx_k]}")

    # 6階テンソル (A, a, B, b, C, c) -> 外側 (A,B,C) × 内側 (a,b,c) の SVD
    U6 = U200.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5).reshape(512, 512)
    s_vals = np.linalg.svd(U6, compute_uv=False)
    log(f"\n外側 (A,B,C) × 内側 (a,b,c) の特異値 (上位 10 個):")
    log("  " + ", ".join(f"{v:.8f}" for v in s_vals[:10]))

    # ランク 1 なら外側・内側を抽出して O8 と比較
    U_svd, S_svd, Vt_svd = np.linalg.svd(U6, full_matrices=False)
    outer = U_svd[:, 0].reshape(8, 8, 8)
    inner = Vt_svd[0, :].reshape(8, 8, 8)
    O8_unit = O8 / np.linalg.norm(O8)
    O8_1 = tri_prod_fast(O8)
    O8_1_unit = O8_1 / np.linalg.norm(O8_1)
    log(f"外側第1特異ベクトルと O8 の射影距離: d(O8)={proj_dist(outer, O8_unit):.3e}, d(Tri(O8))={proj_dist(outer, O8_1_unit):.3e}")
    log(f"外側非ゼロ数 (>1e-6): {int(np.sum(np.abs(outer) > 1e-6))}")
    log(f"内側非ゼロ数 (>1e-6): {int(np.sum(np.abs(inner) > 1e-6))}")

    # 内側 8×8×8 の非ゼロ成分を全表示
    log("\n内側第1特異ベクトル (非ゼロ成分 >1e-6):")
    for a in range(8):
        for b in range(8):
            for c in range(8):
                if abs(inner[a, b, c]) > 1e-6:
                    log(f"  inner[{a},{b},{c}] = {inner[a,b,c]:+.10f}")

    # 内側が 2×2×2 ブロックになっているので、b=0,1 ごとに 2×2 行列を表示
    sub_a = sorted(set(x % 8 for x in idx_i))
    sub_b = sorted(set(x % 8 for x in idx_j))
    sub_c = sorted(set(x % 8 for x in idx_k))
    log(f"\n内側生存添字: a={sub_a}, b={sub_b}, c={sub_c}")
    for b in sub_b:
        M = inner[np.ix_(sub_a, [b], sub_c)][:, 0, :]
        n_rot, n_ref, (xr, yr), (ur, vr) = analyze_2x2_complex(M)
        log(f"  b={b}: M = [[{M[0,0]:+.8f}, {M[0,1]:+.8f}], [{M[1,0]:+.8f}, {M[1,1]:+.8f}]]")
        log(f"       rot_norm={n_rot:.8f} (x={xr:+.8f}, y={yr:+.8f}), ref_norm={n_ref:.8f} (u={ur:+.8f}, v={vr:+.8f})")

    # 3段階エネルギー階層の正体を確認
    ne = node_energies(U200)
    log("\n--- 64 要素のエネルギー階層 (step 200) ---")
    # 軸別エネルギーも見る
    a2 = U200 ** 2
    ei = np.sum(a2, axis=(1, 2))
    ej = np.sum(a2, axis=(0, 2))
    ek = np.sum(a2, axis=(0, 1))
    log("  A  a |   E_avg      E_i        E_j        E_k")
    for A in range(2):  # 周期8なので A=0,1 を詳細表示、残りも一致確認
        for a in range(8):
            idx = 8 * A + a
            if ne[idx] > 1e-10:
                log(f"  {A}  {a} | {ne[idx]:.8f}  {ei[idx]:.8f}  {ej[idx]:.8f}  {ek[idx]:.8f}")
    # 全 A=0..7 で a ごとの値が一様か確認
    log("\n  各 a (0..7) について A=0..7 のエネルギー分散:")
    for a in range(8):
        vals = [ne[8 * A + a] for A in range(8)]
        if max(vals) > 1e-10:
            log(f"    a={a}: mean={np.mean(vals):.10f}, std={np.std(vals):.3e}, min={min(vals):.10f}, max={max(vals):.10f}")

    # 3 段階エネルギーの時間変化 (step 200..210)
    log("\n--- ステップ 200〜210 のエネルギー階層と内側 2×2 回転角の推移 ---")
    Ut = U200.copy()
    for st in range(200, 211):
        net = node_energies(Ut)
        # 内側を SVD で取得
        U6t = Ut.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5).reshape(512, 512)
        Ut_u, St_s, Vt_v = np.linalg.svd(U6t, full_matrices=False)
        inn = Vt_v[0, :].reshape(8, 8, 8)
        M0 = inn[np.ix_(sub_a, [sub_b[0]], sub_c)][:, 0, :]
        M1 = inn[np.ix_(sub_a, [sub_b[1]], sub_c)][:, 0, :]
        nr0, nf0, (x0, y0), (u0, v0) = analyze_2x2_complex(M0)
        nr1, nf1, (x1, y1), (u1, v1) = analyze_2x2_complex(M1)
        th0 = np.degrees(np.arctan2(y0, x0)) if nr0 > nf0 else np.degrees(np.arctan2(v0, u0))
        th1 = np.degrees(np.arctan2(y1, x1)) if nr1 > nf1 else np.degrees(np.arctan2(v1, u1))
        uniq_e = sorted(set(round(float(x), 6) for x in net if x > 1e-8), reverse=True)
        log(f"  step {st}: s1={St_s[0]:.8f}, s2={St_s[1]:.2e}, E_levels={uniq_e}, "
            f"b={sub_b[0]}(rot={nr0:.4f},ref={nf0:.4f},ang={th0:+7.2f}°), "
            f"b={sub_b[1]}(rot={nr1:.4f},ref={nf1:.4f},ang={th1:+7.2f}°)")
        Ut, _ = step_tri(Ut)

    # 2. 内部の三角形視野角 ρ (実効結合定数) の全数解析 (step 200)
    log("\n" + "=" * 78)
    log("2. 512 マス内部の三角形視野角 ρ (= 実効結合定数) の全数分布 (step 200)")
    log("=" * 78)
    # 512 非ゼロ成分すべての寄与三角形を集計
    nz_coords = np.argwhere(np.abs(U200) > 1e-8)
    rhos_all = []
    # 代表として A=0, B=0, C=0 周辺の生存成分と全体を確認
    for (i, j, k) in nz_coords[:16]:  # 対称性が高いので先頭16個を全数計測
        # p, q, r は生存軸のみ見ればよい (他は 0)
        sub_p = idx_j  # U[i,p,q] の第2添字は j 軸 (idx_j), U[r,p,k] の第2添字も j 軸!
        # 正確には全 idx_i, idx_j, idx_k の積集合・和集合を使う
        # U[i, p, q]: p in idx_j, q in idx_k
        # U[q, j, r]: q in idx_i, r in idx_k
        # U[r, p, k]: r in idx_i, p in idx_j
        P_set = idx_j
        Q_set = sorted(set(idx_k) & set(idx_i))
        R_set = sorted(set(idx_k) & set(idx_i))
        a1 = np.abs(U200[i][np.ix_(P_set, Q_set)])  # (lenP, lenQ)
        a2m = np.abs(U200[:, j, :][np.ix_(Q_set, R_set)])  # (lenQ, lenR)
        a3 = np.abs(U200[:, :, k][np.ix_(R_set, P_set)])  # (lenR, lenP)
        # 3次元配列 (lenP, lenQ, lenR)
        e1 = a1[:, :, None]
        e2 = a2m[None, :, :]
        e3 = a3.transpose(1, 0)[:, None, :]  # (lenP, 1, lenR)
        tri_min = np.minimum(np.minimum(e1, e2), e3)
        tri_max = np.maximum(np.maximum(e1, e2), e3)
        valid = tri_min > 1e-12
        rho = (tri_min[valid] / tri_max[valid]).ravel()
        rhos_all.append(rho)
    rhos_cat = np.concatenate(rhos_all)
    uniq_rho, cnt_rho = np.unique(np.round(rhos_cat, 6), return_counts=True)
    log(f"先頭 16 生存成分の有効三角形総数: {len(rhos_cat)} (1成分あたり {len(rhos_cat)//16} 個)")
    log(f"ユニークな視野角 ρ 値 (小数6桁丸め): {len(uniq_rho)} 種類")
    for rv, cv in zip(uniq_rho, cnt_rho):
        log(f"  ρ = {rv:.6f}  (件数: {cv}, 1成分あたり {cv/16:.1f} 個)")

    # 3. 2 次摂動に対する復元性テスト（ヒッグス的硬さ vs ゴールドストーン・モード vs 崩壊）
    log("\n" + "=" * 78)
    log("3. 512 マス状態 (step 200) への 2 次摂動テスト (100 ステップ追跡)")
    log("=" * 78)

    # 基準軌道 (step 200 -> 300) を計算
    ref_traj = [U200.copy()]
    Ur = U200.copy()
    for _ in range(100):
        Ur, _ = step_tri(Ur)
        ref_traj.append(Ur.copy())

    def test_perturb(name, U_init):
        Up = U_init / np.linalg.norm(U_init)
        d0 = proj_dist(Up, ref_traj[0])
        # SVD 第1特異値 (直積性の指標)
        def s2_ratio(X):
            X6 = X.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5).reshape(512, 512)
            sv = np.linalg.svd(X6, compute_uv=False)
            return float(sv[0]), float(sv[1])

        s1_0, s2_0 = s2_ratio(Up)
        # 512 マス外部のエネルギー比
        out_mask = ~mask
        out_e0 = float(np.sum(Up[out_mask] ** 2))
        log(f"\n[実験] {name}")
        log(f"  t=0  : dist_ref={d0:.3e}, out_E={out_e0:.3e}, s1={s1_0:.6f}, s2={s2_0:.3e}")
        checkpoints = {1, 2, 5, 10, 20, 50, 100}
        for t in range(1, 101):
            Up, nrm_p = step_tri(Up)
            if t in checkpoints:
                dt = proj_dist(Up, ref_traj[t])
                out_et = float(np.sum(Up[out_mask] ** 2))
                nz_t = int(np.sum(np.abs(Up) > 1e-5))
                s1_t, s2_t = s2_ratio(Up)
                # 外側 O8 との距離も見る
                X6 = Up.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5).reshape(512, 512)
                Uu, _, _ = np.linalg.svd(X6, full_matrices=False)
                out_t = Uu[:, 0].reshape(8, 8, 8)
                d_O8 = min(proj_dist(out_t, O8_unit), proj_dist(out_t, O8_1_unit))
                log(f"  t={t:3d}: nz={nz_t:4d}, ||P||={nrm_p:.6f}, dist_ref={dt:.3e}, "
                    f"out_E={out_et:.3e}, s2={s2_t:.3e}, d(outer,O8)={d_O8:.3e}")

    # (A) 背景（0領域）への 1 マス摂動 (+1e-3, +1e-1)
    Ua1 = U200.copy(); Ua1[0, 0, 0] += 1e-3
    test_perturb("(A1) 背景0領域への1マス摂動 (0,0,0) +1e-3", Ua1)

    Ua2 = U200.copy(); Ua2[0, 0, 0] += 1e-1
    test_perturb("(A2) 背景0領域への大摂動 (0,0,0) +0.1", Ua2)

    # (B) 内部512マスの1マスだけに微小摂動 (+1e-4, +1e-2) -> 外側 O8 の対称性を破る!
    p_coord = tuple(int(x) for x in nz_coords[0])
    Ub1 = U200.copy(); Ub1[p_coord] += 1e-4
    test_perturb(f"(B1) 内部512マスへの1マス局所摂動 {p_coord} +1e-4", Ub1)

    Ub2 = U200.copy(); Ub2[p_coord] += 1e-2
    test_perturb(f"(B2) 内部512マスへの1マス局所摂動 {p_coord} +1e-2", Ub2)

    # (C) 内部512マスの1マスを完全にゼロ化（穴あけ欠陥！）
    Uc = U200.copy(); Uc[p_coord] = 0.0
    test_perturb(f"(C) 内部512マスの1マス完全ゼロ化 {p_coord} = 0", Uc)

    # (D) 内側2×2複素構造の「非共役・非回転」成分への大域摂動（外側 O8 は保つ）
    inn_pert = inner.copy()
    inn_pert[sub_a[0], sub_b[0], sub_c[0]] += 1e-2  # 2×2 の1成分だけ崩す
    Ud = np.einsum("ijk,abc->iajbkc", outer, inn_pert).reshape(64, 64, 64)
    test_perturb("(D) 外側O8保持・内側2×2複素構造の1成分崩し (+1e-2)", Ud)

    # (E) 内側2×2の「位相だけ」を一様に少し回転させる（ゴールドストーン方向）
    ang = np.radians(1.0)
    R2 = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
    inn_rot = inner.copy()
    for b in sub_b:
        M = inn_rot[np.ix_(sub_a, [b], sub_c)][:, 0, :]
        M_new = R2 @ M
        for ia, a in enumerate(sub_a):
            for ic, c in enumerate(sub_c):
                inn_rot[a, b, c] = M_new[ia, ic]
    Ue = np.einsum("ijk,abc->iajbkc", outer, inn_rot).reshape(64, 64, 64)
    test_perturb("(E) 外側O8保持・内側2×2位相の1°回転摂動", Ue)

    log(f"\n所要時間: {time.time() - t0:.2f} 秒")
    out_path = os.path.join("results", "tri_loop_n64_512melt_anatomy_v0.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    log(f"保存先: {out_path}")


if __name__ == "__main__":
    run()
