#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N = 64 (O8 ⊗ O8, 262,144 成分, 初期非ゼロ 4,096) の実数 1 マス欠陥 (0,1,0)+1e-4 について、
初期の巨大結晶から融解・再結合を経て準安定状態へ至る全過程 (step 0〜500) を
厳密計算し、図化・スペクトル分析用の全データ (.csv, .npz, .txt) を保存する。

乱数なし。BLAS 行列積 2 回の厳密同値実装を使用。
"""
import csv
import os
import time
import numpy as np
from bm_product_algebra_tables_v0 import make_mult_table
from tri_loop_large_crystal_v0 import tri_prod_fast, tensor_prod_3d, proj_dist


def connected_components_from_tensor(U, tol=1e-5):
    """
    要素 0..N-1 を頂点とし、|U[i,j,k]| > tol が存在する組 (i,j), (j,k), (i,k) を
    無向辺で繋いだグラフの連結成分 (非孤立頂点のみ) のサイズ一覧を返す。
    """
    N = U.shape[0]
    mask = np.abs(U) > tol
    # 辺の隣接行列 (対称)
    adj = np.zeros((N, N), dtype=bool)
    e_ij = np.any(mask, axis=2)
    e_jk = np.any(mask, axis=0)
    e_ik = np.any(mask, axis=1)
    adj |= e_ij | e_ij.T | e_jk | e_jk.T | e_ik | e_ik.T
    active = np.any(adj, axis=1) | np.any(mask, axis=(1, 2)) | np.any(mask, axis=(0, 2)) | np.any(mask, axis=(0, 1))

    visited = [False] * N
    sizes = []
    members = []
    for start in range(N):
        if not active[start] or visited[start]:
            continue
        stack = [start]
        visited[start] = True
        comp = []
        while stack:
            v = stack.pop()
            comp.append(v)
            for nb in np.where(adj[v])[0]:
                nb = int(nb)
                if not visited[nb]:
                    visited[nb] = True
                    stack.append(nb)
        sizes.append(len(comp))
        members.append(sorted(comp))
    order = np.argsort(sizes)[::-1]
    sizes = [sizes[m] for m in order]
    members = [members[m] for m in order]
    return sizes, members


def run():
    t0 = time.time()
    os.makedirs("results", exist_ok=True)

    O8 = np.array(make_mult_table(8), dtype=np.float64)
    T0 = tensor_prod_3d(O8, O8)
    T1 = tri_prod_fast(T0)
    T0_unit = T0 / np.linalg.norm(T0)
    T1_unit = T1 / np.linalg.norm(T1)

    N = 64
    max_steps = 500
    snap_steps = [0, 20, 30, 40, 60, 80, 120, 200, 500]
    snap_set = set(snap_steps)

    U = T0.copy()
    U[0, 1, 0] += 1e-4
    U /= np.linalg.norm(U)

    rows = []
    snaps_map_ik = {}
    snaps_node_energy = {}
    snaps_abs_vals = {}
    snaps_comp_members = {}

    U_prev1 = U.copy()
    U_prev2 = U.copy()

    def record(step, U_cur, nrm, d1, d2):
        absU = np.abs(U_cur)
        nz_1e5 = int(np.sum(absU > 1e-5))
        nz_1e3 = int(np.sum(absU > 1e-3))
        nz_1e1 = int(np.sum(absU > 1e-1))
        mask = absU > 1e-5
        di = int(np.sum(np.any(mask, axis=(1, 2))))
        dj = int(np.sum(np.any(mask, axis=(0, 2))))
        dk = int(np.sum(np.any(mask, axis=(0, 1))))
        ipr = float(np.sum(absU ** 4))
        n_eff = (1.0 / ipr) if ipr > 0 else 0.0
        dist_crys = min(proj_dist(T0_unit, U_cur), proj_dist(T1_unit, U_cur))
        comp_sizes, comp_members = connected_components_from_tensor(U_cur, tol=1e-5)
        n_comp = len(comp_sizes)
        max_comp = comp_sizes[0] if comp_sizes else 0
        comp_str = "+".join(str(x) for x in comp_sizes) if comp_sizes else "0"
        rows.append({
            "step": step,
            "nz_1e5": nz_1e5,
            "nz_1e3": nz_1e3,
            "nz_1e1": nz_1e1,
            "dim_i": di,
            "dim_j": dj,
            "dim_k": dk,
            "norm_P": nrm,
            "d1": d1,
            "d2": d2,
            "dist_crys": dist_crys,
            "ipr": ipr,
            "n_eff": n_eff,
            "n_comp": n_comp,
            "max_comp": max_comp,
            "comp_sizes": comp_str,
        })
        if step in snap_set:
            # 64x64 射影強度マップ sqrt(sum_j |U[i,j,k]|^2)
            snaps_map_ik[f"map_ik_{step}"] = np.sqrt(np.sum(absU ** 2, axis=1))
            # 各要素 (0..63) のエネルギー (3 軸平均)
            e_i = np.sum(absU ** 2, axis=(1, 2))
            e_j = np.sum(absU ** 2, axis=(0, 2))
            e_k = np.sum(absU ** 2, axis=(0, 1))
            snaps_node_energy[f"node_e_{step}"] = (e_i + e_j + e_k) / 3.0
            # 非ゼロ (>1e-6) の振幅リスト
            snaps_abs_vals[f"abs_{step}"] = absU[absU > 1e-6]
            snaps_comp_members[step] = comp_members

    # step 0 の norm_P は 1 ステップ後の値を使うため先に計算
    P0 = tri_prod_fast(U)
    nrm0 = float(np.linalg.norm(P0))
    record(0, U, nrm0, 0.0, 0.0)

    for step in range(1, max_steps + 1):
        P = P0 if step == 1 else tri_prod_fast(U)
        nrm = float(np.linalg.norm(P))
        if nrm < 1e-12:
            print(f"step {step}: 0 に消滅 (||P||={nrm:.3e})")
            break
        U_next = P / nrm
        d1 = proj_dist(U, U_next)
        d2 = proj_dist(U_prev1, U_next) if step >= 2 else d1
        record(step, U_next, nrm, d1, d2)
        U_prev2 = U_prev1
        U_prev1 = U
        U = U_next

    # CSV 保存
    csv_path = os.path.join("results", "tri_loop_n64_512melt_timeseries.csv")
    fieldnames = list(rows[0].keys())
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    # NPZ 保存
    npz_path = os.path.join("results", "tri_loop_n64_512melt_snapshots.npz")
    save_dict = {"snap_steps": np.array(snap_steps, dtype=np.int32)}
    save_dict.update(snaps_map_ik)
    save_dict.update(snaps_node_energy)
    save_dict.update(snaps_abs_vals)
    np.savez_compressed(npz_path, **save_dict)

    # テキストサマリー保存
    lines = []
    lines.append("=" * 76)
    lines.append("N = 64 (O8 ⊗ O8) 実数 1 マス欠陥 (0,1,0)+1e-4: step 0〜500 厳密追跡サマリー")
    lines.append("=" * 76)
    show_steps = {0, 10, 20, 30, 35, 40, 50, 60, 70, 80, 100, 120, 150, 200, 300, 400, 500}
    for r in rows:
        if r["step"] in show_steps:
            lines.append(
                f"step {r['step']:3d}: nz(>1e-5)={r['nz_1e5']:5d}, N_eff={r['n_eff']:7.2f}, "
                f"軸=({r['dim_i']:2d},{r['dim_j']:2d},{r['dim_k']:2d}), ||P||={r['norm_P']:.6f}, "
                f"dist(T0)={r['dist_crys']:.4f}, d1={r['d1']:.2e}, 塊={r['n_comp']}個 [{r['comp_sizes']}]"
            )
    lines.append("\n--- スナップショット時点の連結成分メンバー ---")
    for st in snap_steps:
        if st in snaps_comp_members:
            cm = snaps_comp_members[st]
            lines.append(f"  step {st:3d}: {len(cm)} 個 -> {[f'size={len(c)}:{c[:8]}' + ('...' if len(c)>8 else '') for c in cm]}")
    lines.append(f"\n所要時間: {time.time() - t0:.2f} 秒")
    lines.append(f"保存ファイル: {csv_path}, {npz_path}")

    txt_path = os.path.join("results", "tri_loop_n64_512melt_run_v0.txt")
    with open(txt_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    run()
