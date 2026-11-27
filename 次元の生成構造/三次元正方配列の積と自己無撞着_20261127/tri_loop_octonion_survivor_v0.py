#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
三角形ループ積 Tri(T)_{ijk} = sum_{p,q,r} T_{ipq} T_{qjr} T_{rpk} において、
四元数 (4^3)・八元数 (8^3) の積表に 1 マス欠陥を入れたときに出現した
「non_periodic (nz = 8)」の溶け残り軌道を 500 ステップ追跡し、正体を確認する。

確認項目:
  1. 500 ステップ後も 0 に消滅せず回り続けるか (||P|| の推移と最小再帰距離)。
  2. 生き残った 8 マスの位置と値・位相の動き (準周期か周期かカオスか)。
  3. 他の代表位置 (虚数単位どうしの積位置、ゼロ位置など) に 1 マス欠陥を入れたときの挙動。
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


def run():
    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=" * 72)
    out("四元数 (4^3)・八元数 (8^3) 積表 + 三角形ループ積: 溶け残り (nz=8) の 500 ステップ検証")
    out("=" * 72)

    for name, N in [("四元数 (H)", 4), ("八元数 (O)", 8)]:
        out(f"\n=== {name}: N = {N} ===")
        T0 = np.array(make_mult_table(N), dtype=np.complex128)
        T0_unit = T0 / np.linalg.norm(T0)

        defects = [
            ("(0,1,0)+ε [元0位置]", (0, 1, 0), 1e-4),
            ("(0,1,0)+iε [元0位置]", (0, 1, 0), 1e-4j),
            ("(0,1,1)+ε [元+1位置]", (0, 1, 1), 1e-4),
            ("(1,1,0)+ε [元-1位置]", (1, 1, 0), 1e-4),
            ("(1,2,3)+ε [虚数積位置]", (1, 2, 3), 1e-4),
            ("(1,2,3)+iε [虚数積位置]", (1, 2, 3), 1e-4j),
            ("(1,2,0)+ε [元0位置]", (1, 2, 0), 1e-4),
        ]

        for label, pos, eps in defects:
            U = T0.copy()
            U[pos] += eps
            U = U / np.linalg.norm(U)
            hist = [U]
            status = "running"
            for step in range(1, 501):
                P = tri_prod_np(U)
                nrm = float(np.linalg.norm(P))
                if nrm < 1e-12:
                    status = f"zero at step {step}"
                    break
                U_next = P / nrm
                # 周期チェック (step >= 50)
                if step >= 50:
                    for prev_i in range(len(hist) - 1, max(-1, len(hist) - 31), -1):
                        if proj_dist(hist[prev_i], U_next) < 1e-9:
                            p = step - prev_i
                            nz = int(np.sum(np.abs(U_next) > 1e-5))
                            d0 = proj_dist(T0_unit, U_next)
                            status = f"periodic({p}) at step {step} (dist_T0={d0:.4f}, nz={nz}/{N**3}, ||P||={nrm:.4f})"
                            U = U_next
                            break
                    if status != "running":
                        break
                hist.append(U_next)
                U = U_next

            if status == "running":
                nz = int(np.sum(np.abs(U) > 1e-5))
                d0 = proj_dist(T0_unit, U)
                d_prev = proj_dist(hist[-2], U)
                min_rec = min(proj_dist(h, U) for h in hist[-101:-1])
                status = f"ALIVE at step 500 (dist_T0={d0:.4f}, d_prev={d_prev:.3f}, min_rec100={min_rec:.3e}, nz={nz}/{N**3}, ||P||={nrm:.4f})"
            out(f"  [{label:24s}] → {status}")

            # もし step 500 で ALIVE かつ nz <= 16 なら、その位置と値を表示 (最初の 1 例だけ詳細)
            if status.startswith("ALIVE") and label.startswith("(0,1,0)+ε"):
                nz_pos = [(i, j, k, U[i, j, k]) for i in range(N) for j in range(N) for k in range(N) if abs(U[i, j, k]) > 1e-5]
                out(f"    [step 500 の非ゼロ成分 ({len(nz_pos)} 個)]:")
                for (i, j, k, val) in nz_pos[:16]:
                    out(f"      ({i},{j},{k}): {val.real:+.6f} {val.imag:+.6f}i  (|v|={abs(val):.6f})")
                # 直近 10 ステップの ||P|| と d_prev の変動を確認
                d_recent = [round(proj_dist(hist[-m-1], hist[-m]), 4) for m in range(1, 7)]
                out(f"    直近 6 ステップの d_prev: {d_recent}")
                # リャプノフ指数チェック (step 500 の状態に 1e-7 摂動を入れて 20 ステップ追跡)
                V = U.copy()
                V[0, 0, 0] += 1e-7
                V = V / np.linalg.norm(V)
                U_ref = U.copy()
                dists = []
                for k_st in range(1, 21):
                    U_ref = tri_prod_np(U_ref)
                    U_ref = U_ref / np.linalg.norm(U_ref)
                    V = tri_prod_np(V)
                    V = V / np.linalg.norm(V)
                    dists.append(proj_dist(U_ref, V))
                out(f"    step 500 での微小摂動 1e-7 の 1,5,10,20 step 後距離: {[f'{dists[m]:.2e}' for m in (0, 4, 9, 19)]}")

    os.makedirs("results", exist_ok=True)
    out_path = os.path.join("results", "tri_loop_octonion_survivor_v0.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\n結果を保存しました: {out_path}")


if __name__ == "__main__":
    run()
