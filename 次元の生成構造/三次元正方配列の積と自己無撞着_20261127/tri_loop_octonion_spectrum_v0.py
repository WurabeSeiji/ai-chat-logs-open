#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
テスト 2: 八元数積表 (8x8x8, 非ゼロ 64 マス) の全 512 マスすべてに 1 マス欠陥
(+ε および +iε, ε = 1e-4) を入れたときの溶け残りスペクトル (全 1,024 通りの全数分類)。

演算: 三角形ループ積 Tri(T)_{ijk} = sum_{p,q,r} T_{ipq} T_{qjr} T_{rpk}
乱数なし。||P|| < 1e-12 を厳密ゼロ停止とする。
"""
import os
import time
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


def classify_one(T0, T0_unit, T1_unit, pos, eps, max_steps=150):
    U = T0.copy()
    U[pos] += eps
    U = U / np.linalg.norm(U)
    hist = [U]
    for step in range(1, max_steps + 1):
        P = tri_prod_np(U)
        nrm = float(np.linalg.norm(P))
        if nrm < 1e-12:
            return ("zero", 0, 0.0, step)
        U_next = P / nrm
        # 自己修復チェック (元の周期 2 軌道のどちらかへ 1e-6 以内に戻ったか)
        if step >= 10:
            d_orig = min(proj_dist(T0_unit, U_next), proj_dist(T1_unit, U_next))
            if d_orig < 1e-6:
                return ("self_repair", 64, round(nrm, 4), step)
            # 他の固定点・周期軌道チェック (直近 16 ステップ)
            for prev_i in range(len(hist) - 1, max(-1, len(hist) - 17), -1):
                if proj_dist(hist[prev_i], U_next) < 1e-8:
                    p = step - prev_i
                    nz = int(np.sum(np.abs(U_next) > 1e-5))
                    return (f"periodic({p})", nz, round(nrm, 4), step)
        hist.append(U_next)
        U = U_next
    nz = int(np.sum(np.abs(U) > 1e-5))
    return ("alive_nonper", nz, round(nrm, 4), max_steps)


def run():
    t0 = time.time()
    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=" * 72)
    out("テスト 2: 八元数積表 (8^3) 全 512 マス 1 マス欠陥の溶け残りスペクトル (全数分類)")
    out("=" * 72)

    T0_int = make_mult_table(8)
    T0_real = np.array(T0_int, dtype=np.float64)
    T0_cplx = T0_real.astype(np.complex128)
    T1_real = tri_prod_np(T0_real)
    T0_unit_r = T0_real / np.linalg.norm(T0_real)
    T1_unit_r = T1_real / np.linalg.norm(T1_real)
    T0_unit_c = T0_unit_r.astype(np.complex128)
    T1_unit_c = T1_unit_r.astype(np.complex128)

    for mode_name, Tbase, U0, U1, eps in [
        ("実数欠陥 +ε (1e-4)", T0_real, T0_unit_r, T1_unit_r, 1e-4),
        ("虚数欠陥 +iε (1e-4j)", T0_cplx, T0_unit_c, T1_unit_c, 1e-4j),
    ]:
        out(f"\n--- [{mode_name}] 全 512 マス ---")
        # (元値が0か非ゼロか, 結果種別, nz, round(nrm, 2)) で集計
        summary = {}
        examples = {}
        for i in range(8):
            for j in range(8):
                for k in range(8):
                    orig_kind = "元非ゼロ" if T0_int[i][j][k] != 0 else "元ゼロ  "
                    kind, nz, nrm, st = classify_one(Tbase, U0, U1, (i, j, k), eps, max_steps=150)
                    key = (orig_kind, kind, nz, round(nrm, 2))
                    summary[key] = summary.get(key, 0) + 1
                    if key not in examples:
                        examples[key] = ((i, j, k), nrm, st)

        total_check = 0
        for key, cnt in sorted(summary.items(), key=lambda x: (x[0][0], -x[1])):
            orig_kind, kind, nz, nrm2 = key
            ex_pos, ex_nrm, ex_st = examples[key]
            out(f"  {orig_kind} | {kind:14s} | nz={nz:3d}/512 | ||P||≈{ex_nrm:.4f} | 件数={cnt:3d}  (例: {ex_pos}, step={ex_st})")
            total_check += cnt
        assert total_check == 512

    out(f"\n所要時間: {time.time() - t0:.2f} 秒")
    os.makedirs("results", exist_ok=True)
    out_path = os.path.join("results", "tri_loop_octonion_spectrum_v0.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"結果を保存しました: {out_path}")


if __name__ == "__main__":
    run()
