#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
テスト 1:
  (A) N = 16 ケーリー＝ディクソン十六元数 (Sedenion, 16^3 = 4096 成分, 非ゼロ 256 マス) の
      積表に対する三角形ループ積と 1 マス欠陥発展。
  (B) 八元数積表を 2 個並べた 16x16x16 ブロック結晶 (ブロック 0: 0..7, ブロック 1: 8..15) に
      対し、各ブロック内欠陥およびブロック間結合を入れたときの多粒子 (複数溶け残り) 検証。

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


def track_orbit(T0, max_steps=150):
    U = T0 / np.linalg.norm(T0)
    hist = [U]
    for step in range(1, max_steps + 1):
        P = tri_prod_np(U)
        nrm = float(np.linalg.norm(P))
        if nrm < 1e-12:
            return f"zero at step {step}", step, 0, 0.0, U
        U_next = P / nrm
        if step >= 15:
            for prev_i in range(len(hist) - 1, max(-1, len(hist) - 21), -1):
                if proj_dist(hist[prev_i], U_next) < 1e-8:
                    p = step - prev_i
                    nz = int(np.sum(np.abs(U_next) > 1e-5))
                    return f"periodic({p}) at step {step}", step, nz, nrm, U_next
        hist.append(U_next)
        U = U_next
    nz = int(np.sum(np.abs(U) > 1e-5))
    min_rec = min(proj_dist(h, U) for h in hist[-41:-1])
    return f"ALIVE at step {max_steps} (min_rec40={min_rec:.2e})", max_steps, nz, nrm, U


def count_block_nz(U, split=8, tol=1e-5):
    """16x16x16 をブロック 0 (0..7) とブロック 1 (8..15) とブロック間に分けて非ゼロ数を数える"""
    mask = np.abs(U) > tol
    b0 = int(np.sum(mask[:split, :split, :split]))
    b1 = int(np.sum(mask[split:, split:, split:]))
    tot = int(np.sum(mask))
    cross = tot - b0 - b1
    return b0, b1, cross, tot


def run():
    t0 = time.time()
    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=" * 72)
    out("テスト 1: 十六元数 (N=16) 積表 と 2 ブロック八元数 (N=16) の三角形ループ積検証")
    out("=" * 72)

    # --- (A) 十六元数 (N = 16) の積表 ---
    out("\n--- (A) 十六元数 (Sedenion, N = 16, 4096 成分, 初期非ゼロ 256) ---")
    S0 = np.array(make_mult_table(16), dtype=np.float64)
    S1 = tri_prod_np(S0)
    S2 = tri_prod_np(S1)
    nz1 = int(np.sum(np.abs(S1) > 1e-9))
    nz2 = int(np.sum(np.abs(S2) > 1e-9))
    vals1 = sorted({int(round(v)) for v in S1.flat if abs(v) > 1e-9})
    out(f"  欠陥なし step 1: nz={nz1}/4096, 値の集合={vals1}")
    # S2 が S0 に比例するか確認
    ratio_s2 = None
    if nz2 > 0:
        d_s2_s0 = proj_dist(S0, S2)
        out(f"  欠陥なし step 2: nz={nz2}/4096, dist(S2, S0)={d_s2_s0:.3e}, max|S2|={np.max(np.abs(S2)):.1f}")
    st_s, _, nz_s, nrm_s, _ = track_orbit(S0, max_steps=100)
    out(f"  欠陥なし 100 step 軌道: {st_s}, nz={nz_s}/4096, ||P||={nrm_s:.4f}")

    # 十六元数の 1 マス欠陥 6 通り
    S0_unit = S0 / np.linalg.norm(S0)
    sed_defects = [
        ("元0 (0,1,0)+ε", (0, 1, 0), 1e-4, False),
        ("元0 (0,1,0)+iε", (0, 1, 0), 1e-4j, True),
        ("元0 (1,2,0)+ε", (1, 2, 0), 1e-4, False),
        ("元0 (0,9,0)+ε [後半軸]", (0, 9, 0), 1e-4, False),
        ("元非ゼロ (0,1,1)+ε", (0, 1, 1), 1e-4, False),
        ("元非ゼロ (1,9,8)+ε", (1, 9, 8), 1e-4, False),
    ]
    for label, pos, eps, is_c in sed_defects:
        Sd = S0.astype(np.complex128) if is_c else S0.copy()
        Sd[pos] += eps
        st, _, nz, nrm, Uend = track_orbit(Sd, max_steps=150)
        d0 = proj_dist(S0_unit, Uend)
        out(f"  [{label:22s}] → {st}, dist(S0)={d0:.4f}, nz={nz}/4096, ||P||={nrm:.4f}")

    # --- (B) 2 ブロック八元数 (N = 16: ブロック 0 = 0..7, ブロック 1 = 8..15) ---
    out("\n--- (B) 2 ブロック八元数 (16x16x16: ブロック0=0..7, ブロック1=8..15, 初期非ゼロ 128) ---")
    O8 = np.array(make_mult_table(8), dtype=np.float64)
    B0 = np.zeros((16, 16, 16), dtype=np.float64)
    B0[:8, :8, :8] = O8
    B0[8:, 8:, 8:] = O8
    B2 = tri_prod_np(tri_prod_np(B0))
    assert np.allclose(B2, 1296.0 * B0, atol=1e-10), "2ブロック八元数の周期2検算失敗"
    out("[検算合格] 2 ブロック八元数は厳密に周期 2 (B2 = 1296 * B0, 誤差 < 1e-10)。")

    multi_cases = [
        ("各ブロックに同サイズ欠陥 (0,1,0)+1e-4, (8,9,8)+1e-4", [((0, 1, 0), 1e-4), ((8, 9, 8), 1e-4)]),
        ("各ブロックに微差欠陥 (0,1,0)+1e-4, (8,9,8)+1.01e-4", [((0, 1, 0), 1e-4), ((8, 9, 8), 1.01e-4)]),
        ("ブロック0のみ欠陥 (0,1,0)+1e-4 (ブロック1は無傷)", [((0, 1, 0), 1e-4)]),
        ("ブロック間結合1マスのみ (0,8,0)+1e-4", [((0, 8, 0), 1e-4)]),
        ("ブロック間結合対称 (0,8,0)+1e-4, (8,0,8)+1e-4", [((0, 8, 0), 1e-4), ((8, 0, 8), 1e-4)]),
        ("各ブロック欠陥+ブロック間結合 (0,1,0),(8,9,8),(0,8,0) 各1e-4", [((0, 1, 0), 1e-4), ((8, 9, 8), 1e-4), ((0, 8, 0), 1e-4)]),
    ]
    for label, dlist in multi_cases:
        Bd = B0.copy()
        for pos, val in dlist:
            Bd[pos] += val
        st, _, nz, nrm, Uend = track_orbit(Bd, max_steps=150)
        b0_nz, b1_nz, cr_nz, tot_nz = count_block_nz(Uend)
        out(f"  [{label}]")
        out(f"    → {st}, ||P||={nrm:.4f}, 非ゼロ(B0={b0_nz}, B1={b1_nz}, 結合={cr_nz}, 計={tot_nz}/4096)")

    out(f"\n所要時間: {time.time() - t0:.2f} 秒")
    os.makedirs("results", exist_ok=True)
    out_path = os.path.join("results", "tri_loop_sedenion_and_multi_v0.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"結果を保存しました: {out_path}")


if __name__ == "__main__":
    run()
