#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2: 巡回対称な反対称平面波 F_{a,b}(u, v) = ω^{a u + b v} - ω^{a v + b u} の
解析定理の厳密検証と、前回 N = 8 で出た「non_periodic」軌道の正体検証。

解析定理:
  a != b (mod N) に対して
    Prod(F_{a,b})(u, v) = N * F_{2a, 2b}(u, v)
  が厳密に成り立つ。
  - N が奇数なら、2 の乗法位数 ord_N(2) (またはその約数) で厳密に周期的 (U^n = I)。
  - N が 2 の冪 (4, 8, ...) なら、高々 log_2(N) ステップで (2^k a, 2^k b) ≡ (0, 0) mod N
    となり、厳密には 0 に消滅する。

あわせて、前回 N = 4, 5, 8 で non_periodic と出た他の軌道 (2 次位相、R4 のみ、
cos+i sin) について、正規化前ノルム ||P|| の履歴と 200 ステップ追跡により、
「丸めノイズの増幅」「遅い収束」「真の準周期 / カオス」を完全に切り分ける。
"""
import math
import os
import numpy as np


def bm_circulant(F):
    N = F.shape[0]
    F_next = np.zeros_like(F)
    for u in range(N):
        for v in range(N):
            s_val = 0j
            for s in range(N):
                s_val += F[s, v] * F[u, s] * F[(u - s) % N, (v - s) % N]
            F_next[u, v] = s_val
    return F_next


def make_antisym(N, a, b):
    w = np.exp(2j * np.pi / N)
    return np.array(
        [[w ** ((a * u + b * v) % N) - w ** ((a * v + b * u) % N) for v in range(N)] for u in range(N)],
        dtype=np.complex128,
    )


def proj_dist(A, B):
    na = np.linalg.norm(A)
    nb = np.linalg.norm(B)
    if na == 0 or nb == 0:
        return 1.0
    ov = abs(np.vdot(A / na, B / nb))
    return float(np.sqrt(max(0.0, 1.0 - ov * ov)))


def sum_z2_ratio(F):
    denom = float(np.sum(np.abs(F) ** 2))
    return 0.0 if denom == 0 else float(abs(np.sum(F ** 2)) / denom)


def mult_order_2(N):
    """奇数 N に対する 2 の乗法位数 ord_N(2)"""
    k = 1
    val = 2 % N
    while val != 1:
        val = (val * 2) % N
        k += 1
    return k


def run():
    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=" * 72)
    out("2: 反対称波の解析定理の厳密検証 と N=4,5,8 非周期軌道の正体検証")
    out("=" * 72)

    out("\n--- 1. 解析定理 Prod(F_{a,b}) = N * F_{2a, 2b} の数値検算 ---")
    for N in (3, 4, 5, 7, 8, 9, 11, 13):
        a, b = 1, 2
        F = make_antisym(N, a, b)
        P = bm_circulant(F)
        F2 = make_antisym(N, (2 * a) % N, (2 * b) % N)
        err = float(np.max(np.abs(P - N * F2)))
        sz2 = sum_z2_ratio(F)
        if N % 2 == 1:
            ord2 = mult_order_2(N)
            # 実際に ord2 ステップ回して元に戻るか確認
            cur = F.copy()
            for _ in range(ord2):
                cur = bm_circulant(cur) / N
            err_per = float(np.max(np.abs(cur - F)))
            out(f"  N={N:2d} (奇数): max|Prod(F_{{1,2}}) - N*F_{{2,4}}| = {err:.2e}, 周期 ord_N(2)={ord2:2d} 戻り誤差={err_per:.2e}, |ΣF^2|/Σ|F|^2={sz2:.2e}")
        else:
            # N が偶数: 2 の冪乗回で消える様子を確認
            cur = F.copy()
            norms = []
            for st in range(1, 6):
                cur = bm_circulant(cur)
                norms.append(float(np.linalg.norm(cur)))
            out(f"  N={N:2d} (偶数): max|Prod(F_{{1,2}}) - N*F_{{2,4}}| = {err:.2e}, 非正規化 ||T^(k)|| (k=1..5) = {[f'{x:.2e}' for x in norms]}")

    out("\n  [訂正の確認] 前回 N=8 の 反対称 ω^(1u+2v)-対称 が non_periodic と出た原因:")
    F8 = make_antisym(8, 1, 2)
    U = F8 / np.linalg.norm(F8)
    for st in range(1, 5):
        P = bm_circulant(U)
        nrm = float(np.linalg.norm(P))
        out(f"    step {st}: ||Prod(U)|| = {nrm:.3e}")
        U = P / nrm
    out("  → step 3 で厳密には 0 (数値上 ~1e-16 の丸め誤差) になっていたものを正規化して増幅した偽の軌道だった。")

    out("\n--- 2. 前回 non_periodic と出た他の軌道の正体検証 (200 ステップ、||P||<1e-12 停止) ---")
    cases = []
    # (a) N=4 2次位相 p=0, q=1, r=0
    w4 = np.exp(2j * np.pi / 4)
    cases.append(("N=4 2次位相 (p=0,q=1,r=0)", np.array([[w4 ** ((v * v) % 4) for v in range(4)] for u in range(4)], dtype=np.complex128)))
    # (b) N=8 2次位相 p=0, q=1, r=0
    w8 = np.exp(2j * np.pi / 8)
    cases.append(("N=8 2次位相 (p=0,q=1,r=0)", np.array([[w8 ** ((v * v) % 8) for v in range(8)] for u in range(8)], dtype=np.complex128)))
    # (c) N=8 R4のみ ω^(1uv)
    F_r4 = np.zeros((8, 8), dtype=np.complex128)
    for u in range(8):
        for v in range(8):
            if u != 0 and v != 0 and u != v:
                F_r4[u, v] = w8 ** ((u * v) % 8)
    cases.append(("N=8 R4のみ ω^(1uv)", F_r4))
    # (d) N=5 cos(u) + i sin(v)
    cases.append(("N=5 cos(u)+i sin(v)", np.array([[np.cos(2 * np.pi * u / 5) + 1j * np.sin(2 * np.pi * v / 5) for v in range(5)] for u in range(5)], dtype=np.complex128)))
    # (e) N=8 cos(u) + i sin(v)
    cases.append(("N=8 cos(u)+i sin(v)", np.array([[np.cos(2 * np.pi * u / 8) + 1j * np.sin(2 * np.pi * v / 8) for v in range(8)] for u in range(8)], dtype=np.complex128)))

    for label, F0 in cases:
        out(f"\n  [{label}]")
        U = F0 / np.linalg.norm(F0)
        hist = [U]
        for st in range(1, 201):
            P = bm_circulant(U)
            nrm = float(np.linalg.norm(P))
            if nrm < 1e-12:
                out(f"    step {st:3d}: ||P||={nrm:.3e} → 厳密ゼロ (丸めノイズ)")
                break
            U_next = P / nrm
            d_prev = proj_dist(U, U_next)
            d_prev2 = proj_dist(hist[-2], U_next) if len(hist) >= 2 else 1.0
            sz2 = sum_z2_ratio(U_next)
            if st in (1, 10, 40, 100, 200):
                out(f"    step {st:3d}: dist(prev1)={d_prev:.3e}, dist(prev2)={d_prev2:.3e}, |ΣF^2|/Σ|F|^2={sz2:.3e}, ||P||={nrm:.4f}")
            hist.append(U_next)
            U = U_next
        else:
            # 後半 50 ステップの最小再帰距離
            min_rec = min(proj_dist(h, U) for h in hist[-51:-1])
            out(f"    200 step 後: dist(prev1)={proj_dist(hist[-2], U):.3e}, min_rec_50={min_rec:.3e}")

    os.makedirs("results", exist_ok=True)
    out_path = os.path.join("results", "bm_product_antisym_wave_v0.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\n結果を保存しました: {out_path}")


if __name__ == "__main__":
    run()
