#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 1: 複素数 (2x2x2)・四元数 (4x4x4)・八元数 (8x8x8) の積表 f_{abc} そのものに対する
Bhattacharya-Mesner 積 (BM 積) の厳密整数計算。

定義 (三次元正方配列の積_勉強ノート.md §3.1):
    Prod(T, T, T)_{ijk} = sum_{t=0}^{N-1} T_{itk} * T_{ijt} * T_{tjk}

割り算なし・浮動小数なし・乱数なし。すべて厳密な整数で計算する。
あわせて、積表の 3 軸の並べ替え (6 通り) についても同じ計算を行い、
形を保って戻るか・0 に潰れるか・準周期になるかを確かめる。
"""
import itertools
import os
import sys

# --- ケーリー＝ディクソン構成による積表の厳密生成 ---

def cd_mul(x, y):
    """ケーリー＝ディクソン構成: (a, b)(c, d) = (ac - conj(d)*b, d*a + b*conj(c))"""
    n = len(x)
    if n == 1:
        return [x[0] * y[0]]
    h = n // 2
    a, b = x[:h], x[h:]
    c, d = y[:h], y[h:]
    def conj(v):
        return [v[0]] + [-val for val in v[1:]]
    ac = cd_mul(a, c)
    db = cd_mul(conj(d), b)
    da = cd_mul(d, a)
    bc = cd_mul(b, conj(c))
    left = [ac[i] - db[i] for i in range(h)]
    right = [da[i] + bc[i] for i in range(h)]
    return left + right


def make_mult_table(N):
    """N 次元 (N = 2, 4, 8) の積表 T[a][b][c] (e_a * e_b = sum_c T[a][b][c] e_c) を厳密整数で生成"""
    T = [[[0] * N for _ in range(N)] for _ in range(N)]
    for a in range(N):
        ea = [1 if i == a else 0 for i in range(N)]
        for b in range(N):
            eb = [1 if i == b else 0 for i in range(N)]
            prod = cd_mul(ea, eb)
            for c in range(N):
                T[a][b][c] = prod[c]
    # 検算: 各 (a, b) で非ゼロはちょうど 1 個 (±1)、合計 N^2 個
    nz = sum(1 for a in range(N) for b in range(N) for c in range(N) if T[a][b][c] != 0)
    assert nz == N * N, f"非ゼロ個数が不一致: {nz} != {N*N}"
    # 検算: e_0 が単位元、e_i^2 = -e_0 (i >= 1)
    for a in range(N):
        assert T[0][a][a] == 1 and T[a][0][a] == 1
        if a >= 1:
            assert T[a][a][0] == -1
    return T


def permute_axes(T, perm):
    """3 次元配列の軸を perm = (p0, p1, p2) で並べ替える"""
    N = len(T)
    U = [[[0] * N for _ in range(N)] for _ in range(N)]
    for i in range(N):
        for j in range(N):
            for k in range(N):
                idx = (i, j, k)
                U[idx[perm[0]]][idx[perm[1]]][idx[perm[2]]] = T[i][j][k]
    return U


def bm_prod(T):
    """Prod(T, T, T)_{ijk} = sum_t T[i][t][k] * T[i][j][t] * T[t][j][k] (厳密整数)"""
    N = len(T)
    R = [[[0] * N for _ in range(N)] for _ in range(N)]
    for i in range(N):
        for j in range(N):
            for k in range(N):
                s = 0
                for t in range(N):
                    s += T[i][t][k] * T[i][j][t] * T[t][j][k]
                R[i][j][k] = s
    return R


def stats(T):
    N = len(T)
    vals = [T[i][j][k] for i in range(N) for j in range(N) for k in range(N)]
    nz = sum(1 for v in vals if v != 0)
    max_abs = max(abs(v) for v in vals)
    sum_v = sum(vals)
    sum_sq = sum(v * v for v in vals)
    # 値の頻度
    freq = {}
    for v in vals:
        if v != 0:
            freq[v] = freq.get(v, 0) + 1
    return nz, max_abs, sum_v, sum_sq, dict(sorted(freq.items()))


def check_proportional(A, B):
    """B = lam * A となる整数 (または有理) 倍率 lam が存在するか厳密に判定"""
    N = len(A)
    ratio = None  # (num, den)
    for i in range(N):
        for j in range(N):
            for k in range(N):
                a, b = A[i][j][k], B[i][j][k]
                if a == 0 and b == 0:
                    continue
                if a == 0 and b != 0:
                    return None
                if ratio is None:
                    ratio = (b, a)
                else:
                    if b * ratio[1] != ratio[0] * a:
                        return None
    if ratio is None:
        return (0, 1)
    from math import gcd
    g = gcd(abs(ratio[0]), abs(ratio[1]))
    num, den = ratio[0] // g, ratio[1] // g
    if den < 0:
        num, den = -num, -den
    return (num, den)


def check_sum_z2_pairs(T):
    """
    実数配列 T に対し、第 3 軸の異なる 2 成分 (a, b) を組にして
    z_{ij} = T[i][j][a] + i * T[i][j][b] と読んだときの sum_{i,j} z_{ij}^2 の最大絶対値を調べる。
    実部 = sum(T_{ija}^2 - T_{ijb}^2)、虚部 = 2 * sum(T_{ija} * T_{ijb})。
    """
    N = len(T)
    nonzero_pairs = 0
    total_pairs = 0
    max_re = 0
    max_im = 0
    for a in range(N):
        for b in range(a + 1, N):
            total_pairs += 1
            re = 0
            im = 0
            for i in range(N):
                for j in range(N):
                    va, vb = T[i][j][a], T[i][j][b]
                    re += va * va - vb * vb
                    im += 2 * va * vb
            if re != 0 or im != 0:
                nonzero_pairs += 1
            max_re = max(max_re, abs(re))
            max_im = max(max_im, abs(im))
    return total_pairs, nonzero_pairs, max_re, max_im


def run():
    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=" * 72)
    out("Step 1: 複素数 (2^3)・四元数 (4^3)・八元数 (8^3) の積表に対する BM 積 (厳密整数)")
    out("=" * 72)

    for name, N in [("複素数 (C)", 2), ("四元数 (H)", 4), ("八元数 (O)", 8)]:
        out(f"\n--- {name}: N = {N} ---")
        T0 = make_mult_table(N)
        nz, max_abs, sum_v, sum_sq, freq = stats(T0)
        tp, nzp, mre, mim = check_sum_z2_pairs(T0)
        out(f"初期 T^(0): 非ゼロ={nz}/{N**3}, 最大絶対値={max_abs}, 和={sum_v}, 二乗和={sum_sq}, 分布={freq}")
        out(f"  2成分組の Σz^2 検算 (第3軸): {nzp}/{tp} 組が非ゼロ (max|Re|={mre}, max|Im|={mim})")

        # まず元の軸順 (0, 1, 2) で最大 4 ステップ追う
        T = T0
        history = [T0]
        for step in range(1, 5):
            T_next = bm_prod(T)
            nz, max_abs, sum_v, sum_sq, freq = stats(T_next)
            out(f"ステップ {step}: 非ゼロ={nz}/{N**3}, 最大絶対値={max_abs}, 和={sum_v}, 二乗和={sum_sq}")
            if len(freq) <= 10:
                out(f"  分布={freq}")
            else:
                out(f"  分布の種類数={len(freq)}")
            if nz == 0:
                out("  → 0 に完全消滅。")
                break
            # 過去との比例チェック
            found_rel = False
            for prev_idx, prev_T in enumerate(history):
                prop = check_proportional(prev_T, T_next)
                if prop is not None:
                    num, den = prop
                    lam_str = f"{num}" if den == 1 else f"{num}/{den}"
                    out(f"  → T^({step}) = ({lam_str}) * T^({prev_idx})  【厳密一致】")
                    found_rel = True
            tp, nzp, mre, mim = check_sum_z2_pairs(T_next)
            out(f"  2成分組の Σz^2 検算 (第3軸): {nzp}/{tp} 組が非ゼロ (max|Re|={mre}, max|Im|={mim})")
            history.append(T_next)
            T = T_next
            if found_rel:
                break

        # 6 通りの軸の並べ替えで 1 ステップ目の挙動をすべて確認
        out(f"  [軸の並べ替え 6 通りの 1 ステップ目確認]")
        for perm in itertools.permutations((0, 1, 2)):
            P0 = permute_axes(T0, perm)
            P1 = bm_prod(P0)
            nz1, max1, sum1, sq1, freq1 = stats(P1)
            prop = check_proportional(P0, P1)
            prop_str = "なし" if prop is None else (f"{prop[0]}" if prop[1] == 1 else f"{prop[0]}/{prop[1]}")
            out(f"    perm={perm}: 非ゼロ={nz1}, max={max1}, 比例倍率={prop_str}, 分布={freq1 if len(freq1)<=6 else f'{len(freq1)}種'}")

    os.makedirs("results", exist_ok=True)
    out_path = os.path.join("results", "bm_product_algebra_tables_v0.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\n結果を保存しました: {out_path}")


if __name__ == "__main__":
    run()
