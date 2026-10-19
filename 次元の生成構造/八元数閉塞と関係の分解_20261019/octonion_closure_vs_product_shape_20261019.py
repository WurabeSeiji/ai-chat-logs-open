#!/usr/bin/env python3
"""
閉塞の形と積の形（2026-10-19、Claude Code の指摘を受けた追加計算）

`octonion_closure_check_20261019.py` と同じフォルダに置いて実行する。

問い：results.md の次元表の「形の空間」は、何で割った形か。

答え：閉塞の対称性 SO(7) × ℝ⁺（22 次元）で割ったもの。しかし SO(7) は積の対称性ではない
（積を保つのは G₂、14 次元）。読み出し・ホロノミー・関係の分解はすべて積を使うので、
それらが実際に見る形の空間は G₂ × ℝ⁺（15 次元）で割ったもので、こちらの方が大きい。

計算で出た結果：N ≤ 4 では二つが一致し、N = 5 で初めて分かれる（8 と 9）。
N = 5 が閾値だというのは「閉塞の見方と積の見方が初めて食い違う N」という意味。
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, "octonion_closure_check_20261019.py")).read().split("results = {")[0])

def rank_(A, tol=1e-9):
    s = np.linalg.svd(A, compute_uv=False)
    return int(np.sum(s > tol * s[0]))

# g₂ の基底（微分の方程式の核、14 次元）
rows = []
for i in range(8):
    for j in range(8):
        prod = omul(E[i], E[j])
        for c in range(8):
            row = np.zeros(64)
            for m in range(8): row[c*8 + m] += prod[m]
            for m in range(8): row[m*8 + i] -= omul(E[m], E[j])[c]
            for m in range(8): row[m*8 + j] -= omul(E[i], E[m])[c]
            rows.append(row)
u, s, vt = np.linalg.svd(np.array(rows))
g2 = [v.reshape(8, 8) for v in vt[int(np.sum(s > 1e-9 * s[0])):]]

# so(7) の基底（虚部の回転、21 次元）＝ 閉塞の対称性のうち回転の部分
so7 = []
for i in range(1, 8):
    for j in range(i + 1, 8):
        D = np.zeros((8, 8)); D[i, j] = -1; D[j, i] = 1; so7.append(D)
print(f"dim g₂ = {len(g2)}（積の対称性）, dim so(7) = {len(so7)}（閉塞の対称性のうち回転）")

def gen_closure(N, r):
    one = np.ones(N) / np.sqrt(N)
    x0 = r.normal(size=N); x0 -= one * (one @ x0)
    P = np.eye(N) - np.outer(one, one) - np.outer(x0, x0) / (x0 @ x0)
    X = np.array([P @ r.normal(size=N) for _ in range(7)])
    X *= np.sqrt((x0 @ x0) / np.sum(X * X))
    return np.column_stack([x0, X.T])

def constraint_rank(S):
    N = S.shape[0]; J = np.zeros((16, N * 8))
    for n in range(N):
        for k in range(8):
            J[:8, n*8 + k] = E[k]; J[8:, n*8 + k] = 2.0 * ojordan(S[n], E[k])
    return rank_(J)

def orbit_dim(S, alg):
    """スケール ＋ alg（行列のリスト）が配置に作る軌道の次元。"""
    T = [S.reshape(-1)] + [(S @ D.T).reshape(-1) for D in alg]
    return rank_(np.array(T))

print("\n N  解の次元  SO(7)×ℝ 軌道  形（閉塞のみ）  G₂×ℝ 軌道  形（積も見る）  開く虚方向 k")
for N in range(3, 11):
    S = gen_closure(N, np.random.default_rng(4000 + N))
    sol = N * 8 - constraint_rank(S)
    o7 = orbit_dim(S, so7); o2 = orbit_dim(S, g2); k = rank_(S[:, 1:].T)
    print(f"{N:>2}    {sol:>5}       {o7:>6}        {sol-o7:>6}       {o2:>6}       {sol-o2:>6}         {k}")

print("\n複数シードでの確認（各 12 本）")
for N in [4, 5, 6]:
    sols, o7s, o2s = [], [], []
    for t in range(12):
        S = gen_closure(N, np.random.default_rng(7000 + 100*N + t))
        sols.append(N * 8 - constraint_rank(S)); o7s.append(orbit_dim(S, so7)); o2s.append(orbit_dim(S, g2))
    sols = np.array(sols); o7s = np.array(o7s); o2s = np.array(o2s)
    print(f"  N = {N}: 解 {set(sols.tolist())}, 形（閉塞のみ）{set((sols-o7s).tolist())}, "
          f"形（積も見る）{set((sols-o2s).tolist())}")

print("\n四元数部分代数の中に作った N = 5 解との比較")
r = np.random.default_rng(7)
A = r.normal(size=8); A[0] = 0; A /= np.linalg.norm(A)
B = r.normal(size=8); B[0] = 0; B -= (A @ B) * A; B /= np.linalg.norm(B)
C = omul(A, B); C[0] = 0; C /= np.linalg.norm(C)
Bas = np.array([E[0], A, B, C])
N = 5; one = np.ones(N) / np.sqrt(N)
x0 = r.normal(size=N); x0 -= one * (one @ x0)
P = np.eye(N) - np.outer(one, one) - np.outer(x0, x0) / (x0 @ x0)
X = np.array([P @ r.normal(size=N) for _ in range(3)]); X *= np.sqrt((x0 @ x0) / np.sum(X * X))
S = np.column_stack([x0, X.T]) @ Bas
o2 = orbit_dim(S, g2)
print(f"  四元数に乗る N = 5 解: G₂×ℝ 軌道 = {o2}（一般の解は 15）")
print(f"  四元数に乗る解の族は 12 + 8 = 20 次元（24 の中で余次元 4）。その族の中での形 = 20 − {o2} = {20-o2}")
print(f"  一般の解の形 9 との差は 1。残る 3 は、四元数に乗る解が持つ余分な固定部分群（内部 su(2)）に食われている。")
