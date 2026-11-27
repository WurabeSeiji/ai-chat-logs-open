#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
三次元正方配列 (立方体) の二つの積を、十字 N=4 (z_k = i^k) で手計算どおりに確かめる。勉強ノートの検算用。

  二項の積 (t 積, Kilmer-Martin 2011 定義 3.1):
      (A*B)[i,j,s] = Σ_{a+b=s mod n} Σ_t A[i,t,a] B[t,j,b]
  三項の積 (Bhattacharya-Mesner 積, Gnang arXiv:1411.6270 定義 1 の三階の場合):
      Prod(A,B,C)[i,j,k] = Σ_t A[i,t,k] B[i,j,t] C[t,j,k]

乱数・割り算・正規化は使わない。成分はガウス整数 (整数+整数 i) だけで、complex128 で厳密に表せる大きさに収まる。
使い方: python3 cubic_products_N4_cross_v0.py
"""
import numpy as np

n = 4
z = np.array([1, 1j, -1, -1j])          # 十字: z_k = i^k


def tprod(A, B):
    C = np.zeros_like(A)
    for a in range(n):
        for b in range(n):
            C[:, :, (a + b) % n] += A[:, :, a] @ B[:, :, b]
    return C


def bm(A, B, C):
    return np.einsum('itk,ijt,tjk->ijk', A, B, C)


def t_identity():
    I = np.zeros((n, n, n), complex)
    I[:, :, 0] = np.eye(n)
    return I


def show_nonzero(X, name):
    idx = np.argwhere(np.abs(X) > 0)
    print("  %s: 0 でない成分 %d 個" % (name, len(idx)))
    for i, j, k in idx[:16]:
        v = X[i, j, k]
        print("     (%d,%d,%d) = %s" % (i, j, k, fmt(v)))


def fmt(v):
    re, im = int(round(v.real)), int(round(v.imag))
    if im == 0:
        return "%d" % re
    if re == 0:
        return "%di" % im
    return "%d%+di" % (re, im)


print("=== 0. 十字の確認 ===")
print("  z =", [fmt(v) for v in z], " Σz =", fmt(z.sum()), " Σz^2 =", fmt((z**2).sum()),
      " Σz^3 =", fmt((z**3).sum()), " Σz^4 =", fmt((z**4).sum()))

print("\n=== 1. 対角だけに置く: O[k,k,k] = z_k ===")
O = np.zeros((n, n, n), complex)
for k in range(n):
    O[k, k, k] = z[k]
show_nonzero(O, "O")
print(" [二項の積] O, O*O, O*O*O, O*O*O*O")
P = O.copy()
for p in range(2, 9):
    P = tprod(P, O)
    isI = np.array_equal(P, t_identity())
    isO = np.array_equal(P, O)
    print("  O^%d : 単位の立方体に等しい=%s  O に等しい=%s" % (p, isI, isO))
    if p <= 4:
        show_nonzero(P, "O^%d" % p)
print(" [三項の積] Prod(O,O,O) を繰り返す")
Q = O.copy()
for p in range(1, 4):
    Q = bm(Q, Q, Q)
    print("  %d 回目: O に等しい=%s  O の共役に等しい=%s  対角以外は 0=%s"
          % (p, np.array_equal(Q, O), np.array_equal(Q, O.conj()),
             all(abs(Q[i, j, k]) == 0 for i in range(n) for j in range(n) for k in range(n) if not (i == j == k))))

print("\n=== 2. 三つの積を全部並べる: X[i,j,k] = z_i z_j z_k ===")
X = np.einsum('i,j,k->ijk', z, z, z)
vals, cnt = np.unique(np.round(X, 9), return_counts=True)
print("  成分の値と個数:", {fmt(v): int(c) for v, c in zip(vals, cnt)}, " ΣX =", fmt(X.sum()), " ΣX^2 =", fmt((X**2).sum()))
print("  [二項の積] X*X の成分の最大の大きさ =", float(np.abs(tprod(X, X)).max()))
print("  [三項の積] Prod(X,X,X) の成分の最大の大きさ =", float(np.abs(bm(X, X, X)).max()))

print("\n=== 3. 真ん中の添字だけ共役: Y[i,j,k] = z_i conj(z_j) z_k ===")
Y = np.einsum('i,j,k->ijk', z, z.conj(), z)
YY = tprod(Y, Y)
print("  [二項の積] Y*Y = 16 Y か:", np.array_equal(YY, 16 * Y))
B3 = bm(Y, Y, Y)
r = B3[np.abs(Y) > 0] / Y[np.abs(Y) > 0]
print("  [三項の積] Prod(Y,Y,Y) の最大の大きさ =", float(np.abs(B3).max()),
      " Y の定数倍か:", bool(np.allclose(r, r[0])) if np.abs(B3).max() > 0 else "0 なので該当なし")

print("\n=== 4. 三項の積の性質 (決め打ちの整数の配列で確認。乱数ではない) ===")
i, j, k = np.meshgrid(range(n), range(n), range(n), indexing='ij')
A = (1 + i + 2 * j * j + 3 * k * k * k + i * j * k).astype(complex)
Bq = (2 - i * i + j + k * k + 2 * i * k).astype(complex)
Cq = (1 + i * j - 2 * k + j * j * k).astype(complex)
cyc = lambda T: np.transpose(T, (1, 2, 0))          # 軸を一つ回す: T'[j,k,i] = T[i,j,k]
lhs = cyc(bm(A, Bq, Cq))
cands = {"(A',B',C')": bm(cyc(A), cyc(Bq), cyc(Cq)), "(B',C',A')": bm(cyc(Bq), cyc(Cq), cyc(A)),
         "(C',A',B')": bm(cyc(Cq), cyc(A), cyc(Bq))}
for name, rhs in cands.items():
    print("  軸を回した Prod(A,B,C) = Prod%s か: %s" % (name, np.array_equal(lhs, rhs)))
rev = lambda T: np.transpose(T, (2, 1, 0))          # 1 番目と 3 番目の軸の入れ替え
lhs2 = rev(bm(A, Bq, Cq))
for name, rhs in {"(A~,B~,C~)": bm(rev(A), rev(Bq), rev(Cq)), "(C~,B~,A~)": bm(rev(Cq), rev(Bq), rev(A)),
                  "(B~,A~,C~)": bm(rev(Bq), rev(A), rev(Cq)), "(A~,C~,B~)": bm(rev(A), rev(Cq), rev(Bq))}.items():
    print("  軸 1,3 を入れ替えた Prod(A,B,C) = Prod%s か: %s" % (name, np.array_equal(lhs2, rhs)))
P3 = bm(A, A, A)
c1, c2, c3 = bm(P3, A, A), bm(A, P3, A), bm(A, A, P3)
print("  結合性: Prod(Prod(A,A,A),A,A) = Prod(A,Prod(A,A,A),A) か: %s / = Prod(A,A,Prod(A,A,A)) か: %s"
      % (np.array_equal(c1, c2), np.array_equal(c1, c3)))
print("  参考 (二項の積): (A*A)*A = A*(A*A) か: %s" % np.array_equal(tprod(tprod(A, A), A), tprod(A, tprod(A, A))))
