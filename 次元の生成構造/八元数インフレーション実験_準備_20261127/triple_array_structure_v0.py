#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
3 個の八元数の積の配列 X_ijk (3^3 = 27 個) そのものの構造。v0。
演算の繰り返しも、割り算も、倍音も入れない。状態 3 個と、その積だけを調べる。
仕様は「思考実験.md」§17。
使い方: python3 triple_array_structure_v0.py
"""
import itertools
import numpy as np

rng = np.random.default_rng(20261127)
TRI = [(1, 2, 3), (1, 4, 5), (1, 7, 6), (2, 4, 6), (2, 5, 7), (3, 4, 7), (3, 6, 5)]
f = np.zeros((8, 8, 8))
for a, b, c in TRI:
    for (i, j, k), s in [((a, b, c), 1), ((b, c, a), 1), ((c, a, b), 1), ((b, a, c), -1), ((a, c, b), -1), ((c, b, a), -1)]:
        f[i, j, k] = s
T = np.zeros((8, 8, 8)); T[0, 0, 0] = 1
for i in range(1, 8):
    T[0, i, i] = 1; T[i, 0, i] = 1; T[i, i, 0] = -1
T[1:, 1:, 1:] += f[1:, 1:, 1:]
mul = lambda x, y: np.einsum('ijk,i,j->k', T, x, y)
nrm = np.linalg.norm
phi = lambda a, b, c: np.einsum('ijk,i,j,k->', f[1:, 1:, 1:], a, b, c)   # 三つの虚部の 3 形式
IDX = list(itertools.product(range(3), repeat=3))


def arrL(S): return np.array([mul(mul(S[i], S[j]), S[k]) for i, j, k in IDX])      # (O_i O_j) O_k
def arrR(S): return np.array([mul(S[i], mul(S[j], S[k])) for i, j, k in IDX])      # O_i (O_j O_k)
def assoc(S): return mul(mul(S[0], S[1]), S[2]) - mul(S[0], mul(S[1], S[2]))
def wedge2(S): B = S[:, 1:]; return np.linalg.det(B @ B.T)                           # |b1∧b2∧b3|^2


def num_rank(fn, S, h=1e-6):
    x0 = S.ravel(); cols = []
    for i in range(x0.size):
        d = np.zeros(x0.size); d[i] = h
        cols.append((fn((x0 + d).reshape(3, 8)) - fn((x0 - d).reshape(3, 8))).ravel() / (2 * h))
    sv = np.linalg.svd(np.array(cols).T, compute_uv=False)
    return int((sv > 1e-6 * sv[0]).sum())


S = rng.normal(size=(3, 8)); XL, XR = arrL(S), arrR(S)

print("=== (1) 27 個の成分の種類と、括弧 ===")
kinds = {"三つとも同じ (iii)": [t for t in IDX if len(set(t)) == 1], "二つが同じ": [t for t in IDX if len(set(t)) == 2], "三つとも違う": [t for t in IDX if len(set(t)) == 3]}
for name, ts in kinds.items():
    d = max(nrm(XL[IDX.index(t)] - XR[IDX.index(t)]) for t in ts)
    print("  %-16s %2d 個   括弧による差の最大 %.1e" % (name, len(ts), d))
A = assoc(S)
print("  三つとも違う 6 個の括弧の差は ±[O1,O2,O3]: |[O1,O2,O3]| = %.4f、実部 %.1e" % (nrm(A), A[0]))

print("\n=== (2) 大きさ ===")
print("  |X_ijk| と |O_i||O_j||O_k| の差の最大: %.1e" % max(abs(nrm(XL[n]) - nrm(S[i]) * nrm(S[j]) * nrm(S[k])) for n, (i, j, k) in enumerate(IDX)))

print("\n=== (3) 実部 ===")
re = {t: XL[n][0] for n, t in enumerate(IDX)}
print("  添字を巡回させても同じ: 最大の差 %.1e" % max(abs(re[(i, j, k)] - re[(j, k, i)]) for i, j, k in IDX))
print("  二つの添字を入れ替えたときの差 (三つとも違う場合): Re X_123 - Re X_132 = %.6f、 -2φ(b1,b2,b3) = %.6f"
      % (re[(0, 1, 2)] - re[(0, 2, 1)], -2 * phi(S[0, 1:], S[1, 1:], S[2, 1:])))
print("  二つが同じ場合は入れ替えても同じ: Re X_112 - Re X_121 = %.1e" % (re[(0, 0, 1)] - re[(0, 1, 0)]))
print("  27 個の実部のうち、違う値の個数: %d" % len(set(np.round(list(re.values()), 9))))

print("\n=== (4) 独立な量の数 (数値で階数を数える) ===")
gram9 = lambda Q: np.r_[(Q @ Q.T)[np.triu_indices(3)], Q[:, 0]]
print("  状態の実数の数: 24")
print("  配列の全成分 (27×8) が持つ独立な数: %d" % num_rank(arrL, S))
print("  配列の実部 (27 個) が持つ独立な数: %d" % num_rank(lambda Q: arrL(Q)[:, 0], S))
print("  二体の量 (内積 6 個と、実部 a_i 3 個) が持つ独立な数: %d" % num_rank(gram9, S))
print("  二体の量に φ(b1,b2,b3) を足したもの: %d" % num_rank(lambda Q: np.r_[gram9(Q), phi(Q[0, 1:], Q[1, 1:], Q[2, 1:])], S))
print("  配列の実部に、二体の量と φ を足しても増えない: %d" % num_rank(lambda Q: np.r_[arrL(Q)[:, 0], gram9(Q), phi(Q[0, 1:], Q[1, 1:], Q[2, 1:])], S))

print("\n=== (5) 結合子の大きさは、二体の量と φ で決まる ===")
p = phi(S[0, 1:], S[1, 1:], S[2, 1:])
print("  φ^2 + |[O1,O2,O3]|^2 / 4 = %.6f、 |b1∧b2∧b3|^2 = %.6f" % (p ** 2 + nrm(A) ** 2 / 4, wedge2(S)))

print("\n=== (6) 三つの型 ===")
n1 = rng.normal(size=7); n1 /= nrm(n1)
cplx = np.array([np.r_[x, y * n1] for x, y in rng.normal(size=(3, 2))])
quat = np.zeros((3, 8)); quat[:, :4] = rng.normal(size=(3, 4))                 # 1, e1, e2, e3 は四元数
for name, Q in (("同じ複素平面の 3 個", cplx), ("同じ四元数の中の 3 個", quat), ("一般の 3 個", S)):
    X = arrL(Q); comm = max(nrm(mul(Q[i], Q[j]) - mul(Q[j], Q[i])) for i in range(3) for j in range(3))
    pp = phi(Q[0, 1:], Q[1, 1:], Q[2, 1:]); w = max(wedge2(Q), 0.0)
    print("  %-14s 交換子の最大 %.3f  結合子 %.3f  φ^2/|b∧b∧b|^2 = %s  配列の違う成分の数 %d"
          % (name, comm, nrm(assoc(Q)), ("%.4f" % (pp ** 2 / w)) if w > 1e-12 else "(虚部が 2 次元以下)", len({tuple(np.round(x, 8)) for x in X})))

print("\n=== (7) 虚部の和が 0 の 3 個 ===")
Z = rng.normal(size=(3, 8)); Z[:, 1:] -= Z[:, 1:].mean(axis=0)
print("  虚部の張る次元: %d  結合子 %.1e  φ %.1e  (虚部が 2 次元に収まり、必ず一つの四元数の中に入る)"
      % (np.linalg.matrix_rank(Z[:, 1:], tol=1e-9), nrm(assoc(Z)), phi(Z[0, 1:], Z[1, 1:], Z[2, 1:])))

print("\n=== (8) 配列を足して 3 個に戻す三通り (割り算なし) ===")
P = S.sum(axis=0); X3 = XL.reshape(3, 3, 3, 8)
print("  Σ_jk X_ijk と (O_i P)P の差 %.1e | Σ_ik X_ijk と (P O_j)P の差 %.1e | Σ_ij X_ijk と (P P)O_k の差 %.1e"
      % (nrm(X3.sum(axis=(1, 2)) - np.array([mul(mul(o, P), P) for o in S])), nrm(X3.sum(axis=(0, 2)) - np.array([mul(mul(P, o), P) for o in S])), nrm(X3.sum(axis=(0, 1)) - np.array([mul(mul(P, P), o) for o in S]))))
print("  全部足す: Σ_ijk X_ijk と P^3 の差 %.1e" % nrm(XL.sum(axis=0) - mul(mul(P, P), P)))

print("\n=== (9) 配列から状態は決まるか ===")
def cube_roots(c):
    r = nrm(c); b = nrm(c[1:]); th = np.arctan2(b, c[0]); n = c[1:] / b
    return [r ** (1 / 3) * np.r_[np.cos((th + 2 * np.pi * m) / 3), np.sin((th + 2 * np.pi * m) / 3) * n] for m in range(3)]
for name, Q in (("一般の 3 個", S), ("同じ複素平面の 3 個", cplx)):
    X = arrL(Q); roots = [cube_roots(X[IDX.index((i, i, i))]) for i in range(3)]
    ok = [c for c in itertools.product(range(3), repeat=3) if nrm(arrL(np.array([roots[i][c[i]] for i in range(3)])) - X) < 1e-8 * nrm(X)]
    print("  %-14s 対角 X_iii = O_i^3 の 3 乗根の組 27 通りのうち、同じ配列を与えるもの: %d 通り" % (name, len(ok)))
