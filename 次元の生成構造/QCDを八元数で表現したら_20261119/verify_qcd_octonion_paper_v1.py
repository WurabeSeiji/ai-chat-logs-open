#!/usr/bin/env python3
"""「QCD を八元数で表現したら」の検算（v1.0）。

本文の式番号ごとに、乱数で作った値を使って成否を表示する。
実行: python3 verify_qcd_octonion_paper_v1.py   (Python 3 と NumPy のみ)
乱数の種は固定してある。すべての行が True なら本文の式が成り立っている。
"""
import numpy as np

rng = np.random.default_rng(20261119)

# ---- (3.1) 積の表: e_i e_j = -δ_ij + f_ijk e_k ---------------------------------
LINES = [(1, 2, 3), (1, 4, 5), (1, 7, 6), (2, 4, 6), (2, 5, 7), (3, 4, 7), (3, 6, 5)]
f = np.zeros((8, 8, 8))
for (i, j, k) in LINES:
    for (a, b, c), s in (((i, j, k), 1), ((j, k, i), 1), ((k, i, j), 1),
                         ((j, i, k), -1), ((i, k, j), -1), ((k, j, i), -1)):
        f[a, b, c] = s
E = np.eye(8)


def mul(x, y):
    z = np.zeros(8)
    z[0] = x[0] * y[0] - x[1:] @ y[1:]
    z[1:] = x[0] * y[1:] + y[0] * x[1:] + np.einsum("ijk,i,j->k", f[1:, 1:, 1:], x[1:], y[1:])
    return z


def conj(x):
    return np.concatenate([[x[0]], -x[1:]])


x, y, z = (rng.normal(size=8) for _ in range(3))
print("(3.1) ノルムの乗法性 |xy| = |x||y|        :", np.isclose(np.linalg.norm(mul(x, y)), np.linalg.norm(x) * np.linalg.norm(y)))
print("(3.1) 結合則は成り立たない                :", not np.allclose(mul(mul(x, y), z), mul(x, mul(y, z))))

# ---- (3.2) 分解 O = C ⊕ C^3 (i は e7 の左掛け、基底 e1,e2,e3) --------------------
U = [E[1], E[2], E[3]]


def to_c(v):
    return v[0] + 1j * v[7]


def to_c3(v):
    return np.array([v @ u + 1j * (v @ mul(E[7], u)) for u in U])


def emb(a, psi):
    return a.real * E[0] + a.imag * E[7] + sum(psi[k].real * U[k] + psi[k].imag * mul(E[7], U[k]) for k in range(3))


print("(3.2) e7e1 = -e6, e7e2 = e5, e7e3 = e4    :",
      np.allclose(mul(E[7], E[1]), -E[6]) and np.allclose(mul(E[7], E[2]), E[5]) and np.allclose(mul(E[7], E[3]), E[4]))

# ---- (3.3)(3.4) 積の成分表示 ----------------------------------------------------
a, b = rng.normal() + 1j * rng.normal(), rng.normal() + 1j * rng.normal()
psi, chi = rng.normal(size=3) + 1j * rng.normal(size=3), rng.normal(size=3) + 1j * rng.normal(size=3)
p = mul(emb(a, psi), emb(b, chi))
print("(3.3) C 成分  = ab − χ†ψ                  :", np.isclose(to_c(p), a * b - np.vdot(chi, psi)))
print("(3.3) C^3 成分 = aχ + b̄ψ + ψ̄×χ̄            :",
      np.allclose(to_c3(p), a * chi + np.conj(b) * psi + np.cross(np.conj(psi), np.conj(chi))))
print("(3.4) ψ†χ = [χ ψ*]_C                      :", np.isclose(to_c(mul(emb(0j, chi), conj(emb(0j, psi)))), np.vdot(psi, chi)))

av = a.real * E[0] + a.imag * E[7]
print("(3.2) 係数は左から: [aψ]_C^3 = aψ          :", np.allclose(to_c3(mul(av, emb(0j, psi))), a * psi))
print("§9   右から掛けると反線形: ψa = āψ        :", np.allclose(to_c3(mul(emb(0j, psi), av)), np.conj(a) * psi))

# ---- (4.1) 積を保つ無限小変換は 14 次元、e7 を動かさないものは 8 次元 ------------
rows = []
for i in range(8):
    for j in range(8):
        for k in range(8):
            r = np.zeros((8, 8))
            r[k, :] += np.array([mul(E[i], E[j]) @ E[m] for m in range(8)])   # D(e_i e_j)
            r[:, i] -= np.array([mul(E[m], E[j])[k] for m in range(8)])       # (D e_i) e_j
            r[:, j] -= np.array([mul(E[i], E[m])[k] for m in range(8)])       # e_i (D e_j)
            rows.append(r.ravel())
_, s, vt = np.linalg.svd(np.array(rows))
der = vt[np.sum(s > 1e-9):].reshape(-1, 8, 8)
print("(4.1) 積を保つ無限小変換は 14 次元        :", len(der) == 14)
img = np.array([d[:, 7] for d in der]).T
_, s2, vt2 = np.linalg.svd(img)
stab = [np.tensordot(c, der, axes=1) for c in vt2[np.sum(s2 > 1e-9):]]
print("(4.1) e7 を動かさないものは 8 次元        :", len(stab) == 8)
Ts = [np.array([to_c3(D @ u) for u in U]).T for D in stab]
print("(4.1) C^3 上で反エルミート・トレース 0    :",
      all(np.allclose(T, -T.conj().T) and np.isclose(np.trace(T), 0) for T in Ts))
print("(4.1) 1 と e7 には作用しない              :", all(np.allclose(D @ E[0], 0) and np.allclose(D @ E[7], 0) for D in stab))

# ---- (4.2) ゲルマン行列 −iλ^a/2 は積を保つ ---------------------------------------
lam = np.zeros((8, 3, 3), dtype=complex)
lam[0][0, 1] = lam[0][1, 0] = 1
lam[1][0, 1], lam[1][1, 0] = -1j, 1j
lam[2][0, 0], lam[2][1, 1] = 1, -1
lam[3][0, 2] = lam[3][2, 0] = 1
lam[4][0, 2], lam[4][2, 0] = -1j, 1j
lam[5][1, 2] = lam[5][2, 1] = 1
lam[6][1, 2], lam[6][2, 1] = -1j, 1j
lam[7] = np.diag([1, 1, -2]) / np.sqrt(3)


def as_real(T):
    D = np.zeros((8, 8))
    for i in range(1, 7):
        D[:, i] = emb(0j, T @ to_c3(E[i]))
    return D


gens = [as_real(-0.5j * l) for l in lam]


def is_derivation(D):
    return all(np.allclose(D @ mul(E[i], E[j]), mul(D @ E[i], E[j]) + mul(E[i], D @ E[j]))
               for i in range(8) for j in range(8))


print("(4.2) −iλ^a/2 は八元数の積を保つ          :", all(is_derivation(D) for D in gens))

# ---- (5.1) 第 2 項: Σ_α [χ_α (q_α)*]_C = q†χ (χ は γ^0(iγ^μD_μ − m)q にあたる) -----
qs = rng.normal(size=(4, 3)) + 1j * rng.normal(size=(4, 3))
cs = rng.normal(size=(4, 3)) + 1j * rng.normal(size=(4, 3))
print("(5.1) Σ_α [χ_α (q_α)*]_C = q†χ            :",
      np.isclose(sum(to_c(mul(emb(0j, cs[k]), conj(emb(0j, qs[k])))) for k in range(4)), np.vdot(qs, cs)))

# ---- (5.2) トレースの関係 --------------------------------------------------------
Fa = rng.normal(size=8)
F = sum(Fa[n] * gens[n] for n in range(8))
print("(5.2) Tr_O(FF) = −F^a F^a                 :", np.isclose(np.trace(F @ F), -(Fa @ Fa)))

# ---- (6.1)(6.2) 第 1 項の成分展開と 7 本の条件 -----------------------------------
print("(6.1) Σ_ij F^ij F^ij = F^a F^a            :", np.isclose((F[1:7, 1:7] ** 2).sum(), Fa @ Fa))
print("(6.2) f_ijk F^jk = 0 (i = 1..7)           :", np.allclose(np.einsum("ijk,jk->i", f[1:, 1:, 1:], F[1:, 1:]), 0))
c = lambda i, j: F[i, j]
seven = [c(2, 3) + c(4, 5), c(1, 3) - c(4, 6), c(1, 2) - c(5, 6), c(1, 5) + c(2, 6),
         c(1, 4) + c(3, 6), c(2, 4) - c(3, 5), c(1, 6) - c(2, 5) - c(3, 4)]
print("(6.2) 7 本を書き下した形                  :", np.allclose(seven, 0))
idx = [(i, j) for i in range(1, 7) for j in range(i + 1, 7)]
A = np.array([[2 * f[i, p_, q_] for (p_, q_) in idx] for i in range(1, 8)])
print("(6.2) 解の空間は 15 − 7 = 8 次元          :", len(idx) - np.linalg.matrix_rank(A) == 8)

# ---- (7.1) リンク変数とウィルソン作用 --------------------------------------------
def expm(M, terms=40):
    out, term = np.eye(len(M), dtype=M.dtype), np.eye(len(M), dtype=M.dtype)
    for n in range(1, terms):
        term = term @ M / n
        out = out + term
    return out


T3 = sum(Fa[n] * (-0.5j * lam[n]) for n in range(8))
U8, U3 = expm(F), expm(T3)
print("(7.1) U は積を保ち e7 を動かさない        :",
      np.allclose(U8 @ mul(x, y), mul(U8 @ x, U8 @ y)) and np.allclose(U8 @ E[7], E[7]))
print("(4.1) U は C^3 上でユニタリ、行列式 1      :", np.allclose(U3.conj().T @ U3, np.eye(3)) and np.isclose(np.linalg.det(U3), 1))
print("(7.1) Tr_O U = 2 + 2 Re tr_3 U            :", np.isclose(np.trace(U8), 2 + 2 * np.trace(U3).real))
print("(7.1) 6(1 − Re tr_3 U/3) = Tr_O(1 − U)    :", np.isclose(6 * (1 - np.trace(U3).real / 3), np.trace(np.eye(8) - U8)))

# ---- (8.1) 著者の記法 o = q + pℓ との対応 (ℓ = e7 を右から掛ける) ------------------
t, R = rng.normal(), rng.normal()
xs, Qs = rng.normal(size=3), rng.normal(size=3)
q = t * E[0] + xs @ E[1:4]
pq = R * E[0] + Qs @ E[1:4]
o = q + mul(pq, E[7])
print("(8.1) C 成分 = t + iR                     :", np.isclose(to_c(o), t + 1j * R))
print("(8.1) C^3 成分 = x − iQ                   :", np.allclose(to_c3(o), xs - 1j * Qs))
