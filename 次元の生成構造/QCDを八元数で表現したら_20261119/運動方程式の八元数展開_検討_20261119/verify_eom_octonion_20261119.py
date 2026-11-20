#!/usr/bin/env python3
"""「運動方程式の八元数展開」の検算。本文の項目番号 [Vn] ごとに成否を表示する。
実行: python3 verify_eom_octonion_20261119.py   (Python 3, NumPy のみ。乱数の種は固定)
"""
import numpy as np

rng = np.random.default_rng(20261119)
out = []


def check(tag, text, ok):
    out.append(bool(ok))
    print(f"[{tag}] {text:<58}: {bool(ok)}")


# ---------- 八元数 (本文 §1) ----------
LINES = [(1, 2, 3), (1, 4, 5), (1, 7, 6), (2, 4, 6), (2, 5, 7), (3, 4, 7), (3, 6, 5)]
f = np.zeros((8, 8, 8))
for (i, j, k) in LINES:
    for (a, b, c), s in (((i, j, k), 1), ((j, k, i), 1), ((k, i, j), 1), ((j, i, k), -1), ((i, k, j), -1), ((k, j, i), -1)):
        f[a, b, c] = s
E = np.eye(8)


def mul(x, y):
    z = np.zeros(8)
    z[0] = x[0] * y[0] - x[1:] @ y[1:]
    z[1:] = x[0] * y[1:] + y[0] * x[1:] + np.einsum("ijk,i,j->k", f[1:, 1:, 1:], x[1:], y[1:])
    return z


def conj(x):
    return np.concatenate([[x[0]], -x[1:]])


U3 = [E[1], E[2], E[3]]


def to_c(v):
    return v[0] + 1j * v[7]


def to_c3(v):
    return np.array([v @ u + 1j * (v @ mul(E[7], u)) for u in U3])


def emb(p):
    return sum(p[k].real * U3[k] + p[k].imag * mul(E[7], U3[k]) for k in range(3))


def as_real(T):
    """C^3 上の複素 3x3 行列を、八元数 8 成分に働く実 8x8 行列にする (1 と e7 には 0)。"""
    D = np.zeros((8, 8))
    for i in range(1, 7):
        D[:, i] = emb(T @ to_c3(E[i]))
    return D


def to_cmat(D):
    return np.array([to_c3(D @ u) for u in U3]).T


def Dab(a, b):
    """D(a,b)x = [[a,b],x] - 3{(ab)x - a(bx)}"""
    def act(x):
        c = mul(a, b) - mul(b, a)
        return mul(c, x) - mul(x, c) - 3 * (mul(mul(a, b), x) - mul(a, mul(b, x)))
    return np.array([act(E[i]) for i in range(8)]).T


def is_der(D):
    return all(np.allclose(D @ mul(E[i], E[j]), mul(D @ E[i], E[j]) + mul(E[i], D @ E[j])) for i in range(8) for j in range(8))


def rc3():
    return rng.normal(size=3) + 1j * rng.normal(size=3)


def tl(M):
    return M - np.trace(M) / 3 * np.eye(3)


x, y = rng.normal(size=8), rng.normal(size=8)
check("V1", "積: |xy| = |x||y|", np.isclose(np.linalg.norm(mul(x, y)), np.linalg.norm(x) * np.linalg.norm(y)))
a = rng.normal() + 1j * rng.normal()
p, c = rc3(), rc3()
av = a.real * E[0] + a.imag * E[7]
check("V2", "複素数の i = e7 を左から掛ける: [a psi] = a psi", np.allclose(to_c3(mul(av, emb(p))), a * p))
check("V2", "右から掛けると反線形: [psi a] = conj(a) psi", np.allclose(to_c3(mul(emb(p), av)), np.conj(a) * p))
check("V3", "内積: psi^dag chi = [chi psi*]_C", np.isclose(to_c(mul(emb(c), conj(emb(p)))), np.vdot(p, c)))

# ---------- 生成子 (本文 §2) ----------
lam = np.zeros((8, 3, 3), dtype=complex)
lam[0][0, 1] = lam[0][1, 0] = 1
lam[1][0, 1], lam[1][1, 0] = -1j, 1j
lam[2][0, 0], lam[2][1, 1] = 1, -1
lam[3][0, 2] = lam[3][2, 0] = 1
lam[4][0, 2], lam[4][2, 0] = -1j, 1j
lam[5][1, 2] = lam[5][2, 1] = 1
lam[6][1, 2], lam[6][2, 1] = -1j, 1j
lam[7] = np.diag([1, 1, -2]) / np.sqrt(3)
T = lam / 2
G = [as_real(-1j * t) for t in T]            # G_a = -i lambda_a / 2 を 8x8 実行列にしたもの
check("V4", "G_a は積を保つ無限小変換で、1 と e7 を動かさない",
      all(is_der(g) and np.allclose(g @ E[0], 0) and np.allclose(g @ E[7], 0) for g in G))
check("V4", "G_a は反対称な実行列", all(np.allclose(g, -g.T) for g in G))

# G_a を D(e_i,e_j) で書く (本文 表 2.2)
s3 = np.sqrt(3)
D_ = lambda i, j: Dab(E[i], E[j])
G_from_D = [(-D_(1, 5) + D_(2, 6)) / 12, (D_(1, 2) + D_(5, 6)) / 12, (D_(1, 6) + D_(2, 5)) / 12,
            (-D_(1, 4) + D_(3, 6)) / 12, (D_(1, 3) + D_(4, 6)) / 12, (-D_(2, 4) - D_(3, 5)) / 12,
            (D_(2, 3) - D_(4, 5)) / 12, (D_(1, 6) - D_(2, 5) + 2 * D_(3, 4)) / (12 * s3)]
check("V5", "表 2.2: G_a = (1/12) sum ±D(e_i,e_j)", all(np.allclose(g, h) for g, h in zip(G, G_from_D)))

# 構造定数
fabc = np.zeros((8, 8, 8))
for a_ in range(8):
    for b_ in range(8):
        cm = T[a_] @ T[b_] - T[b_] @ T[a_]
        for c_ in range(8):
            fabc[a_, b_, c_] = (np.trace(cm @ T[c_]) * 2 / 1j).real
check("V6", "[G_a, G_b] = f_abc G_c",
      all(np.allclose(G[a_] @ G[b_] - G[b_] @ G[a_], sum(fabc[a_, b_, c_] * G[c_] for c_ in range(8))) for a_ in range(8) for b_ in range(8)))

# ---------- A q (本文 §4.2) ----------
Aa = rng.normal(size=8)
A8 = sum(Aa[n] * G[n] for n in range(8))
A3 = sum(Aa[n] * (-1j * T[n]) for n in range(8))
q1 = rc3()
check("V7", "A q: 実 8x8 行列での作用 = 複素 3x3 行列での作用", np.allclose(to_c3(A8 @ emb(q1)), A3 @ q1))


def D_act(a, b, xx):
    cc = mul(a, b) - mul(b, a)
    return mul(cc, xx) - mul(xx, cc) - 3 * (mul(mul(a, b), xx) - mul(a, mul(b, xx)))


def G_act(n, xx):   # 表 2.2 を八元数の積だけで実行
    d = lambda i, j: D_act(E[i], E[j], xx)
    return [(-d(1, 5) + d(2, 6)) / 12, (d(1, 2) + d(5, 6)) / 12, (d(1, 6) + d(2, 5)) / 12, (-d(1, 4) + d(3, 6)) / 12,
            (d(1, 3) + d(4, 6)) / 12, (-d(2, 4) - d(3, 5)) / 12, (d(2, 3) - d(4, 5)) / 12,
            (d(1, 6) - d(2, 5) + 2 * d(3, 4)) / (12 * s3)][n]


check("V7", "A q: 八元数の積 (交換子と結合子) だけで計算しても同じ",
      np.allclose(sum(Aa[n] * G_act(n, emb(q1)) for n in range(8)), A8 @ emb(q1)))

# ---------- 交換子 (本文 §5.2) ----------
aa, bb, cc_, dd = (rng.normal(size=8) for _ in range(4))
Dm, Dn = Dab(aa, bb), Dab(cc_, dd)
check("V8", "[D(a,b), D(c,d)] = D(D(a,b)c, d) + D(c, D(a,b)d)",
      np.allclose(Dm @ Dn - Dn @ Dm, Dab(Dm @ cc_, dd) + Dab(cc_, Dm @ dd)))
Fa = rng.normal(size=8)
F8 = sum(Fa[n] * G[n] for n in range(8))
Cm = A8 @ F8 - F8 @ A8
check("V8", "[A, F] の成分 = f_abc A_b F_c", np.allclose(Cm, sum(fabc[b_, c_, a_] * Aa[b_] * Fa[c_] * G[a_] for a_ in range(8) for b_ in range(8) for c_ in range(8))))
check("V8", "[A, F] も 7 本の条件を満たす", np.allclose(np.einsum("ijk,jk->i", f[1:, 1:, 1:], Cm[1:, 1:]), 0) and np.allclose(Cm @ E[7], 0))

# ---------- 二つの値から回転 (本文 §5.3) ----------
psi, chi = rc3(), rc3()
Dpc = Dab(emb(psi), emb(chi))


def proj_su3(D):
    """e7 を動かさない部分への射影: P D = D + (1/4) sum_k (D e7)_k D(e_k, e7)"""
    v = D @ E[7]
    return D + 0.25 * sum(v[k] * Dab(E[k], E[7]) for k in range(1, 7))


P = proj_su3(Dpc)
check("V9", "P D(psi,chi) は e7 を動かさず、積を保つ", np.allclose(P @ E[7], 0) and is_der(P))
check("V9", "P D(psi,chi) = -3 tl(psi chi^dag - chi psi^dag)",
      np.allclose(to_cmat(P), -3 * tl(np.outer(psi, chi.conj()) - np.outer(chi, psi.conj()))))
check("V9", "D(psi,chi) e7 = 2i (conj psi x conj chi)", np.allclose(to_c3(Dpc @ E[7]), 2j * np.cross(psi.conj(), chi.conj())))
D1 = Dab(emb(psi), mul(E[7], emb(psi)))
check("V10", "D(psi, e7 psi) は e7 を動かさない", np.allclose(D1 @ E[7], 0))
check("V10", "D(psi, e7 psi) = 6i tl(psi psi^dag)", np.allclose(to_cmat(D1), 6j * tl(np.outer(psi, psi.conj()))))

# ---------- スピノルと色の流れ (本文 §4.1, §5.4) ----------
s0 = np.eye(2)
sig = [np.array([[0, 1], [1, 0]], dtype=complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], dtype=complex)]
g0 = np.block([[s0, 0 * s0], [0 * s0, -s0]]).astype(complex)
gam = [g0] + [np.block([[0 * s0, s], [-s, 0 * s0]]) for s in sig]
eta = np.diag([1, -1, -1, -1.0])
check("V11", "{gamma^mu, gamma^nu} = 2 eta^{mu nu}",
      all(np.allclose(gam[m] @ gam[n] + gam[n] @ gam[m], 2 * eta[m, n] * np.eye(4)) for m in range(4) for n in range(4)))
q = rng.normal(size=(4, 3)) + 1j * rng.normal(size=(4, 3))      # q[alpha, color]
for nu in range(4):
    w = g0 @ gam[nu]                                              # (gamma^0 gamma^nu)_{alpha beta}
    ja = np.array([np.einsum("aj,ab,ji,bi->", q.conj(), w, T[n], q) for n in range(8)])   # qbar gamma^nu T^a q
    N = np.einsum("ab,bi,aj->ij", w, q, q.conj())                 # N_ij = sum w_ab q_{b i} conj(q_{a j})
    J = -0.5j * tl(N)
    ok1 = np.allclose(ja.imag, 0) and np.allclose(sum(-1j * ja[n] * T[n] for n in range(8)), J)
    # 結合項: qbar i gamma^nu A q = A^a j^a
    lhs = np.einsum("aj,ab,ji,bi->", q.conj(), 1j * w, A3, q)
    ok2 = np.isclose(lhs, Aa @ ja.real)
    # 八元数の積だけで J を作る
    Jo = np.zeros((8, 8))
    for al in range(4):
        Jo += -(w[al, al].real / 12) * Dab(emb(q[al]), mul(E[7], emb(q[al])))
        for be in range(al + 1, 4):
            u, v = w[al, be].real, w[al, be].imag
            Jo += -(u / 6) * proj_su3(Dab(emb(q[be]), mul(E[7], emb(q[al])))) - (v / 6) * proj_su3(Dab(emb(q[be]), emb(q[al])))
    ok3 = np.allclose(Jo, as_real(J))
    check("V12", f"nu={nu}: J = -i j^a T^a = -(i/2) tl N、結合項、八元数の積での式", ok1 and ok2 and ok3)

# ---------- 時間方向の生成子 (本文 §4.3) ----------
alpha = [g0 @ gam[k] for k in (1, 2, 3)]
check("V13", "alpha^k = gamma^0 gamma^k と beta = gamma^0 はエルミート", all(np.allclose(m, m.conj().T) for m in alpha + [g0]))
Gen = -np.kron(np.eye(4), A3) - 1j * 0.7 * np.kron(g0, np.eye(3)) - sum(np.kron(alpha[k], sum(rng.normal() * (-1j * T[n]) for n in range(8))) for k in range(3))
check("V13", "-A_0 - i m beta - alpha^k A_k は反エルミート (12x12)", np.allclose(Gen, -Gen.conj().T))


def expm(M, terms=60):
    r, t = np.eye(len(M), dtype=complex), np.eye(len(M), dtype=complex)
    for n in range(1, terms):
        t = t @ M / n
        r = r + t
    return r


v0 = q.reshape(12)
dt = 0.1
check("V13", "exp(dt G) はノルムを保つ", np.isclose(np.linalg.norm(expm(dt * Gen) @ v0), np.linalg.norm(v0)))
check("V13", "一次の足し込み (1 + dt G) はノルムを保たない", not np.isclose(np.linalg.norm((np.eye(12) + dt * Gen) @ v0), np.linalg.norm(v0)))

# ---------- 格子 (本文 §6) ----------
def U_of(coef):
    return expm(sum(coef[n] * G[n] for n in range(8)).astype(complex)).real


Ul = U_of(rng.normal(size=8))
check("V14", "リンク U = exp(A) は直交行列で、積を保ち、e7 を動かさない",
      np.allclose(Ul @ Ul.T, np.eye(8)) and np.allclose(Ul @ mul(x, y), mul(Ul @ x, Ul @ y)) and np.allclose(Ul @ E[7], E[7]))
gx, gy = U_of(rng.normal(size=8)), U_of(rng.normal(size=8))
qx, qy = emb(rc3()), emb(rc3())
hop = lambda a_, U_, b_: to_c(mul(U_ @ b_, conj(a_)))          # a^dag (U b) = [(U b) a*]_C
check("V14", "q(x)^dag U q(y) は点ごとの取り替えで不変", np.isclose(hop(qx, Ul, qy), hop(gx @ qx, gx @ Ul @ gy.T, gy @ qy)))
Us = [U_of(rng.normal(size=8)) for _ in range(4)]
UP = Us[0] @ Us[1] @ Us[2].T @ Us[3].T
gs = [U_of(rng.normal(size=8)) for _ in range(4)]
UPg = (gs[0] @ Us[0] @ gs[1].T) @ (gs[1] @ Us[1] @ gs[2].T) @ (gs[3] @ Us[2] @ gs[2].T).T @ (gs[0] @ Us[3] @ gs[3].T).T
check("V14", "Tr_O(1 - U_P) は点ごとの取り替えで不変", np.isclose(np.trace(np.eye(8) - UP), np.trace(np.eye(8) - UPg)))

print(f"\n合計 {len(out)} 項目、成立 {sum(out)}、不成立 {len(out) - sum(out)}")
