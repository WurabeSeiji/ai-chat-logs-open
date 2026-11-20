"""本質の検算: 八元数の積ひとつから、ヤン＝ミルズの形とディラックの形が出るか。2026-11-20
発言録 [10][12] でアイリスが述べ、引き継ぎ文書で「未検算」としていた二点を確かめる。
"""
import numpy as np

rng = np.random.default_rng(2026112002)
TRIPLES = [(1, 2, 3), (1, 4, 5), (1, 7, 6), (2, 4, 6), (2, 5, 7), (3, 4, 7), (3, 6, 5)]
f = np.zeros((8, 8, 8))
for (i, j, k) in TRIPLES:
    for (a, b, c) in [(i, j, k), (j, k, i), (k, i, j)]:
        f[a, b, c] = 1
        f[b, a, c] = -1
M = np.zeros((8, 8, 8))
for i in range(8):
    M[0, i, i] = 1
    M[i, 0, i] = 1
for i in range(1, 8):
    M[i, i, 0] = -1
M[1:, 1:, 1:] += f[1:, 1:, 1:]
E = np.eye(8)


def mul(x, y):
    return np.einsum('ijk,i,j->k', M, x, y)


def conj(x):
    y = -x.copy()
    y[0] = x[0]
    return y


def cross(b, c):
    return np.einsum('ijk,i,j->k', f[1:, 1:, 1:], b, c)


def oexp(v):  # 虚八元数 v (7 成分) の exp = cos|v| + (v/|v|) sin|v|
    t = np.linalg.norm(v)
    return np.concatenate([[np.cos(t)], np.sin(t) * v / t])


print("[E1] 三つの積の実部 = 1 − ½|閉じ残り|²  (ヤン＝ミルズの形)")
print("     三辺 o = exp(εa), exp(εb), exp(εc)。閉じ残り 𝔉 = ε(a+b+c) + ε²(a×b + a×c + b×c)   (½[a,b] = a×b)")
a, b, c = rng.normal(size=7), rng.normal(size=7), rng.normal(size=7)
print("     (i) 一般の三辺")
for eps in (0.1, 0.05, 0.025, 0.0125):
    S = mul(mul(oexp(eps * a), oexp(eps * b)), oexp(eps * c))[0]
    F = eps * (a + b + c) + eps ** 2 * (cross(a, b) + cross(a, c) + cross(b, c))
    print(f"       ε={eps:<7} 1−Re(ooo) = {1-S:.6e}   ½|𝔉|² = {0.5*F@F:.6e}   比 = {(1-S)/(0.5*F@F):.5f}")
print("     (ii) 一次で閉じた三辺 (c = −a−b)。このとき 𝔉 = ε²(a×b) で、交換子(クロス積)だけが残る")
c2 = -a - b
for eps in (0.1, 0.05, 0.025, 0.0125):
    S = mul(mul(oexp(eps * a), oexp(eps * b)), oexp(eps * c2))[0]
    F = eps ** 2 * cross(a, b)
    print(f"       ε={eps:<7} 1−Re(ooo) = {1-S:.6e}   ½|𝔉|² = (ε⁴/2)|a×b|² = {0.5*F@F:.6e}   比 = {(1-S)/(0.5*F@F):.5f}")
print("     (iii) 同じ複素平面の三辺 (a, b が平行) では a×b = 0 で、一次で閉じれば閉じ残りは 0")
b3 = 0.7 * a
S = mul(mul(oexp(0.1 * a), oexp(0.1 * b3)), oexp(0.1 * (-a - b3)))[0]
print(f"       1−Re(ooo) = {1-S:.3e}")

print("\n[E2] 左からの掛け算は、二乗するとラプラシアンになる (ディラックの形)")
err1 = err2 = 0.0
for _ in range(300):
    x = rng.normal(size=8)
    for m in range(8):
        for n in range(8):
            lhs = mul(E[m], mul(conj(E[n]), x)) + mul(E[n], mul(conj(E[m]), x))
            err1 = max(err1, np.abs(lhs - 2 * (m == n) * x).max())
    p = rng.normal(size=8)          # 𝔇 = Σ e_μ p_μ を「運動量 p」で表したもの
    y = mul(p, mul(conj(p), x))
    err2 = max(err2, np.abs(y - (p @ p) * x).max())
print(f"     e_μ(ē_ν x) + e_ν(ē_μ x) = 2δ_μν x        最大誤差 = {err1:.2e}")
print(f"     𝔇(𝔇̄ x) = |p|² x  (𝔇𝔇̄ = □)              最大誤差 = {err2:.2e}")

print("\n[E3] 同じ三つの積を、物質・つなぎ・物質 と読む (ホッピング項の形)")
print("     Re(o_A o_B o_C) = <conj(o_A), o_B o_C> : conj(o_A) と、o_B を左から掛けた o_C との内積")
err = 0.0
for _ in range(300):
    oA, oB, oC = rng.normal(size=8), rng.normal(size=8), rng.normal(size=8)
    LB = np.einsum('ijk,i->kj', M, oB)        # o_B を左から掛ける 8×8 行列
    err = max(err, abs(mul(mul(oA, oB), oC)[0] - conj(oA) @ (LB @ oC)))
print(f"     Re(o_A o_B o_C) = ψ̄ᵀ (L_oB) ψ、ψ̄ = conj(o_A)、ψ = o_C     最大誤差 = {err:.2e}")
u = rng.normal(size=8)
u /= np.linalg.norm(u)
LB = np.einsum('ijk,i->kj', M, u)
print(f"     |o_B| = 1 のとき L_oB は直交行列 (回転): L Lᵀ − 1 の最大 = {np.abs(LB@LB.T-np.eye(8)).max():.2e}")
