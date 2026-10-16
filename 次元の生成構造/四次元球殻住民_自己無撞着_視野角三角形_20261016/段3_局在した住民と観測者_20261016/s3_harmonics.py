# -*- coding: utf-8 -*-
"""S³ の倍音（調和関数）の係数表現（段 2 以降で共通に使うモジュール）
場 f(q) = Σ_n √(n+1) Σ_{k,l} C_n[k,l] D_n(V(q))[k,l]、q は（複素でもよい）四元数、V(q) は 2×2 行列。
D_n は SL(2,ℂ) の n 次の対称積表現（SU(2) で unitary、D(M₁M₂)=D(M₁)D(M₂)）。
左の並進 f(q)→f(a⁻¹q)：C_n → D_n(V(a⁻¹))ᵀ C_n。右の回転 f(q)→f(q b)：C_n → C_n D_n(V(b))ᵀ。
本当の二乗（S³ 上の f² の平均）：B(C) = Σ_n tr(C_nᵀ J_n⁻¹ C_n J_n)、J_n = D_n(ε)、ε=[[0,1],[-1,0]]（SL(2,ℂ) で不変）。
ノルム（f の絶対値二乗の平均）：H(C) = Σ_n Σ|C_n|²。"""
import numpy as np
from math import comb

def qmul(p, q):
    a1,b1,c1,d1 = p; a2,b2,c2,d2 = q
    return np.array([a1*a2-b1*b2-c1*c2-d1*d2, a1*b2+b1*a2+c1*d2-d1*c2, a1*c2-b1*d2+c1*a2+d1*b2, a1*d2+b1*c2-c1*b2+d1*a2])

def qexp(phi, u):
    """exp(φ u)、u は単位虚四元数の向き（'i','j','k' または 3 成分）。φ は複素でもよい（虚数の角＝ブースト）。"""
    vec = {'i': (1,0,0), 'j': (0,1,0), 'k': (0,0,1)}[u] if isinstance(u, str) else tuple(u)
    return np.array([np.cos(phi)] + [np.sin(phi)*c for c in vec], dtype=complex)

def qinv(q):
    n = q[0]**2 + q[1]**2 + q[2]**2 + q[3]**2      # 双線形のノルム（複素四元数でも可）
    return np.array([q[0], -q[1], -q[2], -q[3]])/n

def V(q):
    q0,q1,q2,q3 = q
    return np.array([[q0+1j*q3, q2+1j*q1], [-q2+1j*q1, q0-1j*q3]]).T

def Dn(M, n):
    a, b, c, d = M[0,0], M[0,1], M[1,0], M[1,1]
    A = np.zeros((n+1, n+1), complex)
    for l in range(n+1):
        for i in range(n-l+1):
            for j in range(l+1):
                A[i+j, l] += comb(n-l, i)*a**(n-l-i)*c**i*comb(l, j)*b**(l-j)*d**j
    s = np.array([np.sqrt(comb(n, k)) for k in range(n+1)])
    return (A*s[None, :])/s[:, None]

EPS = np.array([[0, 1], [-1, 0]], complex)
_J = {}
def J(n):
    if n not in _J: _J[n] = Dn(EPS, n)
    return _J[n]

def translate_left(C, a):
    ai = V(qinv(a)); return [Dn(ai, n).T @ C[n] for n in range(len(C))]

def rotate_right(C, b):
    vb = V(b); return [C[n] @ Dn(vb, n).T for n in range(len(C))]

def true_square(C):
    return sum(np.trace(C[n].T @ np.linalg.inv(J(n)) @ C[n] @ J(n)) for n in range(len(C)))

def norm_h(C):
    return sum(np.sum(np.abs(C[n])**2) for n in range(len(C))).real

def field(C, q):
    vq = V(q); return sum(np.sqrt(n+1)*np.sum(C[n]*Dn(vq, n)) for n in range(len(C)))

def icosians():
    """600 胞体の 120 頂点（単位イコシアン）"""
    import itertools
    phi = (1+5**0.5)/2; P = set()
    for i in range(4):
        for s in (1, -1):
            v = [0.0]*4; v[i] = s; P.add(tuple(v))
    for s in itertools.product((0.5, -0.5), repeat=4): P.add(s)
    base = (phi/2, 0.5, 1/(2*phi), 0.0)
    for perm in itertools.permutations(range(4)):
        if sum(1 for x in range(4) for y in range(x+1, 4) if perm[x] > perm[y]) % 2: continue
        t = tuple(base[perm[k]] for k in range(4))
        for sg in itertools.product((1, -1), repeat=4):
            P.add(tuple(x*s for x, s in zip(t, sg)))
    # 2026-10-16 修正：座標を 9 桁に丸めていたため |N(p)−1| が 8e-10 になり、段 2 の検査 5 が通らなかった。
    # 重複の判定だけを丸めたキーで行い、座標は丸めずに持つ。
    uniq = {}
    for v in P:
        uniq.setdefault(tuple(round(x, 9) for x in v), v)
    return np.array([uniq[k] for k in sorted(uniq)], dtype=float)
