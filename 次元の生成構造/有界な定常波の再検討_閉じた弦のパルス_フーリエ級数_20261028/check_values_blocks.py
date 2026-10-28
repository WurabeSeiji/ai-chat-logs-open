#!/usr/bin/env python3
"""index_blocks.html と同じ式を独立に計算して確認値を出す（確認用）。

鏡 1 枚ごとの半球 ±1 の矩形波をルジャンドル級数で作り、全部の鏡について掛ける。
  S_L(u) = Σ_{l 奇数 ≤ L} c_l P_l(u),  c_l = P_{l−1}(0) − P_{l+1}(0)
  F(x)   = Π_{i=1}^{M} S_L(n_i · x)
鏡の組は index_blocks.html と同じ（二面体 4m、正四面体 24、立方体 48、正二十面体 120）。

1) 係数 c_1 … c_9 を分数の既知値と照合
2) 鏡の組の閉包：n_i を鏡 j で反射した結果が組の中の ±n_k に一致する（鏡映群であることの検算）
3) 乱数点 200000 個で符号ベクトル (sgn n_i·x) の種類を数え、ブロック数 B と各ブロックの面積比 × B を出す。
   いちばん深い点（min_i |n_i·x| が最大の点）の深さを角度で出す
4) L = 9, 49, 99 について、index_blocks.html と同じ格子（θlo = 180°, θla = 0°）で max F・min F、
   いちばん深い点の F、1 因子の行き過ぎ max S_L を出す
"""
import math
import numpy as np
from fractions import Fraction

PHI = (1 + 5 ** 0.5) / 2

def unit(v):
    s = math.sqrt(sum(x * x for x in v)); return [x / s for x in v]

def mirrors(mode, m):
    n = []
    if mode == 'dihedral':
        n.append([0.0, 0.0, 1.0])
        for k in range(m):
            n.append([math.cos(math.pi * k / m), math.sin(math.pi * k / m), 0.0])
    else:
        if mode != '24':
            n += [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        if mode != '120':
            n += [unit(v) for v in ([1, 1, 0], [1, -1, 0], [0, 1, 1], [0, 1, -1], [1, 0, 1], [1, 0, -1])]
        if mode == '120':
            for s1 in (1, -1):
                for s2 in (1, -1):
                    a = [s1 / PHI / 2, s2 * PHI / 2, 0.5]
                    n += [a, [a[2], a[0], a[1]], [a[1], a[2], a[0]]]
    return [[float(x) for x in v] for v in n]

def closure_failures(n):
    bad = 0
    for a in n:
        for b in n:
            d = sum(x * y for x, y in zip(a, b))
            r = [x - 2 * d * y for x, y in zip(a, b)]
            ok = any(max(abs(x - y) for x, y in zip(r, c)) < 1e-9 or max(abs(x + y) for x, y in zip(r, c)) < 1e-9 for c in n)
            if not ok:
                bad += 1
    return bad

def coeffs(L):
    Lo = L if L % 2 else L - 1
    P0 = [0.0] * (Lo + 2); P0[0] = 1.0; P0[1] = 0.0
    for k in range(1, Lo + 1):
        P0[k + 1] = -k * P0[k - 1] / (k + 1)
    c = [0.0] * (Lo + 1)
    for l in range(1, Lo + 1, 2):
        c[l] = P0[l - 1] - P0[l + 1]
    return c, Lo

def S(u, c, Lo):
    p0 = np.ones_like(u); p1 = u; v = c[1] * u
    for l in range(1, Lo):
        p2 = ((2 * l + 1) * u * p1 - l * p0) / (l + 1)
        p0, p1 = p1, p2
        if (l + 1) % 2:
            v = v + c[l + 1] * p2
    return v

def F(X, N, c, Lo):
    f = np.ones(X.shape[:-1])
    for n in N:
        u = np.clip(np.einsum('...j,j->...', X, np.array(n)), -1.0, 1.0)
        f = f * S(u, c, Lo)
    return f

def grid(Lo):
    nL = min(721, max(181, 8 * Lo + 1)); nB = (nL - 1) // 2 + 1
    lam = np.radians(np.linspace(0, 360, nL)); beta = np.radians(np.linspace(-90, 90, nB))
    B, Lm = np.meshgrid(beta, lam, indexing='ij')
    X = np.stack([np.cos(B) * np.cos(Lm), np.cos(B) * np.sin(Lm), np.sin(B)], axis=-1)
    return X, nL, nB

def rotator(lo_deg, la_deg):
    a = math.radians(90 - la_deg); b = math.radians(lo_deg)
    ca, sa, cb, sb = math.cos(a), math.sin(a), math.cos(b), math.sin(b)
    def rot(v):
        x1 = ca * v[0] + sa * v[2]; y1 = v[1]; z1 = -sa * v[0] + ca * v[2]
        return [cb * x1 - sb * y1, sb * x1 + cb * y1, z1]
    return rot

if __name__ == "__main__":
    c9, _ = coeffs(9)
    known = {1: Fraction(3, 2), 3: Fraction(-7, 8), 5: Fraction(11, 16), 7: Fraction(-75, 128), 9: Fraction(133, 256)}
    print("1) coefficients c_l = P_{l-1}(0) - P_{l+1}(0)")
    for l in (1, 3, 5, 7, 9):
        print(f"   c_{l} = {c9[l]:.10f}   known {known[l]} = {float(known[l]):.10f}   diff {abs(c9[l] - float(known[l])):.1e}")

    print()
    print("2)-4) mirror sets")
    rng = np.random.default_rng(0)
    P = rng.normal(size=(200000, 3)); P /= np.linalg.norm(P, axis=1, keepdims=True)
    rot = rotator(180.0, 0.0)
    Ls = (9, 49, 99)
    print("   1D Gibbs: max S_L(u), u in [-1,1] step 1e-5:", end="")
    uu = np.linspace(-1, 1, 200001)
    for L in Ls:
        c, Lo = coeffs(L); print(f"  L={Lo}: {S(uu, c, Lo).max():.4f}", end="")
    print()
    print("mode, m, M mirrors, B expected, closure failures, distinct sign vectors, min/max area*B, deepest point angle [deg]")
    print("    L, nL x nB, max F, min F, F(deepest)")
    for mode, m in (('dihedral', 2), ('dihedral', 6), ('24', 0), ('48', 0), ('120', 0)):
        N = mirrors(mode, m); M = len(N)
        B = 4 * m if mode == 'dihedral' else int(mode)
        bad = closure_failures(N)
        A = np.array(N)
        sg = np.sign(np.einsum('ij,kj->ik', P, A))
        keys, counts = np.unique(sg, axis=0, return_counts=True)
        frac = counts / len(P) * B
        depth = np.abs(np.einsum('ij,kj->ik', P, A)).min(axis=1)
        k = int(depth.argmax()); deep = P[k]; deep_deg = math.degrees(math.asin(depth[k]))
        print(f"{mode}, {m}, {M}, {B}, {bad}, {len(keys)}, {frac.min():.3f}/{frac.max():.3f}, {deep_deg:.2f}")
        Nr = [rot(v) for v in N]
        for L in Ls:
            c, Lo = coeffs(L)
            X, nL, nB = grid(Lo)
            f = F(X, Nr, c, Lo)
            fd = float(F(deep[None, :], N, c, Lo)[0])
            print(f"    {Lo}, {nL} x {nB}, {f.max():.4f}, {f.min():.4f}, {fd:.4f}")
