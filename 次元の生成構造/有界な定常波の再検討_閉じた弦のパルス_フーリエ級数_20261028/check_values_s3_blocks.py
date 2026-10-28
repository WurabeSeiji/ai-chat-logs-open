#!/usr/bin/env python3
"""index_s3_blocks.html と同じ式を独立に計算して確認値を出す（確認用）。

殻：      f = sgn(cos Nz χ)、A(χ) = Σ_{N=1}^{L+1} b_N sin(Nχ)/sin χ、b_N = (2/π)∫ f sin χ sin(Nχ) dχ
ブロック：S_L(u) = Σ_{n 奇数 ≤ L} a_n U_n(u)、a_n = (−1)^{(n−1)/2}·4(n+1)/(π n(n+2))（S³ 上の半球 ±1 の倍音展開）
          F(x) = Π_i S_L(n_i·x)、x ∈ S³ ⊂ ℝ⁴

1) a_n の閉じた式を数値積分 (2/π)∫_0^π sgn(cos χ) sin((n+1)χ) sin χ dχ と突き合わせ。1 因子の行き過ぎ max S_L
2) b_N の閉じた式を数値積分と突き合わせ（Nz = 1, 2, 3）。直流 b_1、A = ±1/2 の交差半径（L = 49）
3) 4 次元の鏡の組の閉包（反射で組の ±法線に戻る）と、乱数点 600000 個の符号ベクトルの種類（＝ブロック数）、面積比 × B、いちばん深い点
4) 各モード、L = 9, 49, 99, 199：球体を覆う 61³ の格子（index_s3_blocks.html と同じ写像、σ = 0）で max F・min F、いちばん深い点の F
"""
import math
import numpy as np

def unit(v): v = np.asarray(v, float); return v / np.linalg.norm(v)

def mirrors4(mode, m=3, n=2):
    N = []
    if mode == 'torus':
        for k in range(m): N.append([math.cos(math.pi*k/m), math.sin(math.pi*k/m), 0, 0])
        for k in range(n): N.append([0, 0, math.cos(math.pi*k/n), math.sin(math.pi*k/n)])
    elif mode == '120':
        s5 = math.sqrt(5)
        V = [[1,1,1,-1/s5],[1,-1,-1,-1/s5],[-1,1,-1,-1/s5],[-1,-1,1,-1/s5],[0,0,0,4/s5]]
        for i in range(5):
            for j in range(i+1, 5): N.append(unit(np.array(V[i]) - np.array(V[j])))
    elif mode == '14400':
        import itertools
        PHI = (1 + 5 ** 0.5) / 2; roots = []
        for i in range(4):
            for s in (1, -1):
                v = [0.0]*4; v[i] = s; roots.append(v)
        for b in range(16): roots.append([-0.5 if (b >> k) & 1 else 0.5 for k in range(4)])
        base = [PHI/2, 0.5, 1/(2*PHI), 0.0]
        def parity(p):
            return sum(1 for i in range(4) for j in range(i+1, 4) if p[i] > p[j]) % 2
        for pm in itertools.permutations(range(4)):
            if parity(pm): continue
            for b in range(8):
                v = [0.0]*4; bit = 0
                for k in range(4):
                    val = base[pm[k]]
                    if val != 0: v[k] = (-1 if (b >> bit) & 1 else 1) * val; bit += 1
                roots.append(v)
        for v in roots:
            k = next(i for i, x in enumerate(v) if abs(x) > 1e-12)
            if v[k] > 0: N.append(v)
    else:
        for i in range(4):
            for j in range(i+1, 4):
                for s in (1, -1):
                    v = [0,0,0,0]; v[i] = 1; v[j] = s; N.append(unit(v))
        if mode in ('384', '1152'):
            for i in range(4):
                v = [0,0,0,0]; v[i] = 1; N.append(v)
        if mode == '1152':
            for s1 in (1,-1):
                for s2 in (1,-1):
                    for s3 in (1,-1): N.append([0.5, 0.5*s1, 0.5*s2, 0.5*s3])
    return np.array(N, float)

def closure_failures(A):
    bad = 0
    for a in A:
        for b in A:
            r = a - 2*np.dot(a, b)*b
            if not (np.min(np.max(np.abs(A - r), axis=1)) < 1e-9 or np.min(np.max(np.abs(A + r), axis=1)) < 1e-9): bad += 1
    return bad

def hemi_coeffs(L):
    Lo = L if L % 2 else L - 1
    a = np.zeros(Lo + 1)
    for n in range(1, Lo + 1, 2): a[n] = (-1)**((n-1)//2) * 4*(n+1) / (math.pi*n*(n+2))
    return a, Lo

def hemi_series(u, a, Lo):
    u = np.asarray(u, float); u0 = np.ones_like(u); u1 = 2*u; v = a[1]*u1
    for n in range(1, Lo):
        u2 = 2*u*u1 - u0; u0, u1 = u1, u2
        if (n+1) % 2: v = v + a[n+1]*u2
    return v

def shell_coeffs(Nz, L):
    Nmax = L + 1; b = np.zeros(Nmax + 1)
    def I(N, x): return x/2 - math.sin(2*x)/4 if N == 1 else math.sin((N-1)*x)/(2*(N-1)) - math.sin((N+1)*x)/(2*(N+1))
    for N in range(1, Nmax + 1):
        s = 0.0; x0 = 0.0
        for j in range(Nz + 1):
            x1 = (j + 0.5)*math.pi/Nz if j < Nz else math.pi
            s += (-1)**j * (I(N, x1) - I(N, x0)); x0 = x1
        b[N] = 2/math.pi * s
    return b, Nmax

def shell_series(x, b, Nmax):
    x = np.asarray(x, float); u0 = np.ones_like(x); u1 = 2*x; v = b[1]*u0
    if Nmax >= 2: v = v + b[2]*u1
    for N in range(3, Nmax + 1):
        u2 = 2*x*u1 - u0; u0, u1 = u1, u2; v = v + b[N]*u2
    return v

def crossings(vals, gs, level):
    out = []
    for i in range(1, len(vals)):
        if (vals[i-1]-level)*(vals[i]-level) < 0:
            out.append(gs[i-1] + (gs[i]-gs[i-1])*(vals[i-1]-level)/(vals[i-1]-vals[i]))
    return out

if __name__ == "__main__":
    chi = np.linspace(0, math.pi, 2000001); w = chi[1]-chi[0]
    print("1) hemisphere coefficients a_n (closed form vs quadrature)")
    f = np.sign(np.cos(chi)); f[len(f)//2] = 0
    for n in (1, 3, 5, 7):
        q = 2/math.pi * np.trapezoid(f*np.sin((n+1)*chi)*np.sin(chi), dx=w)
        c = (-1)**((n-1)//2) * 4*(n+1)/(math.pi*n*(n+2))
        print(f"   a_{n} = {c:.8f}   quadrature {q:.8f}   diff {abs(c-q):.1e}")
    uu = np.linspace(-1, 1, 200001)
    print("   1D Gibbs max S_L:", "  ".join(f"L={Lo}: {hemi_series(uu, a, Lo).max():.4f}" for a, Lo in (hemi_coeffs(L) for L in (9, 49, 99))))

    print()
    print("2) shell coefficients b_N (closed form vs quadrature), DC, crossings at L = 49")
    gs = np.linspace(0, 180, 18001)
    for Nz in (1, 2, 3):
        fz = np.sign(np.cos(Nz*chi)); fz[np.abs(np.cos(Nz*chi)) < 1e-12] = 0
        b, Nmax = shell_coeffs(Nz, 49)
        worst = 0.0
        for N in (1, 2, 3, 4, 5, 6):
            q = 2/math.pi * np.trapezoid(fz*np.sin(chi)*np.sin(N*chi), dx=w)
            worst = max(worst, abs(q - b[N]))
        vals = shell_series(np.cos(np.radians(gs)), b, Nmax)
        cp = crossings(vals, gs, 0.5); cm = crossings(vals, gs, -0.5)
        print(f"   Nz = {Nz}: max|b_N - quad| (N=1..6) = {worst:.1e}, DC b_1 = {b[1]:.4f}, max A = {vals.max():.4f}, min A = {vals.min():.4f}")
        print(f"      zone edges {[round((j+0.5)*180/Nz, 1) for j in range(Nz)]}")
        print(f"      A=+1/2 at {[round(c, 1) for c in cp]}   A=-1/2 at {[round(c, 1) for c in cm]}")

    print()
    print("3)-4) mirror sets in R^4")
    rng = np.random.default_rng(0)
    P = rng.normal(size=(600000, 4)); P /= np.linalg.norm(P, axis=1, keepdims=True)
    Gc = 61; Rb = math.pi
    g1 = np.linspace(-Rb, Rb, Gc); gx, gy, gz = np.meshgrid(g1, g1, g1, indexing='ij')
    r = np.sqrt(gx**2 + gy**2 + gz**2); inside = r <= math.pi
    s = np.where(r > 1e-12, np.sin(r)/np.maximum(r, 1e-12), 0.0)
    X4 = np.stack([np.cos(r), s*gx, s*gy, s*gz], axis=-1)[inside]
    print("mode, M mirrors, B expected, closure failures, distinct sign vectors, min/max area*B, deepest point angle [deg]")
    print("    L, grid points inside ball, max F, min F, F(deepest)")
    for mode, B in (('torus', 24), ('120', 120), ('192', 192), ('384', 384), ('1152', 1152), ('14400', 14400)):
        A = mirrors4(mode); M = len(A)
        bad = closure_failures(A)
        sgn = np.sign(np.einsum('ij,kj->ik', P, A)).astype(np.int8)
        keys, counts = np.unique(sgn, axis=0, return_counts=True)
        frac = counts/len(P)*B
        depth = np.abs(np.einsum('ij,kj->ik', P, A)).min(axis=1); k = int(depth.argmax()); deep = P[k]
        print(f"{mode}, {M}, {B}, {bad}, {len(keys)}, {frac.min():.3f}/{frac.max():.3f}, {math.degrees(math.asin(depth[k])):.2f}")
        for L in (9, 49, 99, 199):
            a, Lo = hemi_coeffs(L)
            F = np.ones(len(X4))
            for nv in A: F = F * hemi_series(np.clip(np.einsum('ij,j->i', X4, nv), -1, 1), a, Lo)
            Fd = float(np.prod([hemi_series(np.clip(np.dot(deep, nv), -1, 1), a, Lo) for nv in A]))
            print(f"    {Lo}, {len(X4)}, {F.max():.4f}, {F.min():.4f}, {Fd:.4f}")
