#!/usr/bin/env python3
"""index_s3.html と同じ式を独立に計算して確認値を出す（確認用）。

S³（半径 R）上の中心位相からの角距離 χ。倍音は Z_n(χ) = sin((n+1)χ) / ((n+1) sin χ) = U_n(cos χ)/(n+1)。
  等振幅      A(χ) = (1/N) Σ_{n=1}^{N} Z_n(χ)
  フェイエール A(χ) = Σ w_n Z_n(χ) / Σ w_n,  w_n = 1 − n/(N+1)
index_s3.html はチェビシェフ第 2 種の三項漸化式で計算する。ここでは sin の式で直接計算し、漸化式とも突き合わせる。

1) χ ∈ [0, 180°] を 0.01° 刻みで評価し、A(0)、A(90°)、A(180°)、max、min、最初の零点、半値全幅、
   S³ 上の平均（重み sin²χ、台形則）を出す。漸化式との最大差も出す
2) 投影の検算：観測者中心の主観座標系の式 cos χ = cos ρ cos σ + sin ρ sin σ (u·ŝ) を、
   4 次元の埋め込み X = (cos ρ, sin ρ u)、S = (cos σ, sin σ ŝ) の内積と乱数点 100000 組で突き合わせる
3) index_s3.html と同じ円盤格子（θlo = 180°, θla = 0°, σ = 60°）で max A・min A を出し、1 次元の値と並べる
"""
import math
import random

def weights(N, fejer):
    w = [0.0] * (N + 1)
    for n in range(1, N + 1):
        w[n] = 1 - n / (N + 1) if fejer else 1.0
    s = sum(w)
    return [x / s for x in w]

def Z_direct(n, chi):
    s = math.sin(chi)
    if abs(s) < 1e-12:
        return 1.0 if math.cos(chi) > 0 else (-1.0) ** n
    return math.sin((n + 1) * chi) / ((n + 1) * s)

def A_direct(chi, N, w):
    return sum(w[n] * Z_direct(n, chi) for n in range(1, N + 1))

def A_recur(x, N, w):
    u0, u1 = 1.0, 2.0 * x
    v = w[1] * u1 / 2
    for n in range(1, N):
        u2 = 2.0 * x * u1 - u0
        u0, u1 = u1, u2
        v += w[n + 1] * u2 / (n + 2)
    return v

def first_cross(gs, vals, level):
    for i in range(1, len(vals)):
        if vals[i - 1] > level >= vals[i]:
            g0, g1, v0, v1 = gs[i - 1], gs[i], vals[i - 1], vals[i]
            return g0 + (g1 - g0) * (v0 - level) / (v0 - v1)
    return None

if __name__ == "__main__":
    step = 0.01
    gs = [i * step for i in range(int(round(180 / step)) + 1)]
    print("1) 1D (chi 0..180 deg, step 0.01 deg)")
    print("mode, N, A(0), A(90), A(180), max A, min A, first zero [deg], FWHM [deg], mean on S3 (sin^2 weight), max|direct-recur|")
    for mode in ("equal", "fejer"):
        for N in (1, 2, 3, 10, 100):
            w = weights(N, mode == "fejer")
            vals = [A_direct(math.radians(g), N, w) for g in gs]
            diff = max(abs(v - A_recur(math.cos(math.radians(g)), N, w)) for g, v in zip(gs[::7], vals[::7]))
            z = first_cross(gs, vals, 0.0); h = first_cross(gs, vals, 0.5)
            zs = f"{z:.2f}" if z is not None else "none"; hs = f"{2 * h:.2f}" if h is not None else "none"
            num = sum((vals[i] * math.sin(math.radians(gs[i])) ** 2 + vals[i + 1] * math.sin(math.radians(gs[i + 1])) ** 2) / 2 for i in range(len(gs) - 1))
            den = sum((math.sin(math.radians(gs[i])) ** 2 + math.sin(math.radians(gs[i + 1])) ** 2) / 2 for i in range(len(gs) - 1))
            print(f"{mode}, {N}, {vals[0]:.4f}, {vals[9000]:.4f}, {vals[-1]:.4f}, {max(vals):.4f}, {min(vals):.4f}, {zs}, {hs}, {num / den:.2e}, {diff:.1e}")

    print()
    print("2) projection check: chart formula vs 4D dot product, 100000 random samples")
    rnd = random.Random(0)
    worst = 0.0
    for _ in range(100000):
        rho = rnd.uniform(0, math.pi); sig = rnd.uniform(0, math.pi)
        u = [rnd.gauss(0, 1) for _ in range(3)]; nu = math.sqrt(sum(c * c for c in u)); u = [c / nu for c in u]
        lo = rnd.uniform(0, 2 * math.pi); la = rnd.uniform(-math.pi / 2, math.pi / 2)
        sh = [math.cos(la) * math.cos(lo), math.cos(la) * math.sin(lo), math.sin(la)]
        c1 = math.cos(rho) * math.cos(sig) + math.sin(rho) * math.sin(sig) * sum(a * b for a, b in zip(u, sh))
        X = [math.cos(rho)] + [math.sin(rho) * c for c in u]
        S = [math.cos(sig)] + [math.sin(sig) * c for c in sh]
        c2 = sum(a * b for a, b in zip(X, S))
        worst = max(worst, abs(c1 - c2))
    print(f"   max |cos chi (chart) - cos chi (4D dot)| = {worst:.1e}")

    print()
    print("3) disc grids as in index_s3.html (theta_lo = 180, theta_la = 0, sigma = 60)")
    print("mode, N, nR x nA, meridian disc max/min, horizontal disc max/min, 1D max/min")
    lo, la, sig = math.radians(180.0), math.radians(0.0), math.radians(60.0)
    h = [math.cos(lo), math.sin(lo), 0.0]
    sh = [math.cos(la) * h[0], math.cos(la) * h[1], math.sin(la)]
    planes = {"meridian": (h, [0.0, 0.0, 1.0]), "horizontal": ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0])}
    for mode in ("equal", "fejer"):
        for N in (10, 100):
            w = weights(N, mode == "fejer")
            nR = min(405, max(181, 4 * (N + 1) + 1)); nA = min(721, max(361, 8 * (N + 1) + 1))
            res = {}
            for name, (e1, e2) in planes.items():
                amax, amin = -math.inf, math.inf
                for j in range(nR):
                    rho = math.pi * j / (nR - 1); cr, sr = math.cos(rho), math.sin(rho)
                    for i in range(nA):
                        al = 2 * math.pi * i / (nA - 1); ca, sa = math.cos(al), math.sin(al)
                        u = [ca * e1[k] + sa * e2[k] for k in range(3)]
                        cx = cr * math.cos(sig) + sr * math.sin(sig) * sum(a * b for a, b in zip(u, sh))
                        cx = max(-1.0, min(1.0, cx))
                        a = A_recur(cx, N, w)
                        if a > amax: amax = a
                        if a < amin: amin = a
                res[name] = (amax, amin)
            vals = [A_direct(math.radians(g), N, w) for g in gs]
            print(f"{mode}, {N}, {nR} x {nA}, {res['meridian'][0]:.4f}/{res['meridian'][1]:.4f}, {res['horizontal'][0]:.4f}/{res['horizontal'][1]:.4f}, {max(vals):.4f}/{min(vals):.4f}")
