#!/usr/bin/env python3
"""index_v3.html と同じ式を独立に計算して確認値を出す（確認用）。

基準点からの角距離 γ。倍音はルジャンドル多項式 P_l(cos γ)、l = 1..N。
  等振幅      A(γ) = (1/N) Σ_{l=1}^{N} P_l(cos γ)
  フェイエール A(γ) = Σ_{l=1}^{N} w_l P_l(cos γ) / Σ w_l,  w_l = 1 − l/(N+1)
P_l は三項漸化式 (l+1) P_{l+1} = (2l+1) x P_l − l P_{l−1}（index_v3.html と同じ）。

1) 軸対称なので γ ∈ [0, 180°] を 0.01° 刻みで評価し、A(0)、A(90°)、A(180°)、max、min、
   最初の零点（符号が変わる γ、線形補間）、半値全幅（A = 0.5 を最初に切る γ の 2 倍）を出す。
2) θlo = 180°, θla = 0° の 3D 格子（index_v3.html と同じ nL = min(721, max(181, 8N+1))、
   緯度は同じ刻み）で cos γ の式を通して A を計算し、max / min を 1 次元の値と並べる。
"""
import math

def weights(N, fejer):
    w = [0.0] * (N + 1)
    for l in range(1, N + 1):
        w[l] = 1 - l / (N + 1) if fejer else 1.0
    s = sum(w)
    return [x / s for x in w]

def legendre_sum(x, N, w):
    p0, p1 = 1.0, x
    v = w[1] * p1
    for l in range(1, N):
        p2 = ((2 * l + 1) * x * p1 - l * p0) / (l + 1)
        v += w[l + 1] * p2
        p0, p1 = p1, p2
    return v

def first_cross(gs, vals, level):
    """vals が level を初めて下回る γ（線形補間）。無ければ None。"""
    for i in range(1, len(vals)):
        if vals[i - 1] > level >= vals[i]:
            g0, g1, v0, v1 = gs[i - 1], gs[i], vals[i - 1], vals[i]
            return g0 + (g1 - g0) * (v0 - level) / (v0 - v1)
    return None

if __name__ == "__main__":
    step = 0.01
    gs = [i * step for i in range(int(round(180 / step)) + 1)]
    print("1D (gamma 0..180 deg, step 0.01 deg)")
    print("mode, N, A(0), A(90), A(180), max A, min A, first zero [deg], FWHM [deg]")
    for mode in ("equal", "fejer"):
        for N in (1, 2, 3, 10, 100):
            w = weights(N, mode == "fejer")
            vals = [legendre_sum(math.cos(math.radians(g)), N, w) for g in gs]
            z = first_cross(gs, vals, 0.0)
            h = first_cross(gs, vals, 0.5)
            zs = f"{z:.2f}" if z is not None else "none"
            hs = f"{2 * h:.2f}" if h is not None else "none"
            print(f"{mode}, {N}, {vals[0]:.4f}, {vals[9000]:.4f}, {vals[-1]:.4f}, {max(vals):.4f}, {min(vals):.4f}, {zs}, {hs}")

    print()
    print("3D grid (theta_lo = 180, theta_la = 0; same grid as index_v3.html)")
    print("mode, N, nL x nB, max A (3D), min A (3D), max A (1D), min A (1D)")
    lo0, la0 = math.radians(180.0), math.radians(0.0)
    sla, cla = math.sin(la0), math.cos(la0)
    for mode in ("equal", "fejer"):
        for N in (1, 10, 100):
            w = weights(N, mode == "fejer")
            nL = min(721, max(181, 8 * N + 1)); nB = (nL - 1) // 2 + 1
            amax, amin = -math.inf, math.inf
            for j in range(nB):
                beta = math.radians(-90 + 180 * j / (nB - 1)); sb, cb = math.sin(beta), math.cos(beta)
                for i in range(nL):
                    lam = math.radians(360 * i / (nL - 1))
                    cg = sb * sla + cb * cla * math.cos(lam - lo0)
                    cg = max(-1.0, min(1.0, cg))
                    a = legendre_sum(cg, N, w)
                    if a > amax: amax = a
                    if a < amin: amin = a
            vals = [legendre_sum(math.cos(math.radians(g)), N, w) for g in gs]
            print(f"{mode}, {N}, {nL} x {nB}, {amax:.4f}, {amin:.4f}, {max(vals):.4f}, {min(vals):.4f}")
