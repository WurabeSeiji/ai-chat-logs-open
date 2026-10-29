#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""index_v2c.html と同じ式を独立に計算する確認用（check_values_v2.py のパルス部分を複素化したもの、2026-10-29）。

複素パルス  Z(ψ) = φ/360 + (2/π) Σ_m [ sin(mπφ/360) / m ] · e^{imπ(ψ−θ)/180}、m は 1..M のうち奇・偶で選んだもの
  実部 = index_v2.html のパルスの高さ z と同じ式、虚部 = 同じ係数の sin の和
  R = √Σ c_m²（直流を含めない）、3D の点は ((R + Re Z) cos ψ, (R + Re Z) sin ψ, Im Z)
出す量：max/min Re（元の max/min z と一致するはず）、max/min Im、R、線形重心、パワー重心、
       c₀/R、P(π)/R²、Σ c_m² e^{−2imθ}/R²、⟨m⟩（平均の波数）、巻き数（|w| = 0.999 の円上）
ψ の刻みは index_v2c.html と同じ（最高波数の 1 周期あたり 8 点以上、最低 1441 点）。
"""
import math

def coeffs(phi, M, odd=True, even=True):
    ms = [m for m in range(1, M + 1) if (m % 2 == 1 and odd) or (m % 2 == 0 and even)]
    cs = [(2 / math.pi) * math.sin(m * math.pi * phi / 360) / m for m in ms]
    return ms, cs, phi / 360

def pulse(psi, theta, ms, cs, c0, r=1.0):
    x = math.pi * (psi - theta) / 180
    re = c0 + sum(c * r ** m * math.cos(m * x) for m, c in zip(ms, cs))
    im = sum(c * r ** m * math.sin(m * x) for m, c in zip(ms, cs))
    return re, im

def run(theta, phi, M, odd, even):
    ms, cs, c0 = coeffs(phi, M, odd, even)
    R = math.sqrt(sum(c * c for c in cs))
    kmax = ms[-1] if ms else 1
    n = min(40001, max(1441, 8 * kmax + 1))
    vals = [pulse(i * 360 / (n - 1), theta, ms, cs, c0) for i in range(n)]
    re = [v[0] for v in vals]; im = [v[1] for v in vals]
    lx = ly = px = py = pw = 0.0
    for i in range(n - 1):
        a = R + re[i]; q = re[i] ** 2 + im[i] ** 2; c = math.cos(2 * math.pi * i / (n - 1)); s = math.sin(2 * math.pi * i / (n - 1))
        lx += a * c; ly += a * s; px += q * R * c; py += q * R * s; pw += q
    lx /= n - 1; ly /= n - 1
    if pw > 0: px /= pw; py /= pw
    thr = math.radians(theta)
    s2 = sum(c * c * complex(math.cos(2 * m * thr), -math.sin(2 * m * thr)) for m, c in zip(ms, cs))
    Ppi = sum((1 if m % 2 == 0 else -1) * c * c for m, c in zip(ms, cs))
    mbar = sum(m * abs(c) for m, c in zip(ms, cs)) / sum(abs(c) for c in cs)
    prev = None; tot = 0.0
    for i in range(n):
        a, b = pulse(i * 360 / (n - 1), theta, ms, cs, c0, 0.999); ang = math.atan2(b, a)
        if prev is not None:
            d = ang - prev
            while d > math.pi: d -= 2 * math.pi
            while d < -math.pi: d += 2 * math.pi
            tot += d
        prev = ang
    wind = round(tot / (2 * math.pi))
    return dict(theta=theta, phi=phi, M=M, odd=odd, even=even, count=len(ms), n=n, R=R, remax=max(re), remin=min(re), immax=max(im), immin=min(im),
                lin=(lx, ly), pw=(px, py), c0R=c0 / R, PpiR2=Ppi / R / R, s2R2=s2 / R / R, mbar=mbar, wind=wind)

if __name__ == "__main__":
    print("theta, phi, M, odd, even, count, n_points, R, max Re, min Re, max Im, min Im, lin(x,y), pow(x,y), c0/R, P(pi)/R2, S2/R2, <m>, wind")
    for (theta, phi, M, odd, even) in [(180, 10, 1, 1, 1), (180, 10, 10, 1, 1), (180, 10, 100, 1, 1),
                                        (180, 10, 10, 1, 0), (180, 10, 10, 0, 1), (180, 10, 100, 1, 0), (180, 10, 100, 0, 1),
                                        (90, 10, 10, 1, 1), (0, 10, 10, 1, 1), (180, 30, 10, 1, 1)]:
        d = run(theta, phi, M, bool(odd), bool(even))
        print(f"{theta}, {phi}, {M}, {int(odd)}, {int(even)}, {d['count']}, {d['n']}, {d['R']:.4f}, {d['remax']:.4f}, {d['remin']:.4f}, {d['immax']:.4f}, {d['immin']:.4f}, "
              f"({d['lin'][0]:+.4f},{d['lin'][1]:+.4f}), ({d['pw'][0]:+.4f},{d['pw'][1]:+.4f}), {d['c0R']:.4f}, {d['PpiR2']:+.4f}, {d['s2R2'].real:+.4f}{d['s2R2'].imag:+.4f}i, {d['mbar']:.4f}, {d['wind']}")
