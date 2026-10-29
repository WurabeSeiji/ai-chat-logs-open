#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""index_v2c2.html と同じ式を独立に計算する確認用（check_values_v2c.py を 2 パルスにしたもの、2026-10-29）。

和 Z = Z₁ + Z₂、Z_i(ψ) = φ_i/360 + (2/π) Σ_m [ sin(mπφ_i/360) / m ] · e^{imπ(ψ−θ_i)/180}（m は 1..M_i から奇・偶で選ぶ）
各パルスに内部位相 φ₀_i（位置とは独立）。和の係数 c_m = c⁽¹⁾_m e^{i(φ₀₁−mθ₁)} + c⁽²⁾_m e^{i(φ₀₂−mθ₂)}、
親のノルム R² = Σ|c_m|² = R₁² + R₂² + 2 Σ c⁽¹⁾c⁽²⁾ cos(mΔ − Φ)、Φ = φ₀₁ − φ₀₂
出す量：max/min Re・Im（和）、R₁、R₂、R、R/(R₁+R₂)、cos Θ、Θ、tan α = R₂/R₁、α、線形重心、パワー重心、
       c₀/R、P(π)/R² = (Σ_偶 |c|² − Σ_奇 |c|²)/R²、Σ c_m²/R²（基準 ψ = 0）、⟨m⟩ = Σ m|c_m|/Σ|c_m|、巻き数（|w| = 0.999）
対照：パルス 2 を空にすると単独（check_values_v2c.py）と一致、同じパルスを Δ = 0 で重ねると和は単独の 2 倍。
ψ の刻みは index_v2c2.html と同じ（n = min(40001, max(1441, 8·kmax + 1))）。
"""
import math, cmath

def coeffs(phi, M, odd=True, even=True):
    ms = [m for m in range(1, M + 1) if (m % 2 == 1 and odd) or (m % 2 == 0 and even)]
    return {m: (2 / math.pi) * math.sin(m * math.pi * phi / 360) / m for m in ms}, phi / 360

def pulse(psi, theta, c, c0, r=1.0):
    x = math.pi * (psi - theta) / 180
    return c0 + sum(v * r ** m * math.cos(m * x) for m, v in c.items()), sum(v * r ** m * math.sin(m * x) for m, v in c.items())

def run(p1, p2, f1=0.0, f2=0.0):
    (th1, phi1, M1, o1, e1), (th2, phi2, M2, o2, e2) = p1, p2
    g1, g2 = math.radians(f1), math.radians(f2)
    c1, c01 = coeffs(phi1, M1, o1, e1); c2, c02 = (coeffs(phi2, M2, o2, e2) if M2 > 0 else ({}, 0.0))
    t1, t2 = math.radians(th1), math.radians(th2)
    ms = sorted(set(c1) | set(c2))
    A = {m: c1.get(m, 0) * cmath.exp(1j * (g1 - m * t1)) + c2.get(m, 0) * cmath.exp(1j * (g2 - m * t2)) for m in ms}
    R1 = math.sqrt(sum(v * v for v in c1.values())); R2 = math.sqrt(sum(v * v for v in c2.values()))
    cross = sum(c1.get(m, 0) * c2.get(m, 0) * math.cos(m * (t1 - t2) - (g1 - g2)) for m in ms)
    R2tot = max(0.0, R1 * R1 + R2 * R2 + 2 * cross); R = math.sqrt(R2tot)   # 完全に打ち消すと丸めで僅かに負になる
    assert abs(R2tot - sum(abs(v) ** 2 for v in A.values())) < 1e-12
    cosT = cross / (R1 * R2) if R1 * R2 > 0 else float('nan'); T = math.degrees(math.acos(max(-1, min(1, cosT)))) if R2 > 0 else float('nan')
    alpha = math.degrees(math.atan2(R2, R1)); c0 = c01 + c02
    kmax = ms[-1]; n = min(40001, max(1441, 8 * kmax + 1))
    pts = [(i * 360 / (n - 1)) for i in range(n)]
    def rot(z, g): return (z[0] * math.cos(g) - z[1] * math.sin(g), z[0] * math.sin(g) + z[1] * math.cos(g))
    vals = [tuple(a + b for a, b in zip(rot(pulse(ps, th1, c1, c01), g1), rot(pulse(ps, th2, c2, c02), g2))) for ps in pts]
    re = [v[0] for v in vals]; im = [v[1] for v in vals]
    lx = ly = px = py = pw = 0.0
    for i in range(n - 1):
        a = R + re[i]; q = re[i] ** 2 + im[i] ** 2; cs = math.cos(2 * math.pi * i / (n - 1)); sn = math.sin(2 * math.pi * i / (n - 1))
        lx += a * cs; ly += a * sn; px += q * R * cs; py += q * R * sn; pw += q
    lx /= n - 1; ly /= n - 1
    if pw > 0: px /= pw; py /= pw
    S2 = sum(v * v for v in A.values()); Ppi = sum((1 if m % 2 == 0 else -1) * abs(v) ** 2 for m, v in A.items())
    sab = sum(abs(v) for v in A.values()); mbar = sum(m * abs(v) for m, v in A.items()) / sab if sab > 0 else float('nan')
    if R2tot < 1e-15:   # 完全に打ち消した場合（Φ = 180°、Δ = 0 の同じパルス）：R で割る量は定義できない
        return dict(n=n, R1=R1, R2=R2, R=0.0, ratio=0.0, cosT=cosT, T=T, alpha=alpha, remax=max(re), remin=min(re), immax=max(im), immin=min(im), lin=(lx, ly), pw=(float('nan'), float('nan')), c0R=float('nan'), PpiR2=float('nan'), s2R2=complex(float('nan'), float('nan')), mbar=mbar, wind=0)
    prev = None; tot = 0.0
    for ps in pts:
        a1, b1 = rot(pulse(ps, th1, c1, c01, 0.999), g1); a2, b2 = rot(pulse(ps, th2, c2, c02, 0.999), g2); ang = math.atan2(b1 + b2, a1 + a2)
        if prev is not None:
            d = ang - prev
            while d > math.pi: d -= 2 * math.pi
            while d < -math.pi: d += 2 * math.pi
            tot += d
        prev = ang
    wind = round(tot / (2 * math.pi))
    return dict(n=n, R1=R1, R2=R2, R=R, ratio=R / (R1 + R2), cosT=cosT, T=T, alpha=alpha, remax=max(re), remin=min(re), immax=max(im), immin=min(im),
                lin=(lx, ly), pw=(px, py), c0R=c0 / R, PpiR2=Ppi / R2tot, s2R2=S2 / R2tot, mbar=mbar, wind=wind)

if __name__ == "__main__":
    print("label | n | R1, R2, R, R/(R1+R2) | cosT, Theta, alpha | max Re, min Re, max Im, min Im | lin, pow | c0/R, P(pi)/R2, S2/R2, <m>, wind")
    cases = [
        ("対照：パルス 2 なし（単独と一致するはず）", (180, 10, 10, 1, 1), (180, 10, 0, 1, 1)),
        ("対照：同じパルス Δ=0（単独の 2 倍）", (180, 10, 10, 1, 1), (180, 10, 10, 1, 1)),
        ("同じパルス Δ=10", (180, 10, 10, 1, 1), (170, 10, 10, 1, 1)),
        ("同じパルス Δ=18", (180, 10, 10, 1, 1), (162, 10, 10, 1, 1)),
        ("同じパルス Δ=90", (180, 10, 10, 1, 1), (90, 10, 10, 1, 1)),
        ("同じパルス Δ=180", (180, 10, 10, 1, 1), (0, 10, 10, 1, 1)),
        ("同じパルス 基準真ん中 Δ=18（θ=189, 171）", (189, 10, 10, 1, 1), (171, 10, 10, 1, 1)),
        ("違うパルス M=10, 100 Δ=0", (180, 10, 10, 1, 1), (180, 10, 100, 1, 1)),
        ("違うパルス M=10, 100 Δ=18", (180, 10, 10, 1, 1), (162, 10, 100, 1, 1)),
        ("奇のみと偶のみ Δ=0", (180, 10, 10, 1, 0), (180, 10, 10, 0, 1)),
        ("奇のみと偶のみ Δ=90", (180, 10, 10, 1, 0), (90, 10, 10, 0, 1)),
    ]
    cases += [
        ("内部位相 Φ=90°、Δ=0（同じパルス）", (180, 10, 10, 1, 1), (180, 10, 10, 1, 1), 90, 0),
        ("内部位相 Φ=180°、Δ=0（同じパルス、打ち消し）", (180, 10, 10, 1, 1), (180, 10, 10, 1, 1), 180, 0),
        ("内部位相 Φ=90°、Δ=18（同じパルス）", (180, 10, 10, 1, 1), (162, 10, 10, 1, 1), 90, 0),
        ("内部位相 Φ=180°、Δ=36（同じパルス：mΔ−Φ=±180 は m=10）", (180, 10, 10, 1, 1), (144, 10, 10, 1, 1), 180, 0),
    ]
    for case in cases:
        label, p1, p2 = case[:3]; f = case[3:] if len(case) > 3 else ()
        d = run(p1, p2, *f)
        print(f"{label} | {d['n']} | {d['R1']:.4f}, {d['R2']:.4f}, {d['R']:.4f}, {d['ratio']:.4f} | {d['cosT']:+.4f}, {d['T']:.3f}, {d['alpha']:.3f} | "
              f"{d['remax']:.4f}, {d['remin']:.4f}, {d['immax']:.4f}, {d['immin']:.4f} | ({d['lin'][0]:+.4f},{d['lin'][1]:+.4f}), ({d['pw'][0]:+.4f},{d['pw'][1]:+.4f}) | "
              f"{d['c0R']:.4f}, {d['PpiR2']:+.4f}, {d['s2R2'].real:+.4f}{d['s2R2'].imag:+.4f}i, {d['mbar']:.4f}, {d['wind']}")
