#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""index_v2c3.html（3 パルス）と同じ式を独立に計算する確認用（check_values_v2c2.py の 3 パルス化、2026-10-29）。

和 Z = Σ_i Z_i、和の係数 c_m = Σ_i c⁽ⁱ⁾_m e^{−imθ_i}、親のノルム R² = Σ_i R_i² + 2 Σ_{i<j} R_i R_j cos Θ_ij
出す量：対ごとの Δ_ij・cos Θ_ij、三角形の向き、R_i、R、R/ΣR_i、max/min Re・Im、c₀/R、P(π)/R²、Σc²/R²、⟨m⟩、巻き数
対照：パルス 2・3 を空にすると単独に一致、同じパルス 3 つを Δ = 0 で重ねると和は単独の 3 倍。
鏡映配置（軸上 1 個 ＋ 鏡映対）で Σc² が実になることを確かめる。
"""
import math, cmath

def coeffs(phi, M, odd=True, even=True):
    if M <= 0: return {}, 0.0
    ms = [m for m in range(1, M + 1) if (m % 2 == 1 and odd) or (m % 2 == 0 and even)]
    return {m: (2 / math.pi) * math.sin(m * math.pi * phi / 360) / m for m in ms}, phi / 360

def pulse(psi, theta, c, c0, r=1.0):
    x = math.pi * (psi - theta) / 180
    return c0 + sum(v * r ** m * math.cos(m * x) for m, v in c.items()), sum(v * r ** m * math.sin(m * x) for m, v in c.items())

def run(ps):
    cs = [coeffs(phi, M, o, e) for (th, phi, M, o, e) in ps]
    ths = [math.radians(th) for (th, *_r) in ps]
    ms = sorted(set().union(*[set(c) for c, _ in cs]))
    A = {m: sum(c.get(m, 0) * cmath.exp(-1j * m * t) for (c, _), t in zip(cs, ths)) for m in ms}
    Rs = [math.sqrt(sum(v * v for v in c.values())) for c, _ in cs]
    pairs = []
    for i in range(len(ps)):
        for j in range(i + 1, len(ps)):
            ci, cj = cs[i][0], cs[j][0]
            cross = sum(ci.get(m, 0) * cj.get(m, 0) * math.cos(m * (ths[i] - ths[j])) for m in ms)
            cosT = cross / (Rs[i] * Rs[j]) if Rs[i] * Rs[j] > 0 else float('nan')
            d = ((ps[i][0] - ps[j][0] + 180) % 360 + 360) % 360 - 180
            pairs.append((i, j, d, cross, cosT))
    R2 = sum(r * r for r in Rs) + 2 * sum(q[3] for q in pairs); R = math.sqrt(R2)
    assert abs(R2 - sum(abs(v) ** 2 for v in A.values())) < 1e-12
    pt = [(math.cos(t), math.sin(t)) for t in ths]
    area2 = (pt[1][0] - pt[0][0]) * (pt[2][1] - pt[0][1]) - (pt[2][0] - pt[0][0]) * (pt[1][1] - pt[0][1])
    orient = 0 if abs(area2) < 1e-9 else (1 if area2 > 0 else -1)
    c0 = sum(c0 for _, c0 in cs)
    kmax = ms[-1]; n = min(40001, max(1441, 8 * kmax + 1))
    pts = [i * 360 / (n - 1) for i in range(n)]
    vals = [tuple(map(sum, zip(*[pulse(p_, th, c, c0i) for (th, *_r), (c, c0i) in zip(ps, cs)]))) for p_ in pts]
    re = [v[0] for v in vals]; im = [v[1] for v in vals]
    S2 = sum(v * v for v in A.values()); Ppi = sum((1 if m % 2 == 0 else -1) * abs(v) ** 2 for m, v in A.items())
    mbar = sum(m * abs(v) for m, v in A.items()) / sum(abs(v) for v in A.values())
    prev = None; tot = 0.0
    for p_ in pts:
        a = sum(pulse(p_, th, c, c0i, 0.999)[0] for (th, *_r), (c, c0i) in zip(ps, cs)); b = sum(pulse(p_, th, c, c0i, 0.999)[1] for (th, *_r), (c, c0i) in zip(ps, cs))
        ang = math.atan2(b, a)
        if prev is not None:
            d = ang - prev
            while d > math.pi: d -= 2 * math.pi
            while d < -math.pi: d += 2 * math.pi
            tot += d
        prev = ang
    return dict(Rs=Rs, R=R, ratio=R / sum(Rs), pairs=pairs, orient=orient, remax=max(re), remin=min(re), immax=max(im), immin=min(im),
                c0R=c0 / R, PpiR2=Ppi / R2, s2R2=S2 / R2, mbar=mbar, wind=round(tot / (2 * math.pi)), n=n)

if __name__ == "__main__":
    cases = [
        ("対照：パルス 2・3 なし（単独と一致）", [(180, 10, 10, 1, 1), (180, 10, 0, 1, 1), (180, 10, 0, 1, 1)]),
        ("対照：同じパルス 3 つ Δ=0（単独の 3 倍）", [(180, 10, 10, 1, 1)] * 3),
        ("三角形 θ = 180, 150, 100", [(180, 10, 10, 1, 1), (150, 10, 10, 1, 1), (100, 10, 10, 1, 1)]),
        ("三角形 θ = 180, 100, 150（向きが逆）", [(180, 10, 10, 1, 1), (100, 10, 10, 1, 1), (150, 10, 10, 1, 1)]),
        ("鏡映配置 θ = 180, 200, 160（軸上 1 個＋鏡映対、基準 ψ=0 は軸）", [(180, 10, 10, 1, 1), (200, 10, 10, 1, 1), (160, 10, 10, 1, 1)]),
        ("鏡映配置 θ = 170, 190, 150（軸が 170、基準 ψ=0 は軸でない）", [(170, 10, 10, 1, 1), (190, 10, 10, 1, 1), (150, 10, 10, 1, 1)]),
        ("違うパルス M = 10, 10, 100、θ = 180, 162, 144", [(180, 10, 10, 1, 1), (162, 10, 10, 1, 1), (144, 10, 100, 1, 1)]),
    ]
    print("label | R1,R2,R3 | R, R/ΣR | (Δ12, cosΘ12), (Δ13, cosΘ13), (Δ23, cosΘ23) | 向き | max Re, min Re, max Im | c0/R, P(pi)/R2, S2/R2, <m>, wind")
    for label, ps in cases:
        d = run(ps)
        pr = ", ".join(f"({q[2]:+.0f}°, {q[4]:+.4f})" for q in d['pairs'])
        print(f"{label} | {', '.join(f'{r:.4f}' for r in d['Rs'])} | {d['R']:.4f}, {d['ratio']:.4f} | {pr} | {d['orient']:+d} | {d['remax']:.4f}, {d['remin']:.4f}, {d['immax']:.4f} | "
              f"{d['c0R']:.4f}, {d['PpiR2']:+.4f}, {d['s2R2'].real:+.4f}{d['s2R2'].imag:+.4f}i, {d['mbar']:.4f}, {d['wind']}")
