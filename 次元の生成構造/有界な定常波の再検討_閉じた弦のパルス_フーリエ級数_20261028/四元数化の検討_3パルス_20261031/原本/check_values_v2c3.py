#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""index_v2c3.html（3 パルス）と同じ式を独立に計算する確認用（check_values_v2c2.py の 3 パルス化、2026-10-29）。

各パルスに内部位相 φ₀_i（位置とは独立）と、Auto の二つのステップ φ_i（位置）・δ_i（内部位相）。
和 Z = Σ_i e^{iφ₀i} Z_i、和の係数 c_m = Σ_i c⁽ⁱ⁾_m e^{i(φ₀i−mθ_i)}、親のノルム R² = Σ_i R_i² + 2 Σ_{i<j} Σ_m c⁽ⁱ⁾c⁽ʲ⁾ cos(mΔ_ij − Φ_ij)
対ごとの速度ベクトル v_ij = (φ_i−φ_j, δ_i−δ_j)、傾き = (δ_i−δ_j)/(φ_i−φ_j)、合成則 v₁₂ = v₁₃ − v₂₃
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

def run(ps, f0=(0, 0, 0), steps=None):
    cs = [coeffs(phi, M, o, e) for (th, phi, M, o, e) in ps]
    ths = [math.radians(th) for (th, *_r) in ps]
    gs = [math.radians(f) for f in f0]
    ms = sorted(set().union(*[set(c) for c, _ in cs]))
    A = {m: sum(c.get(m, 0) * cmath.exp(1j * (g - m * t)) for (c, _), t, g in zip(cs, ths, gs)) for m in ms}
    Rs = [math.sqrt(sum(v * v for v in c.values())) for c, _ in cs]
    pairs = []
    for i in range(len(ps)):
        for j in range(i + 1, len(ps)):
            ci, cj = cs[i][0], cs[j][0]
            cross = sum(ci.get(m, 0) * cj.get(m, 0) * math.cos(m * (ths[i] - ths[j]) - (gs[i] - gs[j])) for m in ms)
            cosT = cross / (Rs[i] * Rs[j]) if Rs[i] * Rs[j] > 0 else float('nan')
            d = ((ps[i][0] - ps[j][0] + 180) % 360 + 360) % 360 - 180
            f = ((f0[i] - f0[j] + 180) % 360 + 360) % 360 - 180
            pairs.append((i, j, d, cross, cosT, f))
    R2 = sum(r * r for r in Rs) + 2 * sum(q[3] for q in pairs); R = math.sqrt(R2)
    assert abs(R2 - sum(abs(v) ** 2 for v in A.values())) < 1e-12
    pt = [(math.cos(t), math.sin(t)) for t in ths]
    area2 = (pt[1][0] - pt[0][0]) * (pt[2][1] - pt[0][1]) - (pt[2][0] - pt[0][0]) * (pt[1][1] - pt[0][1])
    orient = 0 if abs(area2) < 1e-9 else (1 if area2 > 0 else -1)
    c0 = sum(c0 for _, c0 in cs)
    kmax = ms[-1]; n = min(40001, max(1441, 8 * kmax + 1))
    pts = [i * 360 / (n - 1) for i in range(n)]
    def rot(z, g): return (z[0] * math.cos(g) - z[1] * math.sin(g), z[0] * math.sin(g) + z[1] * math.cos(g))
    vals = [tuple(map(sum, zip(*[rot(pulse(p_, th, c, c0i), g) for (th, *_r), (c, c0i), g in zip(ps, cs, gs)]))) for p_ in pts]
    re = [v[0] for v in vals]; im = [v[1] for v in vals]
    S2 = sum(v * v for v in A.values()); Ppi = sum((1 if m % 2 == 0 else -1) * abs(v) ** 2 for m, v in A.items())
    mbar = sum(m * abs(v) for m, v in A.items()) / sum(abs(v) for v in A.values())
    prev = None; tot = 0.0
    for p_ in pts:
        zz = [rot(pulse(p_, th, c, c0i, 0.999), g) for (th, *_r), (c, c0i), g in zip(ps, cs, gs)]; a = sum(z[0] for z in zz); b = sum(z[1] for z in zz)
        ang = math.atan2(b, a)
        if prev is not None:
            d = ang - prev
            while d > math.pi: d -= 2 * math.pi
            while d < -math.pi: d += 2 * math.pi
            tot += d
        prev = ang
    out = dict(Rs=Rs, R=R, ratio=R / sum(Rs), pairs=pairs, orient=orient, remax=max(re), remin=min(re), immax=max(im), immin=min(im),
               c0R=c0 / R, PpiR2=Ppi / R2, s2R2=S2 / R2, mbar=mbar, wind=round(tot / (2 * math.pi)), n=n)
    if steps is not None:   # Auto のステップ (φ_i, δ_i) から対の速度ベクトルと傾き、合成則
        (p1, d1), (p2, d2), (p3, d3) = steps
        v12, v23, v13 = (p1 - p2, d1 - d2), (p2 - p3, d2 - d3), (p1 - p3, d1 - d3)
        sl = lambda v: ('inf' if v[0] == 0 else f"{v[1] / v[0]:.3f}")
        out['vel'] = dict(v12=v12, v23=v23, v13=v13, s12=sl(v12), s23=sl(v23), s13=sl(v13), comp=(v12[0] - (v13[0] - v23[0]), v12[1] - (v13[1] - v23[1])),
                          clock_ratio=('—' if d2 - d3 == 0 else f"{(d1 - d3) / (d2 - d3):.4f}"), vel_ratio=('—' if p2 - p3 == 0 else f"{(p1 - p3) / (p2 - p3):.4f}"), mirror=d2 + d3)
    return out

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
    cases += [
        ("鏡映配置＋対の内部位相 ±45°（θ = 180, 200, 160、φ₀ = 0, +45, −45）：Σc² が実", [(180, 10, 10, 1, 1), (200, 10, 10, 1, 1), (160, 10, 10, 1, 1)], (0, 45, -45)),
        ("同じ配置で対の内部位相が同符号（φ₀ = 0, +45, +45）：Σc² は実でない", [(180, 10, 10, 1, 1), (200, 10, 10, 1, 1), (160, 10, 10, 1, 1)], (0, 45, 45)),
    ]
    print("label | R1,R2,R3 | R, R/ΣR | (Δ12, cosΘ12, Φ12), (Δ13, cosΘ13, Φ13), (Δ23, cosΘ23, Φ23) | 向き | max Re, min Re, max Im | c0/R, P(pi)/R2, S2/R2, <m>, wind")
    for case in cases:
        label, ps = case[0], case[1]; f0 = case[2] if len(case) > 2 else (0, 0, 0)
        d = run(ps, f0)
        pr = ", ".join(f"({q[2]:+.0f}°, {q[4]:+.4f}, {q[5]:+.0f}°)" for q in d['pairs'])
        print(f"{label} | {', '.join(f'{r:.4f}' for r in d['Rs'])} | {d['R']:.4f}, {d['ratio']:.4f} | {pr} | {d['orient']:+d} | {d['remax']:.4f}, {d['remin']:.4f}, {d['immax']:.4f} | "
              f"{d['c0R']:.4f}, {d['PpiR2']:+.4f}, {d['s2R2'].real:+.4f}{d['s2R2'].imag:+.4f}i, {d['mbar']:.4f}, {d['wind']}")
    print()
    print("Auto のステップと対の速度ベクトル（既定：φ = 20, 10, 10、δ = 20, 0, 10）")
    v = run([(180, 10, 10, 1, 1)] * 3, (0, 0, 0), steps=[(20, 20), (10, 0), (10, 10)])['vel']
    print(f"  v12 = {v['v12']} 傾き {v['s12']}（内側）、v23 = {v['v23']} 傾き {v['s23']}（縦軸の上）、v13 = {v['v13']} 傾き {v['s13']}（光円錐の上）")
    print(f"  合成則 v12 − (v13 − v23) = {v['comp']}、内部時計の比 (δ1−δ3)/(δ2−δ3) = {v['clock_ratio']}、速度の比 (φ1−φ3)/(φ2−φ3) = {v['vel_ratio']}、鏡映対の条件 δ2+δ3 = {v['mirror']}")
    v = run([(180, 10, 10, 1, 1)] * 3, (0, 0, 0), steps=[(20, 20), (10, 20), (10, 10)])['vel']
    print(f"  対 13 が光的・対 23 が傾き ∞ のとき対 12：v12 = {v['v12']} 傾き {v['s12']}（光的でない：光円錐は対の性質で基準の取り替えで保たれない）")
