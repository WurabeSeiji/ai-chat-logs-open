#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""index_v2q3.html（四元数 3 パルス）と同じ式を独立に計算する確認用（原本 check_values_v2c3.py の四元数化、2026-10-31）。

四元数は [w, x, y, z]（x ↔ i 円周、y ↔ j 半径、z ↔ k 垂直）。各パルスに虚部の向き n_i（単位純虚四元数、既定は k）。
各パルスは平面 span(1, n_i) ≅ ℂ の中にあるので、式は原本の i を n_i に置き換えただけ：
  和 Q = Σ_i e^{n_i φ₀i} Q_i、和の係数 c_m = Σ_i c⁽ⁱ⁾_m e^{n_i(φ₀i−mθ_i)} ∈ ℍ
  親のノルム R² = Σ_m |c_m|² = Σ_i R_i² + 2 Σ_{i<j} Σ_m c⁽ⁱ⁾c⁽ʲ⁾ Sc[e^{n_i α_i} e^{−n_j α_j}]、α_i = φ₀i − mθ_i
  Sc[e^{n_i a} e^{−n_j b}]  = cos a cos b + (n_i·n_j) sin a sin b                      （n_i = n_j で cos(a − b) = cos(mΔ − Φ)）
  Vec[e^{n_i a} e^{−n_j b}] = sin a cos b n_i − cos a sin b n_j − sin a sin b (n_i×n_j)  （n_i = n_j で sin(a − b) n）
対照：n_i = k（全部）で原本 check_values_v2c3.py（原本/ に無変更コピー、import）の全ケースと機械的に突き合わせる。
新しい量：n_i·n_j、n_i×n_j、対の位置ベクトル Δ⃗_ij（e^{n_iθ_i} e^{−n_jθ_j} の角 × 軸）、交差項のベクトル部 V⃗_ij、Σc² のベクトル部、
位置ベクトルの持ち上げ x⃗_i = Θ_i n_i と ±|x⃗_io|（右パネルの横軸）。巻き数は n が全部同じときだけ定義する。
振幅 A_i（ps の 6 番目、省略時 1）：係数 c⁽ⁱ⁾_m と直流 c⁽ⁱ⁾₀ を A_i 倍する（Q_i → A_i Q_i、R_i → A_i R_i）。共通の A では比の量が不変であることを確かめる。
"""
import math, cmath, os, importlib.util

def coeffs(phi, M, odd=True, even=True):
    if M <= 0: return {}, 0.0
    ms = [m for m in range(1, M + 1) if (m % 2 == 1 and odd) or (m % 2 == 0 and even)]
    return {m: (2 / math.pi) * math.sin(m * math.pi * phi / 360) / m for m in ms}, phi / 360

def pulse(psi, theta, c, c0, r=1.0):
    x = math.pi * (psi - theta) / 180
    return c0 + sum(v * r ** m * math.cos(m * x) for m, v in c.items()), sum(v * r ** m * math.sin(m * x) for m, v in c.items())

# ---- 四元数（index_v2q3.html と同じ。足すのは成分ごと、要るのは掛け算一つ）
K = (0.0, 0.0, 1.0)
def unit(n):
    L = math.sqrt(sum(v * v for v in n)); return K if L < 1e-12 else tuple(v / L for v in n)
def dot(a, b): return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
def crs(a, b): return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
def qmul(p, q):
    return (p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3], p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2],
            p[0] * q[2] - p[1] * q[3] + p[2] * q[0] + p[3] * q[1], p[0] * q[3] + p[1] * q[2] - p[2] * q[1] + p[3] * q[0])
def qexp(n, a): s = math.sin(a); return (math.cos(a), n[0] * s, n[1] * s, n[2] * s)          # e^{n a} = cos a + n sin a
def qconj(q): return (q[0], -q[1], -q[2], -q[3])
def qnorm2(q): return sum(v * v for v in q)
def wrap(d): return ((d + 180) % 360 + 360) % 360 - 180
def v3(v, d=4): return "(" + ", ".join(f"{x:+.{d}f}" for x in v) + ")"

def run(ps, f0=(0, 0, 0), steps=None, ns=None):
    cs = []
    for p in ps:   # (θ, φ, M, 奇, 偶[, A])。振幅 A は係数と直流を A 倍（省略時 1 で原本と同じ）
        c, c0 = coeffs(*p[1:5]); A = p[5] if len(p) > 5 else 1.0
        cs.append(({m: A * v for m, v in c.items()}, A * c0))
    ths = [math.radians(th) for (th, *_r) in ps]
    gs = [math.radians(f) for f in f0]
    nvs = [unit(n) for n in (ns if ns is not None else [K] * len(ps))]
    ms = sorted(set().union(*[set(c) for c, _ in cs]))
    A = {}
    for m in ms:   # c_m = Σ_i c⁽ⁱ⁾_m e^{n_i(φ₀i − mθ_i)}（四元数）
        v = [0.0, 0.0, 0.0, 0.0]
        for (c, _), t, g, nv in zip(cs, ths, gs, nvs):
            cm = c.get(m, 0.0); a = g - m * t; s = cm * math.sin(a)
            v[0] += cm * math.cos(a); v[1] += s * nv[0]; v[2] += s * nv[1]; v[3] += s * nv[2]
        A[m] = tuple(v)
    Rs = [math.sqrt(sum(v * v for v in c.values())) for c, _ in cs]
    pairs = []
    for i in range(len(ps)):
        for j in range(i + 1, len(ps)):
            ci, cj = cs[i][0], cs[j][0]; ni, nj = nvs[i], nvs[j]; nd = dot(ni, nj); nx = crs(ni, nj)
            cross = 0.0; cross2 = 0.0; Vc = [0.0, 0.0, 0.0]; Vc2 = [0.0, 0.0, 0.0]
            for m in ms:
                w = ci.get(m, 0.0) * cj.get(m, 0.0)
                if w == 0.0: continue
                a, b = gs[i] - m * ths[i], gs[j] - m * ths[j]
                g = qmul(qexp(ni, a), qconj(qexp(nj, b)))                      # 掛け算で
                cross += w * g[0]; Vc[0] += w * g[1]; Vc[1] += w * g[2]; Vc[2] += w * g[3]
                cross2 += w * (math.cos(a) * math.cos(b) + nd * math.sin(a) * math.sin(b))   # 閉じた式で（照合）
                for k in range(3): Vc2[k] += w * (math.sin(a) * math.cos(b) * ni[k] - math.cos(a) * math.sin(b) * nj[k] - math.sin(a) * math.sin(b) * nx[k])
            assert abs(cross - cross2) < 1e-12 and max(abs(Vc[k] - Vc2[k]) for k in range(3)) < 1e-12   # スカラー部・ベクトル部の閉じた式
            cosT = cross / (Rs[i] * Rs[j]) if Rs[i] * Rs[j] > 0 else float('nan')
            d = wrap(ps[i][0] - ps[j][0]); f = wrap(f0[i] - f0[j])
            g0 = qmul(qexp(ni, ths[i]), qconj(qexp(nj, ths[j])))               # 対の位置ベクトル：e^{n_iθ_i} e^{−n_jθ_j} の角 × 軸
            Dabs = math.degrees(math.acos(max(-1.0, min(1.0, g0[0])))); gv = math.sqrt(g0[1] ** 2 + g0[2] ** 2 + g0[3] ** 2)
            Dvec = (0.0, 0.0, 0.0) if gv < 1e-12 else tuple(x / gv * Dabs for x in g0[1:])
            pairs.append(dict(i=i, j=j, d=d, cross=cross, cosT=cosT, f=f, nd=nd, nx=nx, Vc=tuple(Vc), Dvec=Dvec, Dabs=Dabs))
    R2 = sum(r * r for r in Rs) + 2 * sum(q['cross'] for q in pairs); R = math.sqrt(R2)
    assert abs(R2 - sum(qnorm2(v) for v in A.values())) < 1e-12                 # R² = Σ|c_m|²（四元数ノルム）と交差項の式の一致
    pt = [(math.cos(t), math.sin(t)) for t in ths]
    area2 = (pt[1][0] - pt[0][0]) * (pt[2][1] - pt[0][1]) - (pt[2][0] - pt[0][0]) * (pt[1][1] - pt[0][1])
    orient = 0 if abs(area2) < 1e-9 else (1 if area2 > 0 else -1)
    c0 = sum(c0 for _, c0 in cs)
    kmax = ms[-1]; n = min(40001, max(1441, 8 * kmax + 1))
    pts = [i * 360 / (n - 1) for i in range(n)]
    def rot(z, g): return (z[0] * math.cos(g) - z[1] * math.sin(g), z[0] * math.sin(g) + z[1] * math.cos(g))
    def Q(p_, r=1.0):   # 和 Q = W + V（各パルスの (re, im) を平面 span(1, n_i) に埋め込んで足す）
        W = 0.0; V = [0.0, 0.0, 0.0]
        for (th, *_r), (c, c0i), g, nv in zip(ps, cs, gs, nvs):
            re, im = rot(pulse(p_, th, c, c0i, r), g); W += re
            for k in range(3): V[k] += im * nv[k]
        return W, V
    vals = [Q(p_) for p_ in pts]
    re = [v[0] for v in vals]; comp = [[v[1][k] for v in vals] for k in range(3)]
    S2 = [0.0, 0.0, 0.0, 0.0]; Ppi = 0.0; sAbs = 0.0; sMAbs = 0.0
    for m, v in A.items():   # c² = (w² − |v|²) + 2 w v
        vv = v[1] ** 2 + v[2] ** 2 + v[3] ** 2; q = v[0] ** 2 + vv; ab = math.sqrt(q)
        S2[0] += v[0] ** 2 - vv; S2[1] += 2 * v[0] * v[1]; S2[2] += 2 * v[0] * v[2]; S2[3] += 2 * v[0] * v[3]
        Ppi += (1 if m % 2 == 0 else -1) * q; sAbs += ab; sMAbs += m * ab
    mbar = sMAbs / sAbs
    same = all(dot(nv, nvs[0]) > 1 - 1e-9 for nv in nvs)
    wind = None
    if same:   # 巻き数：平面を共有するときだけ（原本と同じ計算、虚部は共通の軸の成分）
        prev = None; tot = 0.0
        for p_ in pts:
            W, V = Q(p_, 0.999); ang = math.atan2(dot(V, nvs[0]), W)
            if prev is not None:
                dd = ang - prev
                while dd > math.pi: dd -= 2 * math.pi
                while dd < -math.pi: dd += 2 * math.pi
                tot += dd
            prev = ang
        wind = round(tot / (2 * math.pi))
    out = dict(Rs=Rs, R=R, ratio=R / sum(Rs), pairs=pairs, orient=orient, remax=max(re), remin=min(re),
               vmax=[max(c) for c in comp], vmin=[min(c) for c in comp],
               c0R=c0 / R, PpiR2=Ppi / R2, s2R2=tuple(s / R2 for s in S2), mbar=mbar, wind=wind, n=n, nvs=nvs, same=same)
    if steps is not None:   # Auto のステップ (φ_i, δ_i) から対の速度ベクトルと傾き、合成則（原本と同じ。軸に依らない）
        (p1, d1), (p2, d2), (p3, d3) = steps
        v12, v23, v13 = (p1 - p2, d1 - d2), (p2 - p3, d2 - d3), (p1 - p3, d1 - d3)
        sl = lambda v: ('inf' if v[0] == 0 else f"{v[1] / v[0]:.3f}")
        out['vel'] = dict(v12=v12, v23=v23, v13=v13, s12=sl(v12), s23=sl(v23), s13=sl(v13), comp=(v12[0] - (v13[0] - v23[0]), v12[1] - (v13[1] - v23[1])),
                          clock_ratio=('—' if d2 - d3 == 0 else f"{(d1 - d3) / (d2 - d3):.4f}"), vel_ratio=('—' if p2 - p3 == 0 else f"{(p1 - p3) / (p2 - p3):.4f}"), mirror=d2 + d3)
    return out

def lift(Th, ns, o):
    """右パネルの横軸：位置ベクトル x⃗_i = Θ_i n_i（度）、x⃗_io = x⃗_i − x⃗_o、±|x⃗_io|（符号は n_o への射影）"""
    nvs = [unit(n) for n in ns]
    xs = [tuple(T * v for v in nv) for T, nv in zip(Th, nvs)]
    out = []
    for i in range(len(Th)):
        if i == o: continue
        r = tuple(a - b for a, b in zip(xs[i], xs[o])); L = math.sqrt(dot(r, r))
        out.append((i, r, (-1 if dot(r, nvs[o]) < 0 else 1) * L))
    return out

I, J = (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)
N45 = (math.sin(math.radians(45)), 0.0, math.cos(math.radians(45)))

if __name__ == "__main__":
    base = os.path.dirname(os.path.abspath(__file__))
    cases = [   # 原本 check_values_v2c3.py と同じ 9 例（軸は全部 k）＝対照
        ("対照：パルス 2・3 なし（単独と一致）", [(180, 10, 10, 1, 1), (180, 10, 0, 1, 1), (180, 10, 0, 1, 1)]),
        ("対照：同じパルス 3 つ Δ=0（単独の 3 倍）", [(180, 10, 10, 1, 1)] * 3),
        ("三角形 θ = 180, 150, 100", [(180, 10, 10, 1, 1), (150, 10, 10, 1, 1), (100, 10, 10, 1, 1)]),
        ("三角形 θ = 180, 100, 150（向きが逆）", [(180, 10, 10, 1, 1), (100, 10, 10, 1, 1), (150, 10, 10, 1, 1)]),
        ("鏡映配置 θ = 180, 200, 160（軸上 1 個＋鏡映対、基準 ψ=0 は軸）", [(180, 10, 10, 1, 1), (200, 10, 10, 1, 1), (160, 10, 10, 1, 1)]),
        ("鏡映配置 θ = 170, 190, 150（軸が 170、基準 ψ=0 は軸でない）", [(170, 10, 10, 1, 1), (190, 10, 10, 1, 1), (150, 10, 10, 1, 1)]),
        ("違うパルス M = 10, 10, 100、θ = 180, 162, 144", [(180, 10, 10, 1, 1), (162, 10, 10, 1, 1), (144, 10, 100, 1, 1)]),
        ("鏡映配置＋対の内部位相 ±45°（θ = 180, 200, 160、φ₀ = 0, +45, −45）：Σc² が実", [(180, 10, 10, 1, 1), (200, 10, 10, 1, 1), (160, 10, 10, 1, 1)], (0, 45, -45)),
        ("同じ配置で対の内部位相が同符号（φ₀ = 0, +45, +45）：Σc² は実でない", [(180, 10, 10, 1, 1), (200, 10, 10, 1, 1), (160, 10, 10, 1, 1)], (0, 45, 45)),
    ]
    same3 = [(150, 10, 10, 1, 1)] * 3
    tri = [(180, 10, 10, 1, 1), (150, 10, 10, 1, 1), (100, 10, 10, 1, 1)]
    mir = [(180, 10, 10, 1, 1), (200, 10, 10, 1, 1), (160, 10, 10, 1, 1)]
    qcases = [   # 軸を変えた例（label, ps, f0, ns）
        ("同じパルス 3 つ Δ=0（θ = 150）軸 k, k, k", same3, (0, 0, 0), [K, K, K]),
        ("同じパルス 3 つ Δ=0（θ = 150）軸 i, j, k（直交：交差項は cos·cos だけ）", same3, (0, 0, 0), [I, J, K]),
        ("同じパルス 3 つ Δ=0（θ = 150）軸 k, k, −k（3 は共役の平面）", same3, (0, 0, 0), [K, K, (0, 0, -1)]),
        ("同じパルス 3 つ Δ=0（θ = 150）軸 k, 45°, k", same3, (0, 0, 0), [K, N45, K]),
        ("三角形 θ = 180, 150, 100 軸 i, j, k", tri, (0, 0, 0), [I, J, K]),
        ("三角形 θ = 180, 150, 100 軸 k, 45°, k", tri, (0, 0, 0), [K, N45, K]),
        ("鏡映配置 φ₀ = 0, +45, −45 軸 k, k, k（対照：Σc² 実）", mir, (0, 45, -45), [K, K, K]),
        ("鏡映配置 φ₀ = 0, +45, −45 軸 k, i, −i（対の軸が逆：Σc² は閉じない）", mir, (0, 45, -45), [K, I, (-1, 0, 0)]),
        ("鏡映配置 φ₀ = 0, +45, −45 軸 k, i, i（対が平面を共有：Σc² 実）", mir, (0, 45, -45), [K, I, I]),
        ("鏡映配置 φ₀ = 0, +45, −45 軸 i, j, k（三つとも別の平面）", mir, (0, 45, -45), [I, J, K]),
    ]
    def line(label, d):
        pr = ", ".join(f"({q['d']:+.0f}°, {q['cosT']:+.4f}, {q['f']:+.0f}°)" for q in d['pairs'])
        nn = ", ".join(f"{q['nd']:+.3f}" for q in d['pairs'])
        Vn = ", ".join(v3(tuple(x / (d['Rs'][q['i']] * d['Rs'][q['j']]) for x in q['Vc'])) if d['Rs'][q['i']] * d['Rs'][q['j']] > 0 else "—" for q in d['pairs'])
        Dv = ", ".join(v3(q['Dvec'], 1) for q in d['pairs'])
        s2 = f"{d['s2R2'][0]:+.4f} + {v3(d['s2R2'][1:])}"
        w = "—" if d['wind'] is None else str(d['wind'])
        return (f"{label} | {', '.join(f'{r:.4f}' for r in d['Rs'])} | {d['R']:.4f}, {d['ratio']:.4f} | {pr} | {d['orient']:+d} | "
                f"{d['remax']:.4f}, {d['remin']:.4f}, {v3(d['vmax'])} | {d['c0R']:.4f}, {d['PpiR2']:+.4f}, {s2}, {d['mbar']:.4f}, {w} | n·n: {nn} | V⃗/(RiRj): {Vn} | Δ⃗: {Dv}")
    hdr = "label | R1,R2,R3 | R, R/ΣR | (Δ12, cosΘ12, Φ12), (Δ13, cosΘ13, Φ13), (Δ23, cosΘ23, Φ23) | 向き | max Re, min Re, max(V_i, V_j, V_k) | c0/R, P(pi)/R2, S2/R2 = スカラー + (i, j, k), <m>, wind | n·n (12, 13, 23) | V⃗/(RiRj) | Δ⃗ (度)"
    print("== 軸は全部 k（原本の 9 例と同じ配置）==")
    print(hdr)
    for case in cases:
        label, ps = case[0], case[1]; f0 = case[2] if len(case) > 2 else (0, 0, 0)
        print(line(label, run(ps, f0)))
    print()
    print("== 軸を変えた例 ==")
    print(hdr)
    for label, ps, f0, ns in qcases:
        print(line(label, run(ps, f0, ns=ns)))
    print()
    print("== 振幅を変えた例（A_i：係数と直流を A_i 倍）==")
    print(hdr)
    triA = [(180, 10, 10, 1, 1, 2.0), (150, 10, 10, 1, 1, 1.0), (100, 10, 10, 1, 1, 0.5)]
    acases = [
        ("三角形 θ = 180, 150, 100、A = 2, 1, 0.5、軸 k, k, k", triA, (0, 0, 0), [K, K, K]),
        ("三角形 θ = 180, 150, 100、A = 2, 1, 0.5、軸 i, j, k", triA, (0, 0, 0), [I, J, K]),
        ("同じパルス 3 つ Δ=0（θ = 150）、A = 1, 1, 0（3 つ目を消す ＝ 2 パルス）、軸 k, k, k", [(150, 10, 10, 1, 1, 1.0), (150, 10, 10, 1, 1, 1.0), (150, 10, 10, 1, 1, 0.0)], (0, 0, 0), [K, K, K]),
        ("三角形 θ = 180, 150, 100、A = 2, 2, 2（共通の倍率）", [(180, 10, 10, 1, 1, 2.0), (150, 10, 10, 1, 1, 2.0), (100, 10, 10, 1, 1, 2.0)], (0, 0, 0), [K, K, K]),
    ]
    for label, ps, f0, ns in acases:
        print(line(label, run(ps, f0, ns=ns)))
    d1 = run(tri); d2 = run([(th, phi, M, o, e, 2.0) for (th, phi, M, o, e) in tri])   # 共通の A = 2：R は 2 倍、比の量は不変
    inv = max(abs(d1['ratio'] - d2['ratio']), abs(d1['c0R'] - d2['c0R']), abs(d1['PpiR2'] - d2['PpiR2']), abs(d1['mbar'] - d2['mbar']),
              max(abs(a - b) for a, b in zip(d1['s2R2'], d2['s2R2'])), max(abs(q1['cosT'] - q2['cosT']) for q1, q2 in zip(d1['pairs'], d2['pairs'])))
    print(f"  共通の A = 2：R {d1['R']:.4f} → {d2['R']:.4f}（比 {d2['R'] / d1['R']:.6f}）、R_i も 2 倍、max Re {d1['remax']:.4f} → {d2['remax']:.4f}。"
          f"比の量（R/ΣR、cos Θ、c₀/R、P(π)/R²、Σc²/R²、⟨m⟩）の最大差 {inv:.1e}、巻き数 {d1['wind']} → {d2['wind']}")
    assert abs(d2['R'] / d1['R'] - 2) < 1e-12 and inv < 1e-12 and d1['wind'] == d2['wind']
    print()
    print("Auto のステップと対の速度ベクトル（既定：φ = 20, 10, 10、δ = 20, 0, 10。軸に依らない）")
    v = run([(180, 10, 10, 1, 1)] * 3, (0, 0, 0), steps=[(20, 20), (10, 0), (10, 10)])['vel']
    print(f"  v12 = {v['v12']} 傾き {v['s12']}（内側）、v23 = {v['v23']} 傾き {v['s23']}（縦軸の上）、v13 = {v['v13']} 傾き {v['s13']}（光円錐の上）")
    print(f"  合成則 v12 − (v13 − v23) = {v['comp']}、内部時計の比 (δ1−δ3)/(δ2−δ3) = {v['clock_ratio']}、速度の比 (φ1−φ3)/(φ2−φ3) = {v['vel_ratio']}、鏡映対の条件 δ2+δ3 = {v['mirror']}")
    v = run([(180, 10, 10, 1, 1)] * 3, (0, 0, 0), steps=[(20, 20), (10, 20), (10, 10)])['vel']
    print(f"  対 13 が光的・対 23 が傾き ∞ のとき対 12：v12 = {v['v12']} 傾き {v['s12']}（光的でない：光円錐は対の性質で基準の取り替えで保たれない）")
    print()
    print("位置ベクトルの持ち上げ x⃗_i = Θ_i n_i と右パネルの横軸 ±|x⃗_io|（観測者 o = 3、Θ = 180, 150, 100）")
    for name, ns in (("軸 k, k, k", [K, K, K]), ("軸 i, j, k", [I, J, K]), ("軸 k, 45°, k", [K, N45, K])):
        s = "、".join(f"x⃗_{i + 1}3 = {v3(r, 1)}° → x = {x:+.1f}°" for i, r, x in lift((180, 150, 100), ns, 2))
        print(f"  {name}：{s}（弦に沿った Δ_13 = +80°、Δ_23 = +50°）")
    print()
    # ---- 対照：原本 check_values_v2c3.py（原本/ に無変更コピー）を import し、軸が全部 k のとき全ケースで一致することを機械的に確かめる
    spec = importlib.util.spec_from_file_location("orig_v2c3", os.path.join(base, "原本", "check_values_v2c3.py"))
    orig = importlib.util.module_from_spec(spec); spec.loader.exec_module(orig)
    print("== 対照：原本 check_values_v2c3.py との突き合わせ（軸は全部 k。項目ごとの最大差）==")
    worst = 0.0
    for case in cases:
        label, ps = case[0], case[1]; f0 = case[2] if len(case) > 2 else (0, 0, 0)
        d0 = orig.run(ps, f0, steps=[(20, 20), (10, 0), (10, 10)]); d1 = run(ps, f0, steps=[(20, 20), (10, 0), (10, 10)])
        diffs = {}
        diffs['Rs'] = max(abs(a - b) for a, b in zip(d0['Rs'], d1['Rs']))
        for key in ('R', 'ratio', 'remax', 'remin', 'c0R', 'PpiR2', 'mbar'): diffs[key] = abs(d0[key] - d1[key])
        diffs['immax→V_k max'] = abs(d0['immax'] - d1['vmax'][2]); diffs['immin→V_k min'] = abs(d0['immin'] - d1['vmin'][2])
        diffs['Re Σc²'] = abs(d0['s2R2'].real - d1['s2R2'][0]); diffs['Im Σc²→k 成分'] = abs(d0['s2R2'].imag - d1['s2R2'][3])
        diffs['Σc² i, j 成分 = 0'] = max(abs(d1['s2R2'][1]), abs(d1['s2R2'][2]))
        diffs['orient'] = abs(d0['orient'] - d1['orient']); diffs['wind'] = abs(d0['wind'] - d1['wind']); diffs['n'] = abs(d0['n'] - d1['n'])
        pd = 0.0
        for q0, q1 in zip(d0['pairs'], d1['pairs']):   # (i, j, Δ, cross, cosΘ, Φ)
            assert q0[0] == q1['i'] and q0[1] == q1['j']
            pd = max(pd, abs(q0[2] - q1['d']), abs(q0[3] - q1['cross']), abs(q0[5] - q1['f']))
            if not (math.isnan(q0[4]) and math.isnan(q1['cosT'])): pd = max(pd, abs(q0[4] - q1['cosT']))
            pd = max(pd, abs(q1['Dabs'] - abs(q1['d'])), max(abs(q1['Dvec'][0]), abs(q1['Dvec'][1])), abs(q1['Dvec'][2] - q1['d']) if abs(q1['d']) < 180 else 0.0)   # Δ⃗ = Δ k
        diffs['pairs (Δ, cross, cosΘ, Φ, Δ⃗ = Δk)'] = pd
        assert d0['vel'] == d1['vel']
        m = max(diffs.values()); worst = max(worst, m)
        print(f"  {label}：最大差 {m:.1e}（" + "、".join(f"{k} {v:.0e}" for k, v in diffs.items() if v > 0) + ("）" if any(v > 0 for v in diffs.values()) else "すべて 0）"))
    print(f"  全ケースの最大差 {worst:.1e} → {'一致' if worst < 1e-9 else '不一致'}")
    assert worst < 1e-9
