#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""|Q|² の ψ 平均と R² の関係（2026-11-07）。係数だけの計算と ψ の数値平均で確かめる。
  Q(ψ) = Σ_i e^{n_i φ₀i} [ c₀ᵢ + Σ_m c⁽ⁱ⁾_m e^{n_i m(ψ − θ_i)} ]、各項は cos(mψ + α_i) + n_i sin(mψ + α_i)、α_i = φ₀i − mθ_i
  (1) ψ 平均の交差項（同じ m、i ≠ j）：⟨Sc[q_i q̄_j]⟩_ψ = ½ (1 + n_i·n_j) cos(α_i − α_j)
  (2) R²（index_v2q3 の定義、係数の四元数ノルム）の交差項：Sc[e^{n_i α_i} e^{−n_j α_j}] = cos α_i cos α_j + (n_i·n_j) sin α_i sin α_j
  n_i = n_j のとき (1) = (2) = cos(α_i − α_j)。違うとき (2) は ψ の原点 s（θ_i → θ_i − s、α_i → α_i + ms）に依り、
  (1) は依らない。(2) を s で平均すると (1) になる：⟨R²⟩_s = ⟨|Q|²⟩_ψ − |C₀|²
  直流：C₀ = Σ_i c₀ᵢ e^{n_i φ₀ᵢ} ∈ ℍ、|C₀|² は ψ にも s にも依らない定数項。倍音との交差は平均で 0"""
import math, sys

def coeffs(phi, M=10):
    ms = list(range(1, M + 1)); cs = [(2 / math.pi) * math.sin(m * math.pi * phi / 360) / m for m in ms]
    return ms, cs, phi / 360
def unit(v):
    L = math.sqrt(sum(x * x for x in v)); return [0, 0, 1] if L < 1e-12 else [x / L for x in v]
def qexp(n, a): return [math.cos(a), n[0] * math.sin(a), n[1] * math.sin(a), n[2] * math.sin(a)]
def qmul(p, q): return [p[0]*q[0]-p[1]*q[1]-p[2]*q[2]-p[3]*q[3], p[0]*q[1]+p[1]*q[0]+p[2]*q[3]-p[3]*q[2], p[0]*q[2]-p[1]*q[3]+p[2]*q[0]+p[3]*q[1], p[0]*q[3]+p[1]*q[2]-p[2]*q[1]+p[3]*q[0]]
def qconj(q): return [q[0], -q[1], -q[2], -q[3]]
def dot(a, b): return sum(x * y for x, y in zip(a, b))

def prep(pulses, shift=0.0):
    P = []
    for u in pulses:
        ms, cs, c0 = coeffs(u["phi"]); P.append(dict(u, ms=ms, cs=cs, c0=c0, n=unit(u["n"]), f=math.radians(u["phase0"]), t=math.radians(u["theta"] - shift)))
    return P

def quantities(P, npts=1441):
    """返り値：|C₀|²、R²（係数ノルム）、閉じた式の ψ 平均（直流を除く）、数値の ψ 平均（直流を含む）"""
    C0 = [0, 0, 0, 0]
    for u in P:
        e = qexp(u["n"], u["f"]); C0 = [a + u["c0"] * b for a, b in zip(C0, e)]
    c0abs2 = sum(x * x for x in C0)
    R2 = 0.0; mc = 0.0
    for u in P: R2 += sum(c * c for c in u["cs"]); mc += sum(c * c for c in u["cs"])
    for i in range(len(P)):
        for j in range(i + 1, len(P)):
            a, b = P[i], P[j]; nd = dot(a["n"], b["n"])
            for m, ca in zip(a["ms"], a["cs"]):
                cb = b["cs"][b["ms"].index(m)]; aa, ab = a["f"] - m * a["t"], b["f"] - m * b["t"]
                g = qmul(qexp(a["n"], aa), qconj(qexp(b["n"], ab)))
                R2 += 2 * ca * cb * g[0]; mc += 2 * ca * cb * 0.5 * (1 + nd) * math.cos(aa - ab)
    acc = 0.0
    for k in range(npts - 1):
        psi = 2 * math.pi * k / (npts - 1); Q = [0, 0, 0, 0]
        for u in P:
            re = u["c0"]; im = 0.0
            for m, c in zip(u["ms"], u["cs"]):
                x = m * (psi - u["t"]); re += c * math.cos(x); im += c * math.sin(x)
            w = re * math.cos(u["f"]) - im * math.sin(u["f"]); v = re * math.sin(u["f"]) + im * math.cos(u["f"])
            Q[0] += w; Q[1] += v * u["n"][0]; Q[2] += v * u["n"][1]; Q[3] += v * u["n"][2]
        acc += sum(x * x for x in Q)
    return c0abs2, R2, mc, acc / (npts - 1)

I, J, K = [1, 0, 0], [0, 1, 0], [0, 0, 1]
CASES = {
  "default": [dict(theta=180, phi=20, phase0=0, n=K), dict(theta=180, phi=10, phase0=0, n=K), dict(theta=180, phi=10, phase0=0, n=K)],
  "tri":     [dict(theta=180, phi=10, phase0=0, n=K), dict(theta=150, phi=10, phase0=0, n=K), dict(theta=100, phi=10, phase0=0, n=K)],
  "mirror":  [dict(theta=180, phi=10, phase0=0, n=K), dict(theta=200, phi=10, phase0=45, n=K), dict(theta=160, phi=10, phase0=-45, n=K)],
  "std":     [dict(theta=180, phi=20, phase0=0, n=K), dict(theta=180, phi=10, phase0=0, n=K), dict(theta=90, phi=10, phase0=0, n=K)],
  "tri_ijk": [dict(theta=180, phi=10, phase0=0, n=I), dict(theta=150, phi=10, phase0=0, n=J), dict(theta=100, phi=10, phase0=0, n=K)],
  "mirror_k_i_mi": [dict(theta=180, phi=10, phase0=0, n=K), dict(theta=200, phi=10, phase0=45, n=I), dict(theta=160, phi=10, phase0=-45, n=[-1, 0, 0])],
  "same_pos_ijk": [dict(theta=150, phi=10, phase0=0, n=I), dict(theta=150, phi=10, phase0=0, n=J), dict(theta=150, phi=10, phase0=0, n=K)],
}
def closed_mean(name):
    c0abs2, R2, mc, mn = quantities(prep(CASES[name])); return c0abs2 + mc

def main():
    bad = 0
    print("A. 各配置：|C₀|² + R²、閉じた式の ψ 平均、数値の ψ 平均")
    for name, ps in CASES.items():
        P = prep(ps); same = all(dot(u["n"], P[0]["n"]) > 1 - 1e-9 for u in P)
        c0abs2, R2, mc, mn = quantities(P)
        ok1 = abs(mn - (c0abs2 + mc)) < 1e-9; ok2 = abs(mn - (c0abs2 + R2)) < 1e-9
        print(f"  {name:14s} 軸{'同' if same else '異'}  |C₀|²+R² = {c0abs2 + R2:.6f}  閉じた式 = {c0abs2 + mc:.6f}  数値 = {mn:.6f}  数値=閉じた式 {'一致' if ok1 else '不一致'}  数値=|C₀|²+R² {'一致' if ok2 else '不一致'}")
        bad += (not ok1) + ((not ok2) if same else (ok2 and name != 'default'))
    print("B. ψ の原点をずらす（θ_i → θ_i − s）：R² は軸が違うと s に依り、ψ 平均は依らない。R² の s 平均（0〜360° を 360 点）は ψ 平均に一致する")
    for name in ["tri", "tri_ijk", "mirror_k_i_mi"]:
        row = []; acc = 0.0
        for s in [0, 30, 60, 90]:
            c0abs2, R2, mc, mn = quantities(prep(CASES[name], s), npts=361); row.append(f"s={s:3d}: R²={R2:.6f} ⟨|Q|²⟩={mn:.6f}")
        for s in range(360):
            c0abs2, R2, mc, mn = quantities(prep(CASES[name], s), npts=2); acc += R2
        c0abs2, R2, mc, mn = quantities(prep(CASES[name]))
        ok = abs(acc / 360 - (mn - c0abs2)) < 1e-9
        print(f"  {name:14s} " + "  ".join(row) + f"   ⟨R²⟩_s = {acc / 360:.6f} ＝ ⟨|Q|²⟩ − |C₀|² = {mn - c0abs2:.6f} {'一致' if ok else '不一致'}")
        bad += (not ok)
    print("結論：ψ 平均 ⟨|Q|²⟩ は閉じた式（交差項 ½(1+n_i·n_j) cos(α_i−α_j)）と全例で一致し、ψ の原点に依らない。R²（係数ノルム）は軸が同じなら同じ値だが、軸が違うと ψ の原点 s に依り、その s 平均が ψ 平均になる。")
    return 0 if bad == 0 else 2

if __name__ == "__main__":
    sys.exit(main())
