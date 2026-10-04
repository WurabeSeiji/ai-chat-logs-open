#!/usr/bin/env python3
"""results/*.csv を読んで事実だけを表にする（合否は書かない）。"""
import csv, math, glob, os

def load(p):
    with open(p) as fh:
        return list(csv.DictReader(fh))

def fl(r, k):
    return float(r[k])

print("# 関係配列写像 v1 ― 走行の要約（事実のみ）\n")
for p in sorted(glob.glob(os.path.join("results", "run_*.csv"))):
    rows = load(p)
    r0, rN = rows[0], rows[-1]
    n = int(rN["node"])
    xis = [fl(r, "xi") for r in rows]
    ths = [fl(r, "theta") for r in rows]
    ms = [fl(r, "m") for r in rows]
    es = [fl(r, "eps") for r in rows]
    gam = [fl(r, "gamma") for r in rows]
    turns = fl(rN, "cum_turns")
    print(f"## {os.path.basename(p)}  （節点 {n}、一周数 {turns:.4f}）\n")
    print("| 量 | 初期 | 最終 | 最小 | 最大 |")
    print("|---|---|---|---|---|")
    for name, arr, fmt in [("ξ", xis, ".9f"), ("θ [rad]", ths, ".9e"), ("m", ms, ".12f"), ("ε", es, ".12e"),
                           ("γ = θ/θ_X", gam, ".9f")]:
        print(f"| {name} | {arr[0]:{fmt}} | {arr[-1]:{fmt}} | {min(arr):{fmt}} | {max(arr):{fmt}} |")
    print(f"| R = sin²θ | {fl(r0,'R_sin2theta'):.6e} | {fl(rN,'R_sin2theta'):.6e} | | |")
    print(f"| 2θ/2π（一節点の周の割合） | {fl(r0,'two_theta_over_2pi'):.9f} | {fl(rN,'two_theta_over_2pi'):.9f} | | |")
    print(f"| x₃ᵇ/x₃ᵃ | {fl(r0,'x3b_over_x3a'):.9f} | {fl(rN,'x3b_over_x3a'):.9f} | | |")
    print(f"| F_a, F_b（初期） | {fl(r0,'F_a'):.12f}, {fl(r0,'F_b'):.12f} | | | |")
    print(f"| sin²θ_rad（初期／最終） | {fl(r0,'sin2_theta_rad'):.4e} | {fl(rN,'sin2_theta_rad'):.4e} | | |")
    # 配列の位相 ψ が初期に戻る節点（一周の 1/1000 以内）
    rec = []
    for r in rows[1:]:
        psi = fl(r, "psi")
        d = min(psi, 2 * math.pi - psi)
        if d < 1e-3 * 2 * math.pi:
            rec.append((int(r["node"]), fl(r, "cum_turns"), d))
    print(f"\n配列の位相 ψ が初期値に戻った節点（一周の 1/1000 以内）: {len(rec)} 回。先頭: "
          + ", ".join(f"n={a}（{b:.3f} 周, |Δψ|={c:.2e}）" for a, b, c in rec[:8]) + "\n")
    # 節点あたりの周の割合の有理近似
    frac = fl(r0, "two_theta_over_2pi")
    best = min(((abs(frac - q / p_), p_, q) for p_ in range(1, 200) for q in range(1, p_)), key=lambda t: t[0])
    print(f"初期の 2θ/2π = {frac:.9f}。分母 200 以下の最良有理近似 {best[2]}/{best[1]}（差 {best[0]:.2e}）\n")
