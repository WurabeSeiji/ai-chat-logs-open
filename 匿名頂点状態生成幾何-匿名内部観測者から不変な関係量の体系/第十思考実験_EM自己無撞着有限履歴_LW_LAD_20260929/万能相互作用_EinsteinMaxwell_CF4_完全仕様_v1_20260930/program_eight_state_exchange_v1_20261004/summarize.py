#!/usr/bin/env python3
"""results/*.csv(.gz) を読んで事実だけを表にする（合否は書かない）。"""
import csv, glob, gzip, math, os

def load(p):
    op = gzip.open if p.endswith(".gz") else open
    with op(p, "rt") as fh:
        return list(csv.DictReader(fh))

def fl(r, k):
    return float(r[k])

print("# 八状態の指数表現・固定交換係数 θ₀ = 23π/124（固定点 5）― 走行の要約（事実のみ）\n")
pth = os.path.join("results", "kepler_check.txt")
if os.path.exists(pth):
    print("## 対照：Kepler の閉じた式との比較（--em kepler --f 0）\n\n```")
    print(open(pth).read().strip()); print("```\n")
pth = os.path.join("results", "precession_check.txt")
if os.path.exists(pth):
    print("## 対照：近点移動（E–M の解析値 ∂⟨δε⟩/∂L との比較）\n\n```")
    print(open(pth).read().strip()); print("```\n")

files = sorted(glob.glob(os.path.join("results", "run_*.csv"))) or sorted(glob.glob(os.path.join("results", "run_*.csv.gz")))
for p in files:
    rows = load(p)
    r0, rN = rows[0], rows[-1]
    n = int(rN["node"])
    turns = n * 23.0 / 124.0
    print(f"## {os.path.basename(p)}（節点 {n} ＝ {turns:.1f} 周、T = {fl(rN,'T'):.4e}）\n")
    print("| 量 | 初期 | 最終 | 最小 | 最大 |")
    print("|---|---|---|---|---|")
    for lab, key, fmt in [("ξ", "xi", ".6f"), ("a（長半径）", "a_orb", ".9f"), ("e", "ecc", ".6e"), ("L", "L", ".9f"),
                          ("ε_K", "eps_K", ".9e"), ("ε_full", "eps_full", ".9e"), ("ω", "omega", ".9e"), ("δτ", "dtau", ".6f"),
                          ("|G|", "absG", ".9f"), ("arg G", "argG", ".9f"), ("ΔT（一節点）", "dT_node", ".3f"),
                          ("x₃ᵃ（a の時計の一節点の進み）", "x3a", ".3f"), ("x₃ᵇ", "x3b", ".1f"),
                          ("s_dir（読み a）", "s_dir", ".4e"), ("sin²θ_rad", "sin2_theta_rad", ".4e")]:
        v = [fl(r, key) for r in rows]
        print(f"| {lab} | {v[0]:{fmt}} | {v[-1]:{fmt}} | {min(v):{fmt}} | {max(v):{fmt}} |")
    # 近点移動：離心率ベクトルの向き varpi の積算
    dv = 0.0; prev = fl(r0, "varpi")
    for r in rows[1:]:
        cur = fl(r, "varpi"); dv += (cur - prev + math.pi) % (2 * math.pi) - math.pi; prev = cur
    print(f"\n近点方向 ϖ の変化（積算）：{dv:.6e} rad、一周あたり {dv/turns:.6e} rad")
    # 読み (b)：124 節点ごとの戻り
    rec = [(int(r["node"]), fl(r, "rec124_dz"), fl(r, "rec124_dclock_a"), fl(r, "rec124_dclock_b")) for r in rows if r["rec124_dz"] not in ("nan", "")]
    if rec:
        print(f"読み (b)：124 節点ごとの戻り {len(rec)} 回。関係 z の相対差の最大 {max(d for _, d, _, _ in rec):.2e}。"
              f"時計 a の位相差 [rad] の範囲 {min(c for _, _, c, _ in rec):+.3f}..{max(c for _, _, c, _ in rec):+.3f}、"
              f"時計 b の範囲 {min(c for _, _, _, c in rec):+.3f}..{max(c for _, _, _, c in rec):+.3f}。"
              f"両時計が同時に 0.1 rad 以内に戻った回数 {sum(1 for _, _, ca, cb in rec if abs(ca) < 0.1 and abs(cb) < 0.1)}")
    sd = [fl(r, "s_dir") for r in rows]
    print(f"読み (a)：s_dir の最小 {min(sd):.3e}、平均 {sum(sd)/len(sd):.4f}\n")
