#!/usr/bin/env python3
"""results/*.csv を読んで事実だけを表にする（合否は書かない）。"""
import csv, glob, math, os

def load(p):
    with open(p) as fh:
        return list(csv.DictReader(fh))

def fl(r, k):
    return float(r[k])

print("# 半角位相子の交換（固定点 4）― 走行の要約（事実のみ）\n")
for name in ["kepler_check.txt"]:
    pth = os.path.join("results", name)
    if os.path.exists(pth):
        print("## 対照：Kepler の閉じた式との比較（--em kepler --f 0）\n\n```")
        print(open(pth).read().strip())
        print("```\n")

for p in sorted(glob.glob(os.path.join("results", "run_*.csv"))):
    rows = load(p)
    r0, rN = rows[0], rows[-1]
    n = int(rN["node"])
    # 周回数：双極子の位相の進みの積算
    turns = sum(fl(r, "d_argz") for r in rows[:-1]) / (2 * math.pi)
    print(f"## {os.path.basename(p)}（節点 {n}、t = {fl(rN,'t'):.4e}、双極子の周回 {turns:.3f}）\n")
    print("| 量 | 初期 | 最終 | 最小 | 最大 |")
    print("|---|---|---|---|---|")
    for lab, key, fmt in [("ξ", "xi", ".6f"), ("a（長半径 = N_A+N_B）", "a_orb", ".9f"), ("e", "ecc", ".6e"),
                          ("L", "L", ".9f"), ("ε_K", "eps_K", ".9e"), ("ε_full", "eps_full", ".9e"),
                          ("N_p", "N_p", ".6f"), ("N_q", "N_q", ".6f"), ("Re⟨p|q⟩", "Re_pq", ".6f"), ("Im⟨p|q⟩", "Im_pq", ".6f"),
                          ("θ", "theta", ".9f"), ("R = sin²θ", "R_sin2theta", ".6e"), ("Δt（一節点）", "dt_node", ".4f"),
                          ("x₃ᵇ/x₃ᵃ", "x3b_over_x3a", ".9f"), ("s_dir（読み a）", "s_dir", ".4e"),
                          ("sin²θ_rad", "sin2_theta_rad", ".4e"), ("近点移動率 ∂⟨δε⟩/∂L", "prec_rate", ".4e")]:
        v = [fl(r, key) for r in rows]
        print(f"| {lab} | {v[0]:{fmt}} | {v[-1]:{fmt}} | {min(v):{fmt}} | {max(v):{fmt}} |")
    # 近点移動：一周あたり 2·Σdalpha / 周回数
    dal = sum(fl(r, "dalpha_node") for r in rows[:-1])
    if turns > 0:
        print(f"\n近点移動：共通位相の積算 Σα = {dal:.6e} → 軌道の回転 2Σα = {2*dal:.6e} rad、一周あたり {2*dal/turns:.6e} rad（v7 の 2π(γ_cl−1) = {2*math.pi*2.0e-5:.4e}）")
    # 読み (b)：全配列の初期への戻り
    cfg0 = (fl(r0, "arg_z"), fl(r0, "x1a"), fl(r0, "x1b"))
    rec = []
    for r in rows[1:]:
        d = [abs((fl(r, k) - c0 + math.pi) % (2 * math.pi) - math.pi) for k, c0 in zip(("arg_z", "x1a", "x1b"), cfg0)]
        if max(d) < 1e-3 * 2 * math.pi:
            rec.append((int(r["node"]), fl(r, "t")))
    print(f"読み (b)：双極子の位相と二つの時計が同時に初期へ戻った節点（各 1/1000 周以内）：{len(rec)} 回 {rec[:6]}")
    # 読み (a) の最小
    sd = [fl(r, "s_dir") for r in rows]
    i = sd.index(min(sd))
    print(f"読み (a)：s_dir の最小 {min(sd):.3e}（節点 {rows[i]['node']}）、平均 {sum(sd)/len(sd):.4f}\n")
