#!/usr/bin/env python3
"""results/*.csv(.gz) を読んで事実だけを表にする（合否は書かない）。v2：状態から読んだ時間（ΔT、E_cycle、t_cycle、時計の文字盤）。"""
import csv, glob, gzip, math, os

def load(p):
    op = gzip.open if p.endswith(".gz") else open
    with op(p, "rt") as fh:
        return list(csv.DictReader(fh))

def fl(r, k):
    return float(r[k])

print("# 八状態の指数表現・固定交換係数 θ₀ = 23π/124 v2 ― 走行の要約（事実のみ）\n")
pth = os.path.join("results", "kepler_check.txt")
if os.path.exists(pth):
    print("## 対照：Kepler の閉じた式との比較（--em kepler --f 0）\n\n```")
    print(open(pth).read().strip()); print("```\n")

files = sorted(glob.glob(os.path.join("results", "run_*.csv"))) or sorted(glob.glob(os.path.join("results", "run_*.csv.gz")))
for p in files:
    rows = load(p)
    r0, rN = rows[0], rows[-1]
    n = int(rN["node"])
    turns = n * 23.0 / 124.0
    print(f"## {os.path.basename(p)}（節点 {n} ＝ {turns:.1f} 周）\n")
    print("| 量 | 初期 | 最終 | 最小 | 最大 |")
    print("|---|---|---|---|---|")
    for lab, key, fmt in [("ξ", "xi", ".6f"), ("a（長半径）", "a_orb", ".9f"), ("e", "ecc", ".6e"), ("L", "L", ".9f"),
                          ("ε_K", "eps_K", ".9e"), ("ω", "omega", ".9e"), ("δτ", "dtau", ".6f"), ("ΔT（一節点、状態から）", "dT_node", ".3f"),
                          ("1PN 係数 corr", "corr_1PN", ".9f"),
                          ("a₀（ln 振幅 A）", "a0", ".9f"), ("a₂（ln\|g\|）", "a2", ".9e"), ("a₃（A の時計の進み）", "a3", ".3f"),
                          ("b₀", "b0", ".9f"), ("b₃（B の時計の進み）", "b3", ".1f"),
                          ("s_dir（三位相の差の読み）", "s_dir", ".4e"), ("sin²θ_rad", "sin2_theta_rad", ".4e")]:
        v = [fl(r, key) for r in rows]
        print(f"| {lab} | {v[0]:{fmt}} | {v[-1]:{fmt}} | {min(v):{fmt}} | {max(v):{fmt}} |")
    # 時計の文字盤（状態）が周期内の位置 E と同時に初期配置へ戻る節点（状態の読出しの比較。帳簿なし）
    E0, a10, b10 = fl(r0, "E_cycle"), fl(r0, "a1"), fl(r0, "b1")
    hits = []
    for r in rows[1:]:
        d = [abs((fl(r, k) - c0 + math.pi) % (2 * math.pi) - math.pi) for k, c0 in (("E_cycle", E0), ("a1", a10), ("b1", b10))]
        if max(d) < 0.1:
            hits.append(int(r["node"]))
    print(f"\n周期内の位置 E と二つの時計の文字盤 a₁、b₁ が同時に初期配置へ 0.1 rad 以内で戻った節点：{len(hits)} 回 {hits[:8]}")
    sd = [fl(r, "s_dir") for r in rows]
    print(f"読み (a)：s_dir の最小 {min(sd):.3e}、平均 {sum(sd)/len(sd):.4f}\n")
