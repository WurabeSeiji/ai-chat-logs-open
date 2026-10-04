#!/usr/bin/env python3
"""results/*.csv を読んで事実だけを表にする（合否は書かない）。v3：三位相の差 Δ、源の因子 s_dir、ロック到達。"""
import csv, glob, math, os

def load(p):
    with open(p) as fh:
        return list(csv.DictReader(fh))

def fl(r, k):
    return float(r[k])

print("# 関係配列写像 v3 ― 決定済みの引き込みの走行（事実のみ）\n")
cc = os.path.join("results", "control_compare.txt")
if os.path.exists(cc):
    print("## 対照走行（--lock off を v1 の csv と照合）\n\n```")
    print(open(cc).read().strip())
    print("```\n")

for p in sorted(glob.glob(os.path.join("results", "lock_*.csv"))):
    rows = load(p)
    r0, rN = rows[0], rows[-1]
    n = int(rN["node"])
    print(f"## {os.path.basename(p)}（節点 {n}、一周数 {fl(rN,'cum_turns'):.3f}、K = {r0['K']}）\n")
    print("| 量 | 初期 | 最終 |")
    print("|---|---|---|")
    for name, key, fmt in [("ξ", "xi", ".9f"), ("m", "m", ".12f"), ("θ [rad]", "theta", ".9e"),
                           ("2θ/2π", "two_theta_over_2pi", ".9f"), ("x₃ᵃ", "x3a", ".6f"),
                           ("x₃ᵇ/x₃ᵃ", "x3b_over_x3a", ".9f"),
                           ("Δa = 2θ − x₃ᵃ (mod 2π)", "Delta_a", ".3e"), ("Δb = 2θ − x₃ᵇ (mod 2π)", "Delta_b", ".6e"),
                           ("Δab = x₃ᵇ − x₃ᵃ (mod 2π)", "Delta_ab", ".6e"), ("s_dir（源の因子）", "s_dir", ".6e"),
                           ("sin²θ_rad（実効）", "sin2_theta_rad", ".4e"), ("ε", "eps", ".12e")]:
        print(f"| {name} | {fl(r0,key):{fmt}} | {fl(rN,key):{fmt}} |")
    # ロック到達：s_dir が閾値を初めて切った節点と、その後の最大値
    for thr in (1e-2, 1e-4, 1e-6, 1e-8):
        hit = next((r for r in rows if fl(r, "s_dir") < thr), None)
        if hit:
            after = rows[int(hit["node"]):]
            mx = max(fl(r, "s_dir") for r in after)
            print(f"\n- s_dir < {thr:.0e} に初めて達した節点 n = {hit['node']}（{fl(hit,'cum_turns'):.2f} 周）、"
                  f"そのとき ξ = {fl(hit,'xi'):.6f}、m = {fl(hit,'m'):.9f}、2θ/2π = {fl(hit,'two_theta_over_2pi'):.9f}、"
                  f"Δb = {fl(hit,'Delta_b'):.3e}。以後の s_dir の最大 {mx:.3e}")
        else:
            print(f"\n- s_dir < {thr:.0e} には達していない（最小 {min(fl(r,'s_dir') for r in rows):.3e}）")
    # 10 分割の経過
    step = max(1, n // 10)
    pts = [r for r in rows if int(r["node"]) % step == 0]
    print("\n| n | ξ | m | 2θ/2π | Δb | s_dir | sin²θ_rad |")
    print("|---|---|---|---|---|---|---|")
    for r in pts:
        print(f"| {r['node']} | {fl(r,'xi'):.6f} | {fl(r,'m'):.9f} | {fl(r,'two_theta_over_2pi'):.9f} | "
              f"{fl(r,'Delta_b'):+.4e} | {fl(r,'s_dir'):.3e} | {fl(r,'sin2_theta_rad'):.3e} |")
    print()
