#!/usr/bin/env python3
"""results/lock_*.csv と対照走行を読んで事実だけを表にする（合否は書かない）。"""
import csv, glob, os, re

def load(p):
    with open(p) as fh:
        return list(csv.DictReader(fh))

def fl(r, k):
    return float(r[k])

def nearest_kN(frac, N):
    k = round(frac * N)
    return k, k / N, frac - k / N

print("# 関係配列写像 v2 ― 引き込み二パターンの比較（事実のみ）\n")
if os.path.exists(os.path.join("results", "control_compare.txt")):
    print("## 対照走行（--lock none を v1 の csv と照合）\n")
    print("```")
    print(open(os.path.join("results", "control_compare.txt")).read().strip())
    print("```\n")

print("## 各走行の初期と最終\n")
print("| 走行 | N | 節点 | 2θ/2π 初期 | 2θ/2π 最終 | 最寄りの k/N（最終） | 差 | \\|H\\|² 最終 | 損失率 最終 | m 最終 | x₃ᵃ 最終 | ξ 最終 | 2θ/2π の変化が 10⁻¹²/節点 を切った節点 |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for p in sorted(glob.glob(os.path.join("results", "lock_*.csv"))):
    rows = load(p)
    mobj = re.search(r"lock_([AB])_N(\d+)", p)
    L, N = mobj.group(1), int(mobj.group(2))
    r0, rN = rows[0], rows[-1]
    f0, fN = fl(r0, "two_theta_over_2pi"), fl(rN, "two_theta_over_2pi")
    k, kN, diff = nearest_kN(fN, N)
    # 変化率が閾値を切った最初の節点（1000 節点ごとの差分で評価）
    stop = ""
    step = 1000
    for i in range(step, len(rows), step):
        d = abs(fl(rows[i], "two_theta_over_2pi") - fl(rows[i - step], "two_theta_over_2pi")) / step
        if d < 1e-12:
            stop = str(int(rows[i]["node"]))
            break
    print(f"| {L} | {N} | {int(rN['node'])} | {f0:.9f} | {fN:.9f} | {k}/{N} = {kN:.9f} | {diff:+.2e} | "
          f"{fl(rN,'H2'):.3e} | {fl(rN,'loss_frac'):.3e} | {fl(rN,'m'):.9f} | {fl(rN,'x3a'):.4f} | {fl(rN,'xi'):.6f} | {stop or '切らない'} |")

print("\n## 途中経過（10000 節点ごとの 2θ/2π と \\|H\\|²）\n")
hdr = None
for p in sorted(glob.glob(os.path.join("results", "lock_*.csv"))):
    rows = load(p)
    mobj = re.search(r"lock_([AB])_N(\d+)", p)
    L, N = mobj.group(1), int(mobj.group(2))
    pts = [r for r in rows if int(r["node"]) % 10000 == 0]
    if hdr is None:
        hdr = [int(r["node"]) for r in pts]
        print("| 走行 | " + " | ".join(f"n={h}" for h in hdr) + " |")
        print("|---|" + "---|" * len(hdr))
    print(f"| {L} N={N}: 2θ/2π | " + " | ".join(f"{fl(r,'two_theta_over_2pi'):.8f}" for r in pts) + " |")
    print(f"| {L} N={N}: \\|H\\|² | " + " | ".join(f"{fl(r,'H2'):.2e}" for r in pts) + " |")
