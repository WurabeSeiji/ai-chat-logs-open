#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス10 v2（読出し専用）：非単色状態の倍音・共鳴判定。
v1 の欠陥を修正：(1) 再帰判定は包絡 |W|² の自己相関（うなり周期）に変更（v1 は |·| が位相を殺し T=1 で自明に通過）
(2) 線振動数は整合フィルタ（DTFT 振幅の局所最大化）で精密化し、誤差は前半窓/後半窓の推定差で実測（チューニングなし）
(3) 主判定は共鳴条件（絶対精度で検定可能）：|ω_i+ω_j−2ω_k|、|ω_i+ω_j−ω_k−ω_l| < 誤差和 → パラメトリック共鳴の帳簿
   比 ω_j/ω_base の有理近似は誤差限界つきで併記（限界が粗い場合は判定不能と明記）。"""
import os, csv, json, math
from fractions import Fraction
import numpy as np
from scipy.optimize import minimize_scalar
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
LW = 16384; PAD = 8

def peak_seeds(W, nmax=8):
    n = len(W); h = np.hanning(n)
    F = np.fft.fft(W*h, n*PAD); P = np.abs(F)**2
    fr = np.fft.fftfreq(n*PAD)*2*math.pi
    pk = [(P[i], fr[i]) for i in range(n*PAD)
          if P[i] >= P[(i-1) % (n*PAD)] and P[i] >= P[(i+1) % (n*PAD)]]
    pk.sort(reverse=True)
    pmax = pk[0][0]; out = []
    for p, om in pk:
        if p < 0.01*pmax: break
        if all(abs(om-o) > 4*2*math.pi/n for o, _ in out):
            out.append((om, p/pmax))
        if len(out) >= nmax: break
    return out

def refine(W, om0):
    n = len(W); t = np.arange(n); h = np.hanning(n)
    bin_w = 2*math.pi/n
    f = lambda om: -abs(np.sum(W*h*np.exp(-1j*om*t)))
    r = minimize_scalar(f, bounds=(om0-bin_w, om0+bin_w), method="bounded",
                        options=dict(xatol=1e-9))
    return float(r.x)

def analyze(seg):
    seeds = peak_seeds(seg)
    lines = []
    half = len(seg)//2
    for om0, p in seeds:
        om = refine(seg, om0)
        om_a = refine(seg[:half], om0); om_b = refine(seg[half:], om0)
        err = abs(om_a-om_b) + 1e-9
        lines.append(dict(omega=om, power=p, err=err))
    # 包絡のうなり周期
    env = np.abs(seg)**2; env = env-env.mean()
    n = len(env)
    ac = np.correlate(env, env, "full")[n-1:]
    ac = ac/max(ac[0], 1e-300)
    Tbeat = None
    for T in range(2, n//2):
        if ac[T] >= ac[T-1] and ac[T] >= ac[T+1] and ac[T] > 0.5:
            Tbeat = (T, float(ac[T])); break
    # 共鳴条件（3 波・4 波）
    res3 = []; res4 = []
    L = len(lines)
    for i in range(L):
        for j in range(i, L):
            for k in range(L):
                if k in (i, j): continue
                d = lines[i]["omega"]+lines[j]["omega"]-2*lines[k]["omega"]
                e = lines[i]["err"]+lines[j]["err"]+2*lines[k]["err"]
                if abs(d) < e:
                    res3.append((i, j, k, d, e))
    for i in range(L):
        for j in range(i+1, L):
            for k in range(L):
                for l in range(k+1, L):
                    if len({i, j, k, l}) < 4: continue
                    if (k, l) <= (i, j): continue
                    d = lines[i]["omega"]+lines[j]["omega"]-lines[k]["omega"]-lines[l]["omega"]
                    e = sum(lines[x]["err"] for x in (i, j, k, l))
                    if abs(d) < e:
                        res4.append((i, j, k, l, d, e))
    return lines, Tbeat, res3, res4

targets = []
for r in csv.DictReader(open(os.path.join(ROOT, "results", "composite_wave_summary.csv"))):
    if r["n_lines_late"] not in ("", "None") and int(r["n_lines_late"]) >= 2:
        targets.append((r["tag"], "late", None))
md = ["# 非単色状態の倍音・共鳴判定（パス10 v2）", "",
      "線振動数は整合フィルタで精密化、誤差=前半/後半窓の推定差（実測）。主判定=共鳴条件（絶対精度）。", ""]
rows = []
for tag, win, _ in targets:
    Z = np.load(os.path.join(DATA, tag, "states_treatment.npz"))["Z"]
    W = Z.sum(axis=1); seg = W[-LW:]
    lines, Tbeat, res3, res4 = analyze(seg)
    base = max(lines, key=lambda d: d["power"])["omega"]
    md.append(f"## {tag}（終窓 {LW} step）")
    md.append("")
    md.append("| 線 | ω [rad/step] | ±err | パワー | ω/ω_base | 有理候補 (±限界) |")
    md.append("|---|---|---|---|---|---|")
    for i, d in enumerate(lines):
        r = d["omega"]/base
        err_r = d["err"]/abs(base) + abs(d["omega"])*lines[0]["err"]/base**2
        fr = Fraction(r).limit_denominator(24)
        verdict = f"{fr.numerator}/{fr.denominator}" if abs(r-float(fr)) < err_r and err_r < 0.02 else ("判定不能" if err_r >= 0.02 else "非通約")
        md.append(f"| {i} | {d['omega']:+.7f} | {d['err']:.1e} | {d['power']:.3f} | {r:+.5f} | {verdict} (±{err_r:.4f}) |")
    md.append("")
    if res3:
        md.append("**3 波共鳴 ω_i+ω_j=2ω_k:** " + "; ".join(f"({i},{j};{k}) 残差{d:+.1e}<±{e:.1e}" for i, j, k, d, e in res3))
    if res4:
        md.append("**4 波共鳴 ω_i+ω_j=ω_k+ω_l:** " + "; ".join(f"({i},{j};{k},{l}) 残差{d:+.1e}<±{e:.1e}" for i, j, k, l, d, e in res4))
    md.append(f"**うなり周期（包絡自己相関 >0.5 の最初のピーク）:** {Tbeat if Tbeat else 'なし'}")
    md.append("")
    rows.append(dict(tag=tag, lines=";".join(f"{d['omega']:+.7f}±{d['err']:.1e}:{d['power']:.3f}" for d in lines),
                     n_res3=len(res3), n_res4=len(res4), beat=Tbeat[0] if Tbeat else None))
    print(f"{tag}: " + " ".join(f"{d['omega']:+.6f}±{d['err']:.0e}" for d in lines))
    for i, j, k, d, e in res3: print(f"   3波 ({i},{j};{k}): {d:+.2e} < ±{e:.1e}")
    for i, j, k, l, d, e in res4: print(f"   4波 ({i},{j};{k},{l}): {d:+.2e} < ±{e:.1e}")
    if Tbeat: print(f"   うなり T={Tbeat[0]} (C={Tbeat[1]:.3f})")
with open(os.path.join(ROOT, "results", "harmonic_ladder.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
with open(os.path.join(ROOT, "results", "harmonic_ladder.md"), "w") as f:
    f.write("\n".join(md)+"\n")
print("PASS10 OK")
