#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス12（読出し専用・木原指示）：N=3,4,5,12 の全タグ・全関係波について、
終窓の振動数 ν_e と波長 λ_e=1/|ν_e|（最長でなく最短 λ=1 に規格化）を全数列挙する。分析・解釈なし、列挙のみ。"""
import os, math, csv
import numpy as np
from scipy.optimize import minimize_scalar
import sys
sys.path.insert(0, os.path.dirname(__file__))
from common import edges
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
LW = 16384; PAD = 8

def refine(W, om0):
    n = len(W); t = np.arange(n); h = np.hanning(n)
    bw = 2*math.pi/n
    f = lambda om: -abs(np.sum(W*h*np.exp(-1j*om*t)))
    return float(minimize_scalar(f, bounds=(om0-bw, om0+bw), method="bounded",
                                 options=dict(xatol=1e-9)).x)

def dominant(sig):
    n = len(sig); h = np.hanning(n)
    F = np.fft.fft(sig*h, n*PAD); P = np.abs(F)**2
    fr = np.fft.fftfreq(n*PAD)*2*math.pi
    i = int(np.argmax(P))
    return refine(sig, fr[i])

TAGS = [t for N in (3, 4, 5, 12) for t in (f"mp_N{N}", f"hm_N{N}", f"ne_N{N}", f"rb_N{N}")
        if os.path.exists(os.path.join(DATA, t, "states_treatment.npz"))]
md = ["# 全関係波の振動数と波長の全数列挙（N=3,4,5,12、終窓 16384 step、最短 λ=1 規格化）", ""]
rows = []
for tag in TAGS:
    N = int(tag.split("N")[1]); E = edges(N); M = len(E)
    Z = np.load(os.path.join(DATA, tag, "states_treatment.npz"))["Z"][-LW:]
    amp = np.abs(Z).mean(axis=0); amax = amp.max()
    nus = [dominant(Z[:, e]) for e in range(M)]
    lam = [1.0/abs(v) if abs(v) > 1e-12 else float("inf") for v in nus]
    lmin = min(lam)
    md.append(f"## {tag}")
    md.append("")
    md.append("| 辺 (i,j) | 振幅/max | ν_e [rad/step] | λ_e=1/|ν_e| （最短=1） |")
    md.append("|---|---|---|---|")
    for e in range(M):
        md.append(f"| {E[e]} | {amp[e]/amax:.3f} | {nus[e]:+.7f} | {lam[e]/lmin:.4f} |")
        rows.append(dict(tag=tag, edge=str(E[e]), amp_rel=amp[e]/amax, nu=nus[e], lam_norm=lam[e]/lmin))
    # 波長グループ（0.5% で束ね）
    md.append("")
    groups = []
    for e in range(M):
        v = lam[e]/lmin
        for g in groups:
            if abs(v-g[0])/g[0] < 0.005:
                g[1].append(e); g[0] = (g[0]*(len(g[1])-1)+v)/len(g[1]); break
        else:
            groups.append([v, [e]])
    groups.sort()
    md.append("**波長グループ:** " + " / ".join(f"λ={g[0]:.4f} × {len(g[1])} 本" for g in groups))
    md.append("")
    print(f"{tag}: " + " / ".join(f"λ={g[0]:.4f}×{len(g[1])}本" for g in groups))
with open(os.path.join(ROOT, "results", "wavelength_table.md"), "w") as f:
    f.write("\n".join(md)+"\n")
with open(os.path.join(ROOT, "results", "wavelength_table.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("PASS12 OK")
