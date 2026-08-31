#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス11（読出し専用・実数シンプレックス検討の本体）：辺ごとの波長読み。
鎖（目的の正本 §8）：固有時計 ν0 → 辺 e の振動数 ν_e = m_e·ν0（辺の時系列 z_e(t) のスペクトルで実測）
→ 波長 λ_e = λ0/|m_e| → 長さ L_e = k_e·λ_e/2 = (k_e/|m_e|)·(λ0/2)。測れるのは m_e、未知は k_e。
手順（調整なし）：終窓 16384 step の z_e(t) を整合フィルタで読む。単色度=支配線パワー/全線パワー。
振幅が丸め水準（max辺の 1e-6 未満）の辺は「波なし＝長さ未定義」。
辺振動数の共通測度：全単色辺の ν_e の最小値を仮基音 ν0 とし、r_e=ν_e/ν0 の連分数近似（誤差=分割窓の実測差を伝播）。
実単体テスト（第一歩）：k_e≡1 と置いた L_e ∝ 1/|ν_e| で距離行列を作り、三角不等式と Cayley–Menger（B の PSD）を機械判定。
通れば k≡1 で実単体が成立、通らなければ {k_e} の数え上げが次段。"""
import os, csv, json, math
from fractions import Fraction
import numpy as np
from scipy.optimize import minimize_scalar
import sys
sys.path.insert(0, os.path.dirname(__file__))
from common import edges
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
LW = 16384; PAD = 8
TARGETS = ["hm_N12", "mp_N5", "rb_N5", "mp_N3", "ne_N3", "ne_N4"]

def refine(W, om0):
    n = len(W); t = np.arange(n); h = np.hanning(n)
    bw = 2*math.pi/n
    f = lambda om: -abs(np.sum(W*h*np.exp(-1j*om*t)))
    return float(minimize_scalar(f, bounds=(om0-bw, om0+bw), method="bounded",
                                 options=dict(xatol=1e-9)).x)

def edge_line(sig):
    """辺信号の支配線 (ω, err, 単色度)。"""
    n = len(sig); h = np.hanning(n)
    F = np.fft.fft(sig*h, n*PAD); P = np.abs(F)**2
    fr = np.fft.fftfreq(n*PAD)*2*math.pi
    pk = [(P[i], fr[i]) for i in range(n*PAD)
          if P[i] >= P[(i-1) % (n*PAD)] and P[i] >= P[(i+1) % (n*PAD)]]
    pk.sort(reverse=True)
    pmax = pk[0][0]
    lines = []
    for p, om in pk:
        if p < 0.01*pmax: break
        if all(abs(om-o) > 4*2*math.pi/n for o, _ in lines):
            lines.append((om, p))
        if len(lines) >= 6: break
    tot = sum(p for _, p in lines)
    om0 = lines[0][0]
    om = refine(sig, om0)
    half = n//2
    err = abs(refine(sig[:half], om0)-refine(sig[half:], om0)) + 1e-9
    return om, err, lines[0][1]/tot

def cayley_menger_ok(N, E, L):
    D = np.zeros((N, N))
    for Le, (i, j) in zip(L, E):
        D[i, j] = D[j, i] = Le**2
    # 三角不等式
    tri = True
    for i in range(N):
        for j in range(i+1, N):
            for k in range(N):
                if k in (i, j): continue
                if math.sqrt(D[i, j]) > math.sqrt(D[i, k])+math.sqrt(D[k, j])+1e-12:
                    tri = False
    J = np.eye(N)-np.ones((N, N))/N; B = -0.5*J@D@J
    ev = np.linalg.eigvalsh(B); scale = max(abs(ev).max(), 1e-300)
    return tri, bool(ev.min()/scale > -1e-9), float(ev.min()/scale), int((ev/scale > 1e-8).sum())

md = ["# 辺ごとの波長読み（パス11）：ν_e → λ_e = λ0/|m_e| → L_e=(k_e/|m_e|)·λ0/2", ""]
rows = []
for tag in TARGETS:
    N = int(tag.split("N")[1]); E = edges(N); M = len(E)
    Z = np.load(os.path.join(DATA, tag, "states_treatment.npz"))["Z"][-LW:]
    amp = np.abs(Z).mean(axis=0)
    amax = amp.max()
    md.append(f"## {tag}（終窓、辺 {M} 本）"); md.append("")
    md.append("| 辺 (i,j) | 振幅/max | ν_e [rad/step] | ±err | 単色度 | 状態 |")
    md.append("|---|---|---|---|---|---|")
    freqs = []
    for e in range(M):
        if amp[e] < 1e-6*amax:
            md.append(f"| {E[e]} | {amp[e]/amax:.1e} | — | — | — | 波なし（長さ未定義） |")
            freqs.append(None); continue
        om, err, mono = edge_line(Z[:, e])
        freqs.append((om, err, mono))
        md.append(f"| {E[e]} | {amp[e]/amax:.3f} | {om:+.7f} | {err:.1e} | {mono:.3f} | {'単色' if mono > 0.9 else '混色'} |")
    md.append("")
    live = [(e, f) for e, f in enumerate(freqs) if f]
    if live:
        nu0 = min(abs(f[0]) for _, f in live if abs(f[0]) > 1e-6)
        rats = []
        for e, (om, err, mono) in live:
            r = abs(om)/nu0; er = err/nu0
            fr = Fraction(r).limit_denominator(12)
            lock = abs(r-float(fr)) < max(er, 1e-4)
            rats.append((e, r, er, fr, lock))
        md.append("**振動数比 |ν_e|/ν0（ν0=最小の辺振動数）と有理判定:** " +
                  "; ".join(f"{E[e]}:{r:.4f}±{er:.4f}" + (f"={fr.numerator}/{fr.denominator}" if lock else "(非有理域)")
                            for e, r, er, fr, lock in rats))
        # k≡1 の実単体テスト（波なし辺は除外できないので、全辺に波がある場合のみ）
        if all(f is not None for f in freqs) and all(abs(f[0]) > 1e-6 for f in freqs):
            L = [1.0/abs(f[0]) for f in freqs]
            Lmin = min(L); L = [x/Lmin for x in L]
            tri, psd, mineig, rank = cayley_menger_ok(N, E, L)
            md.append(f"**実単体テスト（k≡1、L_e∝1/|ν_e|）:** 三角不等式={'成立' if tri else '破れ'}、"
                      f"CM 半正定値={'成立' if psd else '破れ'}（min eig/scale={mineig:+.2e}）、rank={rank}")
            rows.append(dict(tag=tag, all_edges_wave=True, tri=tri, psd=psd, min_eig=mineig, rank=rank))
        else:
            nwave = sum(1 for f in freqs if f)
            md.append(f"**実単体テスト:** 波を持つ辺 {nwave}/{M} 本のみ → 波なし辺の長さが未定義のため k≡1 テスト不可"
                      "（星型＝衛星間の関係が消えた縮退幾何）")
            rows.append(dict(tag=tag, all_edges_wave=False, tri=None, psd=None, min_eig=None, rank=None))
    md.append("")
with open(os.path.join(ROOT, "results", "edge_wavelengths.md"), "w") as f:
    f.write("\n".join(md)+"\n")
with open(os.path.join(ROOT, "results", "edge_wavelengths_summary.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("PASS11 OK")
