# -*- coding: utf-8 -*-
"""
検証：Σx_n^2 = R^2（本当の二乗）を光円錐として読み、速さ・1/γ・速度合成を場から読み出す（円周）
成分   x_n(φ) = iA (a e^{inφ/2} + b e^{-inφ/2})   n = 1..N    （a, b ≥ 0：± 回転項の重み）
二乗   x_n^2 = -A^2 (a^2 e^{inφ} + b^2 e^{-inφ} + 2ab)        （交差項 2ab は φ によらない＝空間に張り付く）
読出し R^2(φ) をフーリエ分解し（e^{inφ} の係数）
  P+ = Σ_{n>0} |c_n|, P- = Σ_{n<0} |c_n|, P0 = |c_0|
  方向成分 k = P+ − P-、張り付いた成分 m = P0、全体 R = P+ + P-
  v = k/R、1/γ = m/R、検証：v^2 + (1/γ)^2 = 1（その場の光速が 1）
速度合成：重み比 r = a/b を掛け合わせる（r_total = r1·r2）と v が相対論的合成則に従うか
"""
import numpy as np, json
A, N, M = 1.0, 18, 4096
phi = 2*np.pi*np.arange(M)/M
ns = np.arange(1, N+1)

def R2_field(a, b):
    return sum(-A**2*(a*np.exp(1j*n*phi/2) + b*np.exp(-1j*n*phi/2))**2 for n in ns)  # 本当の二乗の和

def readout(a, b):
    c = np.fft.fft(R2_field(a, b))/M                 # c[n] は e^{inφ} の係数
    P0 = abs(c[0]); Pp = np.abs(c[1:N+1]).sum(); Pm = np.abs(c[-N:]).sum()
    R = Pp + Pm
    return dict(a=a, b=b, Pplus=Pp, Pminus=Pm, P0=P0, v=(Pp-Pm)/R, inv_gamma=P0/R,
                cone_residual=((Pp+Pm)**2 - (Pp-Pm)**2 - P0**2)/R**2)

# 1) 重み比を走査
scan = []
for lr in np.linspace(-3, 3, 61):                    # ln(a/b)
    r = np.exp(lr); a, b = r/(1+r), 1/(1+r)
    d = readout(a, b); d["ln_ratio"] = float(lr); scan.append(d)
# 光子的（片方だけ）
photon = readout(1.0, 0.0)
# 2) 速度合成
adds = []
for l1, l2 in [(0.3, 0.5), (0.8, 0.8), (1.5, -0.4), (2.0, 1.0), (-1.2, 0.7)]:
    d1 = readout(np.exp(l1), 1.0); d2 = readout(np.exp(l2), 1.0); d12 = readout(np.exp(l1+l2), 1.0)
    rel = (d1["v"]+d2["v"])/(1+d1["v"]*d2["v"])
    adds.append(dict(ln_r1=l1, ln_r2=l2, v1=d1["v"], v2=d2["v"], v_composed_field=d12["v"],
                     v_relativistic_formula=rel, v_galilean=d1["v"]+d2["v"]))

out = dict(settings=dict(A=A, N=N, M=M, doc=__doc__), scan=scan, photon=photon, velocity_addition=adds)
json.dump(out, open("results_lightcone_readout_s1.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1, default=float)
print("max |cone residual| over scan: %.2e" % max(abs(d["cone_residual"]) for d in scan))
print("v vs tanh(ln r): max diff %.2e" % max(abs(d["v"]-np.tanh(d["ln_ratio"])) for d in scan))
print("photon: v=%.6f 1/gamma=%.6f" % (photon["v"], photon["inv_gamma"]))
print("rest (a=b):", [ (round(d["v"],6), round(d["inv_gamma"],6)) for d in scan if abs(d["ln_ratio"])<1e-9])
for d in adds:
    print("v1=%.4f v2=%.4f  field=%.6f  relativistic=%.6f  galilean=%.4f" % (d["v1"], d["v2"], d["v_composed_field"], d["v_relativistic_formula"], d["v_galilean"]))
