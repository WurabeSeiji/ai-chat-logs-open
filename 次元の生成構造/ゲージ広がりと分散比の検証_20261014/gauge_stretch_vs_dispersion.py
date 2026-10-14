# -*- coding: utf-8 -*-
"""
検証：ゲージの広がり s と ω/k = sqrt(1 + m^2/k^2) は一致するか
定義（このファイル内で使うものだけ）
  R^2 = Σ x_n^2（本当の二乗）、x_n = iA×成分  -> |R| = sqrt(-R^2)
  ゲージ目盛 = 最高次成分の節（位相で等間隔）
  s_int = ピークを含む目盛区間の単位位相あたり測地線長 / 全体平均（測地線長 = ∫|R| dθ）
  s_pt  = ピーク点の |R| / 位相平均の |R|
  比較対象 ω/k = sqrt(1 + 1/κ_N^2)
     m = 1/|R|（双曲側の固有値ギャップ (k^2+1)/|R|^2 から）、k = κ_N/|R|
     κ_N：S1 の cos(nφ/2) は N/2、S3 の cos(nψ) は N、S3 の U_N は N+1
"""
import numpy as np, json
A = 1.0
def R_abs_S1(phi, ns):  return np.sqrt(A**2*sum(np.cos(n*phi/2)**2 for n in ns))
def R_abs_S3a(psi, ns): return np.sqrt(A**2*sum(np.cos(n*psi)**2 for n in ns))
def R_abs_S3b(psi, ns):
    s = np.sin(psi); s = np.where(np.abs(s) < 1e-12, 1e-12, s)
    return np.sqrt(A**2*sum((np.sin((n+1)*psi)/s/(n+1))**2 for n in ns))

def stretch(Rf, ns, marks, lo, hi, peak=0.0, M=400000):
    ds = []
    for a, b in zip(marks[:-1], marks[1:]):
        t = np.linspace(a, b, max(200, M//len(marks)))
        ds.append(np.trapezoid(Rf(t, ns), t))
    ds = np.array(ds)/np.diff(marks); c = 0.5*(marks[:-1]+marks[1:])   # 単位位相あたり（端の半区間も公平に）
    ip = np.argmin(np.abs(c-peak))
    ds_mean = np.sum(ds*np.diff(marks))/(marks[-1]-marks[0])
    t = np.linspace(lo, hi, M); mean_R = np.trapezoid(Rf(t, ns), t)/(hi-lo)
    return float(ds[ip]/ds_mean), float(Rf(np.array([peak+1e-9]), ns)[0]/mean_R)

rows = []
for N in [2, 4, 6, 8, 10, 12, 14, 16, 18]:
    ns = list(range(1, N+1))
    # S1：一周、節 (2k+1)π/N
    k = np.arange(N+1) - (N+1)//2; m1 = (2*k+1)*np.pi/N
    s1i, s1p = stretch(R_abs_S1, ns, m1, -np.pi, np.pi)
    # S3(a)：ピーク ψ=0 から対蹠点 π、節 (2k+1)π/(2N)、ピーク区間は [0, π/(2N)]
    m3a = np.concatenate([[0.0], (2*np.arange(N)+1)*np.pi/(2*N), [np.pi]])
    s3ai, s3ap = stretch(R_abs_S3a, ns, m3a, 0, np.pi)
    # S3(b)：節 kπ/(N+1)
    m3b = np.arange(N+2)*np.pi/(N+1)
    s3bi, s3bp = stretch(R_abs_S3b, ns, m3b, 1e-6, np.pi-1e-6)
    rows.append(dict(N=N,
        S1=dict(s_int=s1i, s_pt=s1p, w_over_k=float(np.sqrt(1+1/(N/2)**2))),
        S3a=dict(s_int=s3ai, s_pt=s3ap, w_over_k=float(np.sqrt(1+1/N**2))),
        S3b=dict(s_int=s3bi, s_pt=s3bp, w_over_k=float(np.sqrt(1+1/(N+1)**2)))))
json.dump(dict(definitions=__doc__, results=rows), open("results_gauge_stretch_vs_dispersion.json","w",encoding="utf-8"), ensure_ascii=False, indent=1)
print(" N | S1 s_int s_pt  w/k | S3a s_int s_pt  w/k | S3b s_int s_pt  w/k")
for r in rows:
    print("%2d | %.3f %.3f %.3f | %.3f %.3f %.3f | %.3f %.3f %.3f" % (r["N"],
      r["S1"]["s_int"], r["S1"]["s_pt"], r["S1"]["w_over_k"],
      r["S3a"]["s_int"], r["S3a"]["s_pt"], r["S3a"]["w_over_k"],
      r["S3b"]["s_int"], r["S3b"]["s_pt"], r["S3b"]["w_over_k"]))
