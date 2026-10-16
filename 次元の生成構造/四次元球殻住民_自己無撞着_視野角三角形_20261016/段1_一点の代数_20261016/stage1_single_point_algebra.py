# -*- coding: utf-8 -*-
"""
段 1：一点の代数（設計書 v2 §4、判定は v1 §2 のとおり事前固定）
一つの複素四元数 v ∈ ℂ⁴ に、右の交換の回転 exp(θ R_u) と左のブースト exp(iη L_w) を同時に掛ける。
生成子 G = θ R_u + iη L_w（複素反対称）。初期状態は G の固有ベクトル（非零固有値）。exp(G) を 2n 回掛ける。
記録：Σz²（本当の二乗 vᵀv）、ノルム v†v、初期状態との重なり vᵀv₀、右の回転の固有成分（e^{±iθ} 側）の重み a, b と
v_read=(a²−b²)/(a²+b²)、invgamma_read=2ab/(a²+b²)。
判定（事前固定）：(a) |ΔΣz²|/ノルム ≤ 1e-13、(b) η=0 で n 回後に戻る（誤差 ≤ 1e-12）・η≠0 で戻らない、
(c) 右の回転の固有位相（周期 n）が η によらない、(d) v_read・invgamma_read と η の関係は記録のみ。
"""
import json, itertools
import numpy as np
from scipy.linalg import expm

def qmul(p, q):
    a1,b1,c1,d1 = p; a2,b2,c2,d2 = q
    return np.array([a1*a2-b1*b2-c1*c2-d1*d2, a1*b2+b1*a2+c1*d2-d1*c2, a1*c2-b1*d2+c1*a2+d1*b2, a1*d2+b1*c2-c1*b2+d1*a2])
E = np.eye(4)
AX = {'i': E[1], 'j': E[2], 'k': E[3]}
Lm = {n: np.array([qmul(u, E[c]) for c in range(4)]).T for n, u in AX.items()}   # 左の掛け算
Rm = {n: np.array([qmul(E[c], u) for c in range(4)]).T for n, u in AX.items()}   # 右の掛け算

CASES_THETA = {'23/124': (23, 124), '1/4': (1, 4), '1/6': (1, 6)}
ETAS = [0.0, 0.3, 0.8, 1.5]
out = []
for (tn, (m, n)), uR, uL, eta in itertools.product(CASES_THETA.items(), 'ijk', 'ijk', ETAS):
    theta = 2*np.pi*m/n
    G = theta*Rm[uR] + 1j*eta*Lm[uL]
    U = expm(G)
    UR = expm(theta*Rm[uR])                      # 右の回転だけ（固有成分の分解に使う）
    wR, VR = np.linalg.eig(UR)
    mu, V = np.linalg.eig(G)
    for e in range(4):
        if abs(mu[e]) < 1e-12: continue
        v0 = V[:, e]; v0 = v0/np.sqrt(np.vdot(v0, v0).real)
        v = v0.copy(); s0 = v0 @ v0; nrm0 = np.vdot(v0, v0).real
        dz = 0.0; ret_n = None
        for step in range(1, 2*n+1):
            v = U @ v
            dz = max(dz, abs(v @ v - s0) / max(nrm0, 1e-300))
            if step == n: ret_n = float(np.linalg.norm(v - v0))
        # 右の回転の固有成分：固有値 e^{+iθ} 側と e^{-iθ} 側への射影の重み（初期状態で読む）
        coef = np.linalg.solve(VR, v0)
        plus = sum(abs(coef[k])**2 for k in range(4) if np.angle(wR[k]) > 1e-9)
        minus = sum(abs(coef[k])**2 for k in range(4) if np.angle(wR[k]) < -1e-9)
        a, b = np.sqrt(plus), np.sqrt(minus)
        vr = (a*a - b*b)/(a*a + b*b) if a*a + b*b > 0 else float('nan')
        ig = 2*a*b/(a*a + b*b) if a*a + b*b > 0 else float('nan')
        out.append(dict(theta=tn, uR=uR, uL=uL, eta=eta, eig=[float(mu[e].real), float(mu[e].imag)],
                        sumz2_init=[float(s0.real), float(s0.imag)], max_rel_dsumz2=float(dz),
                        norm_after_2n=float(np.vdot(v, v).real), return_err_at_n=ret_n,
                        right_phase_per_step=float(theta), v_read=float(vr), invgamma_read=float(ig)))
json.dump(out, open('stage1_results.json', 'w'), ensure_ascii=False, indent=1)
# 判定の集計
a_ok = all(r['max_rel_dsumz2'] <= 1e-13 for r in out)
b_eta0 = [r for r in out if r['eta'] == 0.0]
b_ok0 = all(r['return_err_at_n'] <= 1e-12 for r in b_eta0)
b_okn = all(r['return_err_at_n'] > 1e-6 for r in out if r['eta'] != 0.0)
summary = dict(n_records=len(out),
               judge_a_sumz2_preserved=a_ok, max_rel_dsumz2=max(r['max_rel_dsumz2'] for r in out),
               judge_b_returns_eta0=b_ok0, max_return_err_eta0=max(r['return_err_at_n'] for r in b_eta0),
               judge_b_not_return_eta_nonzero=b_okn,
               min_return_err_eta_nonzero=min(r['return_err_at_n'] for r in out if r['eta'] != 0.0),
               sumz2_init_abs_max=max(abs(complex(*r['sumz2_init'])) for r in out),
               eigenvalues_seen=sorted({(round(r['eig'][0], 6), round(r['eig'][1], 6)) for r in out if r['theta'] == '23/124' and r['uR'] == 'i' and r['uL'] == 'i'}),
               v_invgamma_read_examples=[(r['theta'], r['uR'], r['uL'], r['eta'], round(r['v_read'], 6), round(r['invgamma_read'], 6)) for r in out if r['theta'] == '23/124' and r['uR'] == 'i' and r['uL'] in 'ij'][:16])
json.dump(summary, open('stage1_summary.json', 'w'), ensure_ascii=False, indent=1, default=str)
print(json.dumps(summary, ensure_ascii=False, indent=1, default=str))
