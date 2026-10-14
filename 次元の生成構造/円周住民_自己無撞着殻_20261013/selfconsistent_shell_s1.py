# -*- coding: utf-8 -*-
"""
円周住民の自己無撞着殻（二次元版）  2026-10-13

設定
  成分     x_n = iA cos(n (phi - phi0)/2)
  殻の形   sum_n x_n^2 + (iR)^2 = 0  ->  R^2 = A^2 sum_n cos^2(n (phi - phi0)/2)
  住民位相 phi = 2*pi * (ピーク点からの弧長) / (全周長)   背景角 theta は住民に見えない
  重心     殻（一様な線密度の閉曲線）の弧長重心。各方向の R はこの重心から測り直す
  反復     R(theta) が変わらなくなるまで（減衰係数 beta で混合）

並進を除いて形だけを追う: 毎回、殻を「重心から測った R」で組み直すと、
新しい殻の重心は構成中心から d だけずれる。形が収束したときの d を記録する。
d != 0 なら、重心から測り直すたびに殻全体が d だけ並進する（形は不変）。

系: odd (n=1,3,..,N) / even (n=2,4,..,N) / all (n=1..N)、A=1
"""
import numpy as np, json, sys

def components(system, N):
    return {'odd': np.arange(1, N+1, 2), 'even': np.arange(2, N+1, 2), 'all': np.arange(1, N+1)}[system]

def R_of_phase(dphi, ns, A=1.0):
    return A*np.sqrt(np.sum(np.cos(np.outer(ns, dphi)/2.0)**2, axis=0))

def arclength(P):
    d = np.diff(np.vstack([P, P[:1]]), axis=0)
    seg = np.hypot(d[:, 0], d[:, 1]); L = seg.sum()
    s = np.concatenate([[0.0], np.cumsum(seg)[:-1]])
    mid = 0.5*(P + np.roll(P, -1, axis=0))
    return s, L, (mid*seg[:, None]).sum(0)/L

def step(R, th, ns):
    """構成中心=原点の殻 R(theta) から、重心で測り直した新しい R を返す"""
    P = R[:, None]*np.c_[np.cos(th), np.sin(th)]
    s, L, cen = arclength(P)
    phi = 2*np.pi*(s % L)/L                         # theta=0 の点（ピーク）を位相の基準
    a = np.unwrap(np.arctan2(P[:, 1]-cen[1], P[:, 0]-cen[0]))
    fold = int(np.sum(np.diff(a) <= 0))             # 重心から見て方向が逆行する点の数（>0 なら殻が重心に対して星形でない）
    am = np.mod(a, 2*np.pi); o = np.argsort(am); a_s = am[o]; p_s = np.unwrap(phi[o])
    w = 2*np.pi*np.sign(p_s[-1]-p_s[0])
    phi_th = np.interp(th, np.r_[a_s-2*np.pi, a_s, a_s+2*np.pi], np.r_[p_s-w, p_s, p_s+w])
    return R_of_phase(phi_th, ns), cen, L, fold, phi

def solve(system, N, M=4096, beta=0.2, itmax=8000, tol=1e-11):
    ns = components(system, N)
    th = np.linspace(0, 2*np.pi, M, endpoint=False)
    R = np.ones(M); folds = 0; hist = []
    for it in range(itmax):
        Rn, cen, L, fold, phi = step(R, th, ns)
        e = float(np.max(np.abs(Rn-R))); hist.append(e); folds += fold
        R = (1-beta)*R + beta*Rn
        if e < tol: break
    Rn, cen, L, fold, phi = step(R, th, ns)
    return dict(system=system, N=N, ns=ns.tolist(), M=M, beta=beta, iterations=it+1,
                converged=bool(e < tol), final_residual=e, fold_points_final=fold, fold_points_total=folds,
                d=cen.tolist(), L=float(L), d_over_L=float(cen[0]/L),
                R_max=float(R.max()), R_min=float(R.min()), residual_history=hist), th, R, phi

CASES = {'odd_N1': ('odd', 1, 0.2), 'odd_N3': ('odd', 3, 0.2), 'odd_N17': ('odd', 17, 0.2),
         'even_N2': ('even', 2, 0.05), 'even_N4': ('even', 4, 0.05), 'even_N18': ('even', 18, 0.05),
         'all_N2': ('all', 2, 0.2), 'all_N4': ('all', 4, 0.2), 'all_N18': ('all', 18, 0.2)}

if __name__ == '__main__':
    for key in sys.argv[1:]:
        system, N, beta = CASES[key]
        res, th, R, phi = solve(system, N, beta=beta)
        np.savez(f'shell_{key}.npz', theta=th, R=R, phi=phi, d=np.array(res['d']))
        json.dump(res, open(f'result_{key}.json', 'w'), ensure_ascii=False, indent=1)
        print(key, {k: res[k] for k in ('iterations', 'converged', 'final_residual', 'fold_points_final', 'd', 'L', 'd_over_L', 'R_max', 'R_min')})
