# -*- coding: utf-8 -*-
"""
四次元球殻（S^3）の住民の自己無撞着殻  2026-10-13

ピーク一つの波はピーク軸の周りに対称なので、殻は R^4 内の回転超曲面
  (x, rho*omega),  omega in S^2,  x = r(chi) cos chi,  rho = r(chi) sin chi,  chi in [0, pi]
とし、子午線 r(chi) だけを求める（三つの theta0 は軸の向きを決めるだけ）。

  住民位相   psi = pi * s(chi) / s(pi)   s はピーク（chi=0）からの子午線弧長＝測地距離
  殻の形     sum x_n^2 + (iR)^2 = 0 -> R^2 = sum_n |x_n|^2 / A^2 * A^2
     (a) x_n = iA cos(n psi)                         円周の cos(n phi/2) で phi/2 -> psi
     (b) x_n = iA U_n(cos psi)/(n+1)                 S^3 の固有モード、ピーク値を A にそろえる
                 U_n(cos psi) = sin((n+1) psi)/sin psi
  重心       体積要素 dV = 4 pi rho^2 ds による重心（軸上）。各方向の R はこの重心から測り直す
  反復       r(chi) が変わらなくなるまで（減衰係数 beta）。重心のずれ d を記録

系: odd (n=1,3,..,N) / even (n=2,4,..,N) / all (n=1..N)、A=1
"""
import numpy as np, json, sys

def components(system, N):
    return {'odd': np.arange(1, N+1, 2), 'even': np.arange(2, N+1, 2), 'all': np.arange(1, N+1)}[system]

def R_of_psi(psi, ns, form):
    psi = np.clip(psi, 0, np.pi)
    if form == 'a':
        v = np.cos(np.outer(ns, psi))
    else:
        s = np.sin(psi); n1 = (ns+1)[:, None]
        with np.errstate(invalid='ignore', divide='ignore'):
            v = np.sin(n1*psi)/(n1*s)
        # 端点の極限 U_n(±1)/(n+1) = (±1)^n
        v[:, psi < 1e-12] = 1.0
        v[:, psi > np.pi-1e-12] = ((-1.0)**ns)[:, None]
    return np.sqrt(np.sum(v**2, axis=0))

def geometry(R, chi):
    x, rho = R*np.cos(chi), R*np.sin(chi)
    dx, dr = np.diff(x), np.diff(rho); ds = np.hypot(dx, dr)
    s = np.concatenate([[0.0], np.cumsum(ds)]); S = s[-1]
    xm, rm = 0.5*(x[1:]+x[:-1]), 0.5*(rho[1:]+rho[:-1])
    w = rm**2*ds
    return x, rho, s, S, float((xm*w).sum()/w.sum())

def step(R, chi, ns, form):
    x, rho, s, S, xb = geometry(R, chi)
    psi = np.pi*s/S
    a = np.arctan2(rho, x-xb)                     # 重心から見た極角（0..pi）
    fold = int(np.sum(np.diff(a) <= 0))
    o = np.argsort(a)
    psi_chi = np.interp(chi, a[o], psi[o])
    return R_of_psi(psi_chi, ns, form), xb, S, fold, psi

def solve(form, system, N, M=4096, beta=0.2, itmax=8000, tol=1e-11):
    ns = components(system, N)
    chi = np.linspace(0, np.pi, M)
    R = np.ones(M); hist = []
    for it in range(itmax):
        Rn, xb, S, fold, psi = step(R, chi, ns, form)
        e = float(np.max(np.abs(Rn-R))); hist.append(e)
        R = (1-beta)*R + beta*Rn
        if e < tol: break
    Rn, xb, S, fold, psi = step(R, chi, ns, form)
    return dict(form=form, system=system, N=N, ns=ns.tolist(), M=M, beta=beta, iterations=it+1,
                converged=bool(e < tol), final_residual=e, fold_points_final=fold,
                d=xb, S_meridian=float(S), d_over_S=float(xb/S),
                R_peak=float(R[0]), R_equator=float(np.interp(np.pi/2, psi, R)), R_antipode=float(R[-1]),
                R_max=float(R.max()), R_min=float(R.min()), residual_history=hist), chi, R, psi

CASES = [(f, s, N) for f in ('a', 'b') for s, N in (('odd', 1), ('odd', 3), ('odd', 17), ('even', 2), ('even', 4), ('even', 18), ('all', 2), ('all', 4), ('all', 18))]

if __name__ == '__main__':
    sel = sys.argv[1] if len(sys.argv) > 1 else 'ab'
    for f, s, N in CASES:
        if f not in sel: continue
        key = f'{f}_{s}_N{N}'
        res, chi, R, psi = solve(f, s, N)
        np.savez(f'shell3_{key}.npz', chi=chi, R=R, psi=psi, d=res['d'])
        json.dump(res, open(f'result3_{key}.json', 'w'), ensure_ascii=False)
        print(key, 'it', res['iterations'], 'conv', res['converged'], 'res %.1e' % res['final_residual'], 'fold', res['fold_points_final'],
              'd %+.2e' % res['d'], 'S %.4f' % res['S_meridian'], 'R0 %.4f Req %.4f Rpi %.4f' % (res['R_peak'], res['R_equator'], res['R_antipode']))
