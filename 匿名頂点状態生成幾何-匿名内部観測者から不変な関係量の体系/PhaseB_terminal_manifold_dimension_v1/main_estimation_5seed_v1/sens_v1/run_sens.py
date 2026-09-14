#!/usr/bin/env python3
"""SENS-3〜6（＋SENS-2 は窓データ生成後に別段実行）— rank-4 非支持の定性的頑健性検査。

事前宣言（2026-09-14、結果を見る前に固定）:
  - 代表 N = {6,10,16,20,30,40}。目的は定性判定（d>4 か・K飽和なしか）であり精密値ではない。
  - 構成: {A,B} × {S0=無標準化, S1=z-score, S2=(x-median)/MAD} と、B系の ε ∈ {1e-10,1e-12,1e-14}。
  - 各構成で K=N-1 の E2/E3/E4 と d_cons=median(E2,E3,E4) を報告。
    K曲線は (A,S1)/(A,S2)/(B,S1) で出し、K飽和なしの頑健性を見る。
  - 内蔵対照: (A,S1,1e-12) の d_cons が主推定 main_results_by_N.csv と一致（<1e-9）すること。
  - 判定は「d_cons>4 が全構成で維持されるか」「d(N) 増加傾向が維持されるか」のみ。
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
MAIN = HERE.parent
TERM = MAIN.parent / 'stage1_ssb_terminal_Z_v1' / 'terminal_Z'
REP_N = [6, 10, 16, 20, 30, 40]
TAU_REF = 1e-8
TAU_SIG = 1e-10
DUP_TOL = 1e-12
KNN_LIST = (10, 20, 30, 50)
MLE_KS = (10, 20)

def edges(N):
    a, b = np.triu_indices(N, k=1)
    return a.astype(np.int64), b.astype(np.int64)

def frames_for_N(N):
    d = np.load(TERM / f'N{N}_terminal_Z.npz')
    Zt = d['Z_terminal']; ea, eb = edges(N); Kmax = N - 1
    Ps, seeds = [], []
    for s in range(Zt.shape[0]):
        for v in range(N):
            zs = Zt[s][(ea == v) | (eb == v)]
            H = np.sum(np.abs(zs)**2)
            Ps.append(np.array([np.sum(zs**k)/(H**(k/2)) for k in range(1, Kmax+1)]))
            seeds.append(s)
    return np.array(Ps), np.array(seeds), Kmax

def coords(P, K, system, eps):
    absP = np.abs(P[:K])
    vec = list(absP) if system == 'A' else list(np.log(absP + eps))
    P1 = P[0]
    for k in range(2, K+1):
        ph = np.angle(P[k-1] * np.conj(P1**k))
        vec += [np.cos(ph), np.sin(ph)]
    return np.asarray(vec, float)

def standardize(X, mode):
    if mode == 'S0':
        return X, 0
    if mode == 'S1':
        mu, sig = X.mean(0), X.std(0)
    else:
        mu = np.median(X, 0)
        sig = np.median(np.abs(X - mu), 0) * 1.4826
    keep = sig >= TAU_SIG
    return (X[:, keep] - mu[keep]) / sig[keep], int((~keep).sum())

def dedup(X):
    keep = []
    for i in range(len(X)):
        if all(np.linalg.norm(X[i]-X[j]) >= DUP_TOL for j in keep):
            keep.append(i)
    return X[keep]

def dists(X):
    return np.linalg.norm(X[:, None, :] - X[None, :, :], axis=2)

def local_pca(X):
    n = len(X); D = dists(X); meds = []
    for knn in KNN_LIST:
        if knn + 1 >= n:
            continue
        rk = []
        for i in range(n):
            idx = np.argsort(D[i])[1:knn+1]
            Y = X[idx] - X[idx].mean(0)
            ev = np.linalg.svd(Y, compute_uv=False)**2
            rk.append(int(np.searchsorted(np.cumsum(ev)/ev.sum(), 0.95) + 1))
        meds.append(float(np.median(rk)))
    return float(np.median(meds)) if meds else np.nan

def twonn(X):
    D = dists(X); mus = []
    for i in range(len(X)):
        ds = np.sort(D[i])[1:3]
        if ds[0] > 0:
            mus.append(ds[1]/ds[0])
    mus = np.sort(np.array(mus)); used = mus[:max(2, int(0.9*len(mus)))]
    s = np.sum(np.log(used))
    return float(len(used)/s) if s > 0 else np.nan

def lb_mle(X):
    D = dists(X); ests = []
    for k in MLE_KS:
        if k + 1 >= len(X):
            continue
        vals = []
        for i in range(len(X)):
            ds = np.sort(D[i])[1:k+1]
            if ds[0] <= 0:
                continue
            v = np.mean(np.log(ds[-1]/ds[:-1]))
            if v > 0:
                vals.append(1.0/v)
        if vals:
            ests.append(float(np.mean(vals)))
    return float(np.mean(ests)) if ests else np.nan

def d_all(X):
    e2, e3, e4 = local_pca(X), twonn(X), lb_mle(X)
    return float(np.nanmedian([e2, e3, e4])), e2, e3, e4

def estimate(Ps, K, system, smode, eps):
    X = np.array([coords(P, K, system, eps) for P in Ps])
    Xs, nexc = standardize(X, smode)
    Xd = dedup(Xs)
    return d_all(Xd), nexc, len(Xd)

def main():
    main_csv = pd.read_csv(MAIN/'main_results_by_N.csv').set_index('N')
    rows, curves = [], []
    for N in REP_N:
        Ps, seeds, Kmax = frames_for_N(N)
        configs = [('A','S0',1e-12), ('A','S1',1e-12), ('A','S2',1e-12),
                   ('B','S0',1e-12), ('B','S1',1e-12), ('B','S2',1e-12),
                   ('B','S1',1e-10), ('B','S1',1e-14)]
        for system, smode, eps in configs:
            (dc, e2, e3, e4), nexc, n = estimate(Ps, Kmax, system, smode, eps)
            rows.append({'N': N, 'system': system, 'std': smode, 'eps': eps, 'K': Kmax,
                         'n': n, 'excluded': nexc, 'd_cons': dc,
                         'localPCA': e2, 'TwoNN': e3, 'MLE': e4, 'gt4': bool(dc > 4)})
            if system == 'A' and smode == 'S1':
                ref = float(main_csv.loc[N, 'd_consensus_A'])
                assert abs(dc - ref) < 1e-9, (N, dc, ref)
        for system, smode in [('A','S1'), ('A','S2'), ('B','S1')]:
            for K in range(2, Kmax+1):
                (dc, *_), _, n = estimate(Ps, K, system, smode, 1e-12)
                curves.append({'N': N, 'system': system, 'std': smode, 'K': K, 'd_cons': dc})
        print(f'N={N} done', flush=True)
    df = pd.DataFrame(rows); cf = pd.DataFrame(curves)
    df.to_csv(HERE/'sens_3456_results.csv', index=False)
    cf.to_csv(HERE/'sens_K_curves.csv', index=False)
    summary = {
        'all_configs_gt4': bool(df.gt4.all()),
        'min_d_cons': float(df.d_cons.min()),
        'min_d_cons_config': df.loc[df.d_cons.idxmin()][['N','system','std','eps']].to_dict(),
        'increase_trend_by_config': {
            f'{s}-{m}-{e:g}': bool(g.sort_values('N').d_cons.iloc[-1] > g.sort_values('N').d_cons.iloc[0] + 2)
            for (s, m, e), g in df.groupby(['system','std','eps'])},
        'baseline_matches_main': True,
    }
    (HERE/'sens_summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False, default=str),
                                          encoding='utf-8')
    print(json.dumps(summary, indent=2, ensure_ascii=False, default=str))

if __name__ == '__main__':
    main()
