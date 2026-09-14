#!/usr/bin/env python3
"""Phase B 主推定（5seed・凍結仕様 v1.2 §13）— d_est(N,K) と κ_B null 対照。

凍結済み仕様（変更しない）:
  標本: N固定 × 最終state（stage1 terminal_Z 正本）× 5seed × N Bob
  chart: r = min{k: |P_k| > 1e-8}
  座標: A系 = (|P_1|,|P_k|,cosφ_k,sinφ_k)、B系 = (log(|P_k|+1e-12), cosφ_k, sinφ_k)
  標準化: S1（成分ごと平均0分散1）
  重複除去: 標準化前距離 < 1e-12
  K: 2..N-1 全走査
  推定器: E1 PCA / E2 local PCA / E3 TwoNN / E4 Levina-Bickel MLE
  CI: seed-cluster bootstrap 優先

実行前に事前宣言する実装細部（2026-09-14、結果を見る前に固定）:
  - τ_σ = 1e-10: 標準化時、σ < τ_σ の成分は除外し除外数を記録
  - E2 local PCA: 近傍 k_nn ∈ {10,20,30,50}（k_nn < n のもののみ）、局所固有値の累積寄与 95% での
    rank、全点の中央値を d_E2 とし k_nn 中央値を報告
  - E3 TwoNN: μ=r2/r1、μ の大きい方 10% を棄却し d = n_used / Σ log μ
  - E4 MLE: k ∈ {10,20}、Levina-Bickel 逆平均を点平均し k 平均
  - K*: 最小の K で |d_cons(K) − d_cons(K_max)| ≤ 0.25 が K 以降すべて成立
  - d_cons(K) = median(E2, E3, E4)（A系）。§20 のとおり PCA 単独では判定しない
  - bootstrap: seed 5個を復元抽出→ユニーク化（≥2 seed を要求、満たさなければ引き直し）→
    重複除去→ K=K_max の d_cons を再計算、B=200、percentile 2.5/97.5
  - κ_B null 対照: P(κ ≤ κ_obs) = 1 − (1−κ²)^{(m−1)/2}（m=N−1）の log10 を報告
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
TERM = HERE.parent / 'stage1_ssb_terminal_Z_v1' / 'terminal_Z'
TAU_REF = 1e-8
EPS = 1e-12
TAU_SIG = 1e-10
DUP_TOL = 1e-12
KNN_LIST = (10, 20, 30, 50)
MLE_KS = (10, 20)
BOOT = 200
RNG = np.random.default_rng(20260914)

def edges(N):
    a, b = np.triu_indices(N, k=1)
    return a.astype(np.int64), b.astype(np.int64)

def star_P(zs, Kmax):
    H = np.sum(np.abs(zs)**2)
    return np.array([np.sum(zs**k) / (H**(k/2)) for k in range(1, Kmax+1)])

def frames_for_N(N):
    """全 (seed,Bob) の P ベクトル（K_max=N-1）と r を返す。"""
    d = np.load(TERM / f'N{N}_terminal_Z.npz')
    Zt = d['Z_terminal']
    ea, eb = edges(N)
    Kmax = N - 1
    Ps, rs, seeds, kappaB = [], [], [], []
    for s in range(Zt.shape[0]):
        for v in range(N):
            zs = Zt[s][(ea == v) | (eb == v)]
            P = star_P(zs, Kmax)
            absP = np.abs(P)
            r = next((k for k in range(1, Kmax+1) if absP[k-1] > TAU_REF), None)
            Ps.append(P); rs.append(r); seeds.append(s); kappaB.append(absP[1])
    return np.array(Ps), np.array(rs, object), np.array(seeds), np.array(kappaB), Kmax

def coords(P, r, K, system):
    """1標本の座標ベクトル（K 打切り）。"""
    absP = np.abs(P[:K])
    if system == 'A':
        vec = list(absP)
    else:
        vec = list(np.log(absP + EPS))
    Pr = P[r-1]
    for k in range(1, K+1):
        if k == r:
            continue
        ph = np.angle((P[k-1]**r) * np.conj(Pr**k))
        vec += [np.cos(ph), np.sin(ph)]
    return np.asarray(vec, float)

def build_matrix(Ps, rs, K, system):
    X = np.array([coords(P, r, K, system) for P, r in zip(Ps, rs)])
    mu = X.mean(0); sig = X.std(0)
    keep = sig >= TAU_SIG
    Xs = (X[:, keep] - mu[keep]) / sig[keep]
    return Xs, int((~keep).sum()), X.shape[1]

def dedup(X, seeds):
    keep = []
    for i in range(len(X)):
        if all(np.linalg.norm(X[i] - X[j]) >= DUP_TOL for j in keep):
            keep.append(i)
    return X[keep], seeds[np.array(keep)]

def pca_ranks(X):
    sv = np.linalg.svd(X - X.mean(0), compute_uv=False)
    en = sv**2; cum = np.cumsum(en) / en.sum()
    return {f'pca_{p}': int(np.searchsorted(cum, p) + 1) for p in (0.95, 0.99, 0.999)}, sv

def local_pca(X):
    n = len(X)
    D = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=2)
    med_ranks = []
    for knn in KNN_LIST:
        if knn + 1 >= n:
            continue
        ranks = []
        for i in range(n):
            idx = np.argsort(D[i])[1:knn+1]
            Y = X[idx] - X[idx].mean(0)
            ev = np.linalg.svd(Y, compute_uv=False)**2
            cum = np.cumsum(ev) / ev.sum()
            ranks.append(int(np.searchsorted(cum, 0.95) + 1))
        med_ranks.append(float(np.median(ranks)))
    return float(np.median(med_ranks)) if med_ranks else np.nan

def twonn(X):
    n = len(X)
    D = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=2)
    mus = []
    for i in range(n):
        ds = np.sort(D[i])[1:3]
        if ds[0] > 0:
            mus.append(ds[1] / ds[0])
    mus = np.sort(np.array(mus))
    used = mus[:max(2, int(0.9*len(mus)))]
    s = np.sum(np.log(used))
    return float(len(used) / s) if s > 0 else np.nan

def lb_mle(X):
    n = len(X)
    D = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=2)
    ests = []
    for k in MLE_KS:
        if k + 1 >= n:
            continue
        vals = []
        for i in range(n):
            ds = np.sort(D[i])[1:k+1]
            if ds[0] <= 0:
                continue
            Tk = ds[-1]
            v = np.mean(np.log(Tk / ds[:-1]))
            if v > 0:
                vals.append(1.0 / v)
        if vals:
            ests.append(float(np.mean(vals)))
    return float(np.mean(ests)) if ests else np.nan

def d_consensus(X):
    e2, e3, e4 = local_pca(X), twonn(X), lb_mle(X)
    return float(np.nanmedian([e2, e3, e4])), e2, e3, e4

def main():
    curve_rows, main_rows, kap_rows = [], [], []
    for N in range(6, 41):
        Ps, rs, seeds, kappaB, Kmax = frames_for_N(N)
        assert all(r == 1 for r in rs)
        m = N - 1
        # κ_B null 対照（主次元推定と独立）
        med = float(np.median(kappaB))
        logp = ((m-1)/2) * np.log10(max(1.0 - med**2, 1e-300)) if med < 1 else 0.0
        logp = float(np.log10(1.0 - (1.0 - med**2)**((m-1)/2))) if med < 0.1 else logp
        kap_rows.append({'N': N, 'kappa_B_median': med, 'kappa_B_max': float(kappaB.max()),
                         'null_exact_mean': float(np.sqrt(np.pi)/2 *
                                                  np.exp(__import__('math').lgamma((m+1)/2) -
                                                         __import__('math').lgamma((m+2)/2))),
                         'log10_P_null_kappa_le_median': logp})
        # K 走査（A/B 両系）
        dA_last = None
        for system in ('A', 'B'):
            for K in range(2, Kmax+1):
                Xs, nexc, Damb = build_matrix(Ps, rs, K, system)
                Xd, sd = dedup(Xs, seeds)
                dcons, e2, e3, e4 = d_consensus(Xd)
                pr, _ = pca_ranks(Xd)
                curve_rows.append({'N': N, 'system': system, 'K': K, 'n': len(Xd),
                                   'D_ambient': Damb, 'n_over_D': len(Xd)/Damb,
                                   'excluded_const': nexc, 'd_consensus': dcons,
                                   'localPCA': e2, 'TwoNN': e3, 'MLE': e4, **pr})
            if system == 'A':
                dA_last = [r for r in curve_rows if r['N'] == N and r['system'] == 'A']
        # K*（A系 d_consensus 基準）
        dmaxv = dA_last[-1]['d_consensus']
        Kstar = None
        for row in dA_last:
            if all(abs(r2['d_consensus'] - dmaxv) <= 0.25 for r2 in dA_last if r2['K'] >= row['K']):
                Kstar = row['K']; break
        # bootstrap（K=Kmax、A系）
        boots = []
        XsA, _, _ = build_matrix(Ps, rs, Kmax, 'A')
        for _ in range(BOOT):
            while True:
                pick = RNG.integers(0, 5, 5)
                u = np.unique(pick)
                if len(u) >= 2:
                    break
            mask = np.isin(seeds, u)
            Xb, _ = dedup(XsA[mask], seeds[mask])
            if len(Xb) > 12:
                boots.append(d_consensus(Xb)[0])
        lo, hi = (np.percentile(boots, [2.5, 97.5]) if boots else (np.nan, np.nan))
        last = dA_last[-1]
        lastB = [r for r in curve_rows if r['N'] == N and r['system'] == 'B'][-1]
        main_rows.append({'N': N, 'n': last['n'], 'Kstar': Kstar,
                          'd_consensus_A': last['d_consensus'], 'd_consensus_B': lastB['d_consensus'],
                          'localPCA_A': last['localPCA'], 'TwoNN_A': last['TwoNN'], 'MLE_A': last['MLE'],
                          'pca95_A': last['pca_0.95'],
                          'boot_lo': float(lo), 'boot_hi': float(hi)})
        print(f"N={N:3d} n={last['n']:4d} K*={Kstar} dA={last['d_consensus']:.2f} "
              f"dB={lastB['d_consensus']:.2f} (E2={last['localPCA']:.2f} E3={last['TwoNN']:.2f} "
              f"E4={last['MLE']:.2f} pca95={last['pca_0.95']}) CI=[{lo:.2f},{hi:.2f}]", flush=True)
    pd.DataFrame(curve_rows).to_csv(HERE/'d_est_N_K_curves.csv', index=False)
    pd.DataFrame(main_rows).to_csv(HERE/'main_results_by_N.csv', index=False)
    pd.DataFrame(kap_rows).to_csv(HERE/'kappaB_null_contrast.csv', index=False)
    print('DONE')

if __name__ == '__main__':
    main()
