#!/usr/bin/env python3
"""Phase B pilot — long10000 N=6 read-only readout validation.

目的（設計書 v1.2 §17/§23、木原指示 2026-09-14）:
  Bob-star 抽出・辺index対応・P_k・基準次数 r・|P1|,|P2|,|P3|,|P4|・
  Bob 間重複率の検査のみを行う。次元推定は行わず、主張しない。

対象データ（読み取り専用・変更しない）:
  N3_N40_long10000_20260905/results/hm_N6_den_6_states_10000.npz
  辺規約: np.triu_indices(N, k=1)（生成プログラム run_N3_N40_steps10000_v1.py の edges() と同一）
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = ('/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/'
        'マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/電子の反跳実験/'
        'seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831')
NPZ = Path(ROOT) / 'N3_N40_long10000_20260905/results/hm_N6_den_6_states_10000.npz'
OUT = Path(__file__).resolve().parent
TAU_REF = 1e-8
DUP_TOL = 1e-12

def edges(N):
    a, b = np.triu_indices(N, k=1)
    return a.astype(np.int64), b.astype(np.int64)

def Pk(z, k):
    H = np.sum(np.abs(z)**2)
    return np.sum(z**k) / (H**(k/2))

def star_invariants(z, K):
    """|P_1..P_K|, selected r, phi_{k|r} (charge-balanced)."""
    P = [Pk(z, k) for k in range(1, K+1)]
    absP = [abs(p) for p in P]
    r = next((k for k in range(1, K+1) if absP[k-1] > TAU_REF), None)
    phis = {}
    if r is not None:
        Pr = P[r-1]
        for k in range(1, K+1):
            if k == r:
                continue
            phis[k] = float(np.angle((P[k-1]**r) * np.conj(Pr**k)))
    return absP, r, phis

def main():
    d = np.load(NPZ, allow_pickle=False)
    Z = d['Z']; N = int(d['N']); steps = int(d['steps'])
    M = N*(N-1)//2
    assert Z.shape == (steps+1, M), (Z.shape, steps, M)
    ea, eb = edges(N)

    # --- 辺規約検査: 各頂点の star がちょうど N-1 本、全 star の合併が各辺を2回被覆
    counts = np.zeros(M, int)
    star_idx = {}
    for v in range(N):
        idx = np.flatnonzero((ea == v) | (eb == v))
        assert len(idx) == N-1, (v, len(idx))
        star_idx[v] = idx
        counts[idx] += 1
    assert np.all(counts == 2)

    zfin = Z[-1]; z0 = Z[0]
    K = N - 1
    rows = []
    vecs = {}
    for v in range(N):
        zs = zfin[star_idx[v]]
        absP, r, phis = star_invariants(zs, K)
        H_B = float(np.sum(np.abs(zs)**2))
        row = {'Bob': v, 'H_B': H_B, 'r_selected': r, 'kappa_B': absP[1]}
        for k in range(1, K+1):
            row[f'absP{k}'] = absP[k-1]
        for k, ph in phis.items():
            row[f'phi_{k}_ref{r}'] = ph
        rows.append(row)
        # Bob 間比較ベクトル（A系 raw: |P_k| と (cos,sin) of phi_{k|r}）
        vec = list(absP)
        for k in sorted(phis):
            vec += [np.cos(phis[k]), np.sin(phis[k])]
        vecs[v] = (r, np.asarray(vec, float))

    df = pd.DataFrame(rows)

    # --- Bob 間距離（同一 r の Bob 同士のみ直接比較可能）
    dist_rows = []
    for i in range(N):
        for j in range(i+1, N):
            ri, xi = vecs[i]; rj, xj = vecs[j]
            comparable = (ri == rj) and (len(xi) == len(xj))
            dd = float(np.linalg.norm(xi - xj)) if comparable else np.nan
            dist_rows.append({'Bob_i': i, 'Bob_j': j, 'same_r': comparable, 'distance': dd})
    ddf = pd.DataFrame(dist_rows)
    comp = ddf[ddf.same_r]

    # --- 大域量（sanity）: 力学は Z^T Z と ‖Z‖ を厳密保存するはず
    def kap_global(z):
        return float(abs(np.sum(z*z)) / np.sum(np.abs(z)**2))
    summary = {
        'N': N, 'M': M, 'steps': steps,
        'edge_convention': 'np.triu_indices(N,k=1) row-major, verified: each star has N-1 edges, union covers each edge twice',
        'tau_ref': TAU_REF, 'dup_tol': DUP_TOL,
        'global_kappa_step0': kap_global(z0),
        'global_kappa_final': kap_global(zfin),
        'H_total_step0': float(np.sum(np.abs(z0)**2)),
        'H_total_final': float(np.sum(np.abs(zfin)**2)),
        'r_distribution': {str(k): int((df.r_selected == k).sum()) for k in sorted(df.r_selected.dropna().unique())},
        'absP1_range': [float(df.absP1.min()), float(df.absP1.max())],
        'absP2_range': [float(df.absP2.min()), float(df.absP2.max())],
        'absP3_range': [float(df.absP3.min()), float(df.absP3.max())],
        'absP4_range': [float(df.absP4.min()), float(df.absP4.max())],
        'kappa_B_range': [float(df.kappa_B.min()), float(df.kappa_B.max())],
        'bob_pairs_total': int(len(ddf)),
        'bob_pairs_comparable_same_r': int(len(comp)),
        'bob_pair_distance_min': float(comp.distance.min()) if len(comp) else None,
        'bob_pair_distance_median': float(comp.distance.median()) if len(comp) else None,
        'bob_pair_distance_max': float(comp.distance.max()) if len(comp) else None,
        'duplicate_pairs_lt_tol': int((comp.distance < DUP_TOL).sum()) if len(comp) else 0,
    }

    df.to_csv(OUT/'pilot_bobstar_invariants_N6.csv', index=False)
    ddf.to_csv(OUT/'pilot_bob_pair_distances_N6.csv', index=False)
    (OUT/'pilot_summary_N6.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding='utf-8')
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print()
    print(df.to_string(index=False))

if __name__ == '__main__':
    main()
