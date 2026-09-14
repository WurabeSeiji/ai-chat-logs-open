#!/usr/bin/env python3
"""5seed 有効標本評価（承認済み手順の最終段・読出しのみ・次元推定なし）。

各 N の保存済み終端 Z（第一段 Tier-0 合格版）に対して:
  1. Bob-star 不変量（|P_k|・基準次数 r・φ_{k|r}）を全 5seed×N Bob で読出す
  2. seed 内 Bob 対距離 vs seed 間距離を測り、Bob 多重度が実効標本かを判定
  3. n_raw=5N、重複（<1e-12）除去後の n_unique を実測
  4. 設計書 v1.2 §8.3 の事前登録予測（|P1|,|P2| 縮退・|P4| 相対大・r=4 浮上）を検定
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
TAU_REF = 1e-8
DUP_TOL = 1e-12
K_EVAL = 8  # 距離比較用の共通打切り（記録は k<=8; 完全フレームは主推定で使用）

def edges(N):
    a, b = np.triu_indices(N, k=1)
    return a.astype(np.int64), b.astype(np.int64)

def Pk(z, k):
    H = np.sum(np.abs(z)**2)
    return np.sum(z**k) / (H**(k/2))

def star_frame(zs, Kmax):
    P = [Pk(zs, k) for k in range(1, Kmax+1)]
    absP = np.array([abs(p) for p in P])
    r = next((k for k in range(1, Kmax+1) if absP[k-1] > TAU_REF), None)
    vec = list(absP)
    if r is not None:
        for k in range(1, Kmax+1):
            if k == r:
                continue
            ph = np.angle((P[k-1]**r) * np.conj(P[r-1]**k))
            vec += [np.cos(ph), np.sin(ph)]
    return r, absP, np.asarray(vec, float)

def main():
    rows = []
    star_rows = []
    for N in range(6, 41):
        d = np.load(HERE/'terminal_Z'/f'N{N}_terminal_Z.npz')
        Zt = d['Z_terminal']  # (5, M)
        ea, eb = edges(N)
        Keff = min(K_EVAL, N-1)
        frames = {}
        for s in range(Zt.shape[0]):
            for v in range(N):
                zs = Zt[s][(ea == v) | (eb == v)]
                r, absP, vec = star_frame(zs, Keff)
                frames[(s, v)] = (r, vec)
                star_rows.append({'N': N, 'seed': s, 'Bob': v, 'r': r,
                                  'absP1': absP[0], 'absP2': absP[1],
                                  'absP3': absP[2] if Keff >= 3 else np.nan,
                                  'absP4': absP[3] if Keff >= 4 else np.nan,
                                  'kappa_B': absP[1]})
        keys = list(frames)
        rs = [frames[k][0] for k in keys]
        r_mode = max(set(rs), key=rs.count)
        within, across = [], []
        dup = 0; pairs = 0
        for i in range(len(keys)):
            ri, xi = frames[keys[i]]
            for j in range(i+1, len(keys)):
                rj, xj = frames[keys[j]]
                if ri != rj or len(xi) != len(xj):
                    continue
                dd = float(np.linalg.norm(xi - xj)); pairs += 1
                if dd < DUP_TOL:
                    dup += 1
                (within if keys[i][0] == keys[j][0] else across).append(dd)
        # n_unique: 重複クラスタを1つに数える（貪欲）
        used = set(); n_unique = 0
        for i in range(len(keys)):
            if i in used:
                continue
            n_unique += 1
            ri, xi = frames[keys[i]]
            for j in range(i+1, len(keys)):
                rj, xj = frames[keys[j]]
                if ri == rj and len(xi) == len(xj) and np.linalg.norm(xi-xj) < DUP_TOL:
                    used.add(j)
        rows.append({
            'N': N, 'n_raw': 5*N, 'n_unique': n_unique,
            'r_mode': r_mode, 'r_all_same': len(set(rs)) == 1,
            'within_seed_dist_median': float(np.median(within)) if within else np.nan,
            'within_seed_dist_min': float(np.min(within)) if within else np.nan,
            'across_seed_dist_median': float(np.median(across)) if across else np.nan,
            'within_over_across': (float(np.median(within))/float(np.median(across))
                                   if within and across else np.nan),
            'duplicate_pairs': dup, 'comparable_pairs': pairs,
        })
        print(f'N={N:3d} n_raw={5*N:4d} n_unique={n_unique:4d} r_mode={r_mode} '
              f'within_med={rows[-1]["within_seed_dist_median"]:.3g} '
              f'across_med={rows[-1]["across_seed_dist_median"]:.3g}', flush=True)
    df = pd.DataFrame(rows); sf = pd.DataFrame(star_rows)
    df.to_csv(HERE/'effective_samples_by_N.csv', index=False)
    sf.to_csv(HERE/'terminal_star_invariants_all.csv', index=False)
    pred = {
        'prediction_registered': 'terminal: |P1|,|P2| small; |P4| relatively large; r=4 may emerge',
        'median_absP1': float(sf.absP1.median()), 'median_absP2': float(sf.absP2.median()),
        'median_absP3': float(sf.absP3.median()), 'median_absP4': float(sf.absP4.median()),
        'median_kappa_B': float(sf.kappa_B.median()),
        'max_kappa_B': float(sf.kappa_B.max()),
        'r_distribution_overall': {str(k): int((sf.r == k).sum()) for k in sorted(sf.r.dropna().unique())},
        'n_stars_total': int(len(sf)),
    }
    (HERE/'effective_samples_summary.json').write_text(
        json.dumps(pred, indent=2, ensure_ascii=False), encoding='utf-8')
    print(json.dumps(pred, indent=2, ensure_ascii=False))

if __name__ == '__main__':
    main()
