#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""90度床の解析的生成器・辺集合一般化版・走行2 版（2026-10-25）: 隣接重みづけなし（adjacency が全対 1）。それ以外は走行1 の floor_analytic_v1.py と同一。

論文3（2026-09-09/12）の `make_90deg_floor_analytic_v1.py`（同梱 original_reference_…）の `integer_K` / `releq_residual` /
`build_analytic_floor` を文字どおり写し、変更は次の 2 点だけ：
  (1) 辺集合を引数で選ぶ：triu（上三角 i<j、M=N(N−1)/2、論文3と同一）／ all（対角込みの全順序対、M=N²）
  (2) adjacency を辺リスト (ea, eb) から作る（規則は正本と同一：添字を共有すれば 1、対角 0。triu では正本 adjacency(N) と一致）
ラベル則 k_e=(i+j) mod 4、整数 K、iK の σ_max 固有ベクトル、離散ラベル反復（4象限へ再量子化、大域位相 4 通りで前ラベルに最一致、
固定点／ゲージ等価サイクルで停止、最小残差の床を採用）は論文3と同じ。one_step は 2026-09-05 正本から SHA 照合で抽出（論文3と同じ）。
種は置かない（論文3 の床 npz と同じく g=0、Z0=z）。
"""
import hashlib
import json
import math
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
CANON = ('/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/'
         '次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831/'
         'N3_N40_stage123_sweep_20260905/run_N3_N40_stage123_v1.py')
PROG_SHA = '1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567'

src = open(CANON, encoding='utf-8').read()
assert hashlib.sha256(src.encode()).hexdigest() == PROG_SHA, '正本SHA不一致'
g = {'__name__': 'canonical_funcs', '__file__': CANON}
exec(compile(src.split('rows=[]; summaries=[]')[0], CANON, 'exec'), g)
one_step = g['one_step']
canonical_adjacency = g['adjacency']


def edges_triu(N):
    a, b = np.triu_indices(N, k=1)
    return a.astype(np.int64), b.astype(np.int64)


def edges_all(N):
    a, b = np.meshgrid(np.arange(N), np.arange(N), indexing='ij')
    return a.ravel().astype(np.int64), b.ravel().astype(np.int64)


def adjacency_from_edges(ea, eb):
    """走行2（2026-10-25）: 隣接重みづけをやめる。全対 1、対角 0（対角は integer_K の sin(0)=0 と one_step の対角ゼロ化で元々効かない）。ea, eb は M を決めるためだけに使う。"""
    M = len(ea)
    A = np.ones((M, M), dtype=np.float64)
    np.fill_diagonal(A, 0.0)
    return A


def integer_K(labels, A):
    diff = (labels[None, :] - labels[:, None]) % 4
    s = np.zeros_like(diff, dtype=float)
    s[diff == 1] = 1.0
    s[diff == 3] = -1.0
    return A * s


def releq_residual(z, A, N):
    w = one_step(z, A, N)
    ip = np.vdot(z, w)
    return float(np.linalg.norm(w - np.exp(1j * np.angle(ip)) * z) / np.linalg.norm(z))


def build_analytic_floor(N, edges_fn, max_iter=500):
    ea, eb = edges_fn(N)
    A = adjacency_from_edges(ea, eb)
    labels = ((ea + eb) % 4).astype(int)          # 決定論的種
    best = None                                    # (residual, z, labels)
    seen = {}
    for it in range(max_iter):
        K = integer_K(labels, A)
        w, V = np.linalg.eigh(1j * K)
        v = V[:, -1]                               # σ_max
        z = v / np.linalg.norm(v)
        res = releq_residual(z, A, N)
        if best is None or res < best[0]:
            best = (res, z.copy(), labels.copy())
        newlab = (np.round(np.degrees(np.angle(v)) % 360 / 90).astype(int)) % 4
        cand = None; agmax = -1
        for sh in range(4):
            c = (newlab + sh) % 4; ag = int(np.sum(c == labels))
            if ag > agmax: agmax = ag; cand = c
        key = tuple(cand.tolist())
        if np.array_equal(cand, labels) or key in seen:
            break
        seen[key] = it
        labels = cand
    res, z, labels = best
    a, b = z.real, z.imag
    ph = np.degrees(np.angle(z * np.exp(-1j * np.angle(z[np.argmax(np.abs(z))])))) % 360
    Kf = integer_K(labels, A)
    chk = dict(N=N, M=len(z), iters=it,
               norm2=float(np.vdot(z, z).real), centroid_abs=float(abs(z.sum())),
               re_closure=float(a @ a - b @ b), im_closure_cross=float(a @ b),
               closure_abs=float(abs(np.sum(z * z))), releq_residual=res,
               on_90deg_axes=bool(all(min(abs((p - x + 180) % 360 - 180) for x in (0, 90, 180, 270)) < 1e-8 for p in ph)),
               n_quadrant=len(set((np.round(ph / 90).astype(int) % 4).tolist())),
               n_distinct_amp=len(set(np.round(np.abs(z), 9).tolist())),
               sigma=float(np.max(np.linalg.eigvalsh(1j * Kf))),
               rank_K=int(np.linalg.matrix_rank(Kf, tol=1e-9)),
               centroid_exact_zero=bool(N % 4 == 0),
               seed_labels_kept=bool(np.array_equal(labels, ((ea + eb) % 4).astype(int))))
    return z, labels, chk, ea, eb


if __name__ == '__main__':
    structure, outdir = sys.argv[1], sys.argv[2]
    Ns = [int(x) for x in sys.argv[3:]]
    fn = {'triu': edges_triu, 'all': edges_all}[structure]
    os.makedirs(outdir, exist_ok=True)
    manifest = {}
    print('N    M | 象限 軸上 振幅種 | |Σz| |Σz²| | 残差 σ rank | 反復 種ラベル保持')
    for N in Ns:
        z, labels, chk, ea, eb = build_analytic_floor(N, fn)
        if structure == 'triu':
            assert np.array_equal(adjacency_from_edges(ea, eb), canonical_adjacency(N))
        name = f'floor_90deg_{structure}_N{N:05d}_analytic.npz'
        np.savez_compressed(os.path.join(outdir, name), v=z, g=np.zeros(len(z), dtype=np.complex128), Z0=z,
                            sigma=np.array([chk['sigma']]), residual=np.array(chk['releq_residual']),
                            n=np.array(N), seed=np.array(-1), delta=np.array(0.0), tol=np.array(0.0),
                            iters=np.array(chk['iters']), family=labels.astype(np.int64), theta=np.angle(z),
                            rule=np.array(f'theoretical_floor_90deg_integerK_sigmamax_analytic_{structure}'),
                            ea=ea, eb=eb, M=np.array(len(z)), structure=np.array(structure))
        manifest[name] = chk
        print(f"{N:>2} {chk['M']:>4} | {chk['n_quadrant']:>3} {str(chk['on_90deg_axes']):>5} {chk['n_distinct_amp']:>3} | "
              f"{chk['centroid_abs']:.3f} {chk['closure_abs']:.1e} | {chk['releq_residual']:.1e} {chk['sigma']:.4f} {chk['rank_K']:>3} | "
              f"{chk['iters']:>2} {chk['seed_labels_kept']}", flush=True)
    json.dump(manifest, open(os.path.join(outdir, f'manifest_floor_{structure}.json'), 'w'), indent=2, ensure_ascii=False)
    print('FLOORS DONE')
