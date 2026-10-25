#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""初期値の生成（走行4 やり直し、2026-10-25）：走行1 の解析床（N²）から配列の最終要素を 1 つ除いた床 npz を作る。
走行系には触れない。Z0・v・theta・family・g と、辺リスト ea, eb の同じ末尾を落とす（M → M−1）。残る値は床のまま（再構成しない）。
N² の並び e=i·N+j で最終要素は自己対 (N−1, N−1)。
使い方: python make_initial_droplast_v1.py SRC_DIR OUT_DIR N1 N2 ..."""
import os, sys
import numpy as np

src, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
for N in [int(x) for x in sys.argv[3:]]:
    name = f'floor_90deg_all_N{N:05d}_analytic.npz'
    d = dict(np.load(os.path.join(src, name)))
    removed = (int(d['ea'][-1]), int(d['eb'][-1]), complex(d['Z0'][-1]))
    for k in ('Z0', 'v', 'theta', 'family', 'g', 'ea', 'eb'):
        d[k] = d[k][:-1]
    d['M'] = np.array(d['Z0'].size)
    d['removed_edge'] = np.array([removed[0], removed[1]])
    d['removed_value'] = np.array(removed[2])
    d['rule'] = np.array(str(d['rule']) + '_LAST_ELEMENT_REMOVED_initial_value')
    np.savez_compressed(os.path.join(out, name), **d)
    print(f'N={N}: M={d["Z0"].size} (removed edge ({removed[0]},{removed[1]}), |z|^2={abs(removed[2])**2:.3e})', flush=True)
print('INITIAL VALUES DONE')
