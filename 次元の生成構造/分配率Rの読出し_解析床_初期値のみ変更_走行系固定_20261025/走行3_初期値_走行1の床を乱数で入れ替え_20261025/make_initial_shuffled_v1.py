#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""初期値の生成（走行3 やり直し、2026-10-25）：走行1 の解析床（N²）の配列を乱数で入れ替えた床 npz を作る。
走行系には触れない。ea, eb（位置の構造）はそのまま、Z0・v・theta・family だけを同じ置換で並べ替える。置換は npz に記録。
使い方: python make_initial_shuffled_v1.py SRC_DIR OUT_DIR N1 N2 ..."""
import os, sys
import numpy as np

src, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
for N in [int(x) for x in sys.argv[3:]]:
    name = f'floor_90deg_all_N{N:05d}_analytic.npz'
    d = dict(np.load(os.path.join(src, name)))
    perm = np.random.default_rng(20261025 + N).permutation(d['Z0'].size)   # 旧走行3 と同じ乱数の種
    for k in ('Z0', 'v', 'theta', 'family'):
        d[k] = d[k][perm]
    d['g'] = d['g'][perm]
    d['perm'] = perm
    d['rule'] = np.array(str(d['rule']) + '_SHUFFLED_initial_value')
    np.savez_compressed(os.path.join(out, name), **d)
    print(f'N={N}: M={d["Z0"].size} shuffled (seed {20261025 + N}), ea/eb unchanged', flush=True)
print('INITIAL VALUES DONE')
