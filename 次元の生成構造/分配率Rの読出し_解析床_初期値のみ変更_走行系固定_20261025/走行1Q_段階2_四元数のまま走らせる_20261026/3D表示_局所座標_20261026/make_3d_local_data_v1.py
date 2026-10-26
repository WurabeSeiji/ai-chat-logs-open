#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""段階2 本走行の状態を「局所 3 次元座標＋振幅」の 3D 表示用データに書き出す（読み出しのみ、2026-10-26）。

各歩・各関係の四元数 q = a + b i + c j + d k を極形式 q = |q| e^(θ n) で読む：
  |q| = √(a²+b²+c²+d²)                       振幅
  θ  = atan2(|(b,c,d)|, a) ∈ [0, π]          1 からの角度
  n  = (b,c,d)/|(b,c,d)|                      虚軸の向き（単位ベクトル）
  位置 = θ n ∈ 半径 π の球の中（i, j, k 成分の向きに沿う 3 次元座標）
n が決まらない点（(b,c,d) がほぼ 0）：θ≈0（q≈|q|）は原点、θ≈π（q≈−|q|）は (π, 0, 0) に置く（−1 は半径 π の球面全体が同じ点。
1–i 面の組では −1 = e^(iπ) なので i 軸の向きを代表に選ぶ）。
大きさ用に |q| / max|q|（その走行の全歩での最大）。群：初期状態で 1–i 面（終点 w 偶数）＝0、j–k 面（w 奇数）＝1。
出力：data/run_N{N}_den{den}.js（座標と振幅を int16 に量子化して base64、座標の分解能 π/32767）
使い方: python make_3d_local_data_v1.py STATES_DIR STATE0_DIR OUT_DIR N1 N2 ..."""
import base64, json, math, os, sys
import numpy as np

SCALE = 32767.0 / math.pi

if __name__ == '__main__':
    sdir, s0dir, odir = sys.argv[1:4]
    Ns = [int(x) for x in sys.argv[4:]]
    os.makedirs(odir, exist_ok=True)
    index = []
    for N in Ns:
        w = np.load(os.path.join(s0dir, f'floor_90deg_all_N{N:05d}_analytic.npz'))['eb']
        group = (w % 2).astype(int).tolist()
        for den in (N, 124):
            Q = np.load(os.path.join(sdir, f'hm_N{N}_den_{den}_states_500.npz'))['Q']   # (S+1, M, 4)
            amp = np.linalg.norm(Q, axis=2)
            V = Q[..., 1:]
            vn = np.linalg.norm(V, axis=2)
            theta = np.arctan2(vn, Q[..., 0])
            small = vn <= 1e-12 * amp
            with np.errstate(divide='ignore', invalid='ignore'):
                X = V / vn[..., None] * theta[..., None]
            X[small] = 0.0
            X[small & (Q[..., 0] < 0)] = np.array([math.pi, 0.0, 0.0])
            A = amp / amp.max()
            arr = np.concatenate([np.round(X * SCALE), np.round(A[..., None] * 32767.0)], axis=2).astype('<i2')
            key = f'N{N}_den{den}'
            payload = dict(N=N, den=den, M=int(Q.shape[1]), steps=int(Q.shape[0] - 1), scale=SCALE, group=group,
                           n_undefined_direction=int(small.sum()), b64=base64.b64encode(arr.tobytes()).decode())
            with open(os.path.join(odir, f'run_{key}.js'), 'w') as f:
                f.write('window.RUNS = window.RUNS || {};\nwindow.RUNS["%s"] = %s;\n' % (key, json.dumps(payload)))
            index.append(dict(key=key, N=N, den=den, M=payload['M'], n_undefined_direction=payload['n_undefined_direction'],
                              undefined_at_step0=int(small[0].sum())))
            print(key, 'M', payload['M'], 'n undefined (all steps)', payload['n_undefined_direction'], 'at step0', int(small[0].sum()), flush=True)
    json.dump(index, open(os.path.join(odir, 'index.json'), 'w'), indent=1)
    print('DATA DONE')
