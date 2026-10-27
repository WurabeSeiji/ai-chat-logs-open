#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""段階2 本走行の状態（results_main/hm_N{N}_den_{den}_states_500.npz）を、複素図形式の 3D 表示用データに書き出す（読み出しのみ、2026-10-26）。

各歩・各関係の四元数 q = a + b i + c j + d k について、位置
    X = |q| · (b, c, d) / a
を倍精度で計算して float32 で保存する（中心投影の座標 (b,c,d)/a に振幅 |q| を掛けたもの）。
±q は同じ点なので、a < 0 なら q を −q に置き換えてから計算する（a ≥ 0 に揃える）。
|a| ≤ 1e-14·|q| の関係は実部が丸めの範囲で 0（純虚）とみなし、X を Infinity にする（無限遠。表示側で個数を題に出す）。
群：初期状態で 1–i 面（終点 w 偶数）＝0、j–k 面（w 奇数）＝1。
出力：data/run_N{N}_den{den}.js（1 点 3 値の float32 を base64）
使い方: python make_3d_central_data_v1.py STATES_DIR STATE0_DIR OUT_DIR N1 N2 ..."""
import base64, json, os, sys
import numpy as np

ZERO_REL = 1e-14

if __name__ == '__main__':
    sdir, s0dir, odir = sys.argv[1:4]
    Ns = [int(x) for x in sys.argv[4:]]
    os.makedirs(odir, exist_ok=True)
    index = []
    for N in Ns:
        w = np.load(os.path.join(s0dir, f'floor_90deg_all_N{N:05d}_analytic.npz'))['eb']
        group = (w % 2).astype(int).tolist()
        for den in (N, 124):
            Q = np.load(os.path.join(sdir, f'hm_N{N}_den_{den}_states_500.npz'))['Q'].astype(np.float64)   # (S+1, M, 4)
            amp = np.linalg.norm(Q, axis=2)
            sgn = np.where(Q[..., 0] < 0, -1.0, 1.0)
            Qs = Q * sgn[..., None]
            a = Qs[..., 0]
            zero = a <= ZERO_REL * amp
            with np.errstate(divide='ignore', invalid='ignore'):
                X = amp[..., None] * Qs[..., 1:] / a[..., None]
            X[zero] = np.inf
            arr = X.astype('<f4')
            key = f'N{N}_den{den}'
            payload = dict(N=N, den=den, M=int(Q.shape[1]), steps=int(Q.shape[0] - 1), group=group,
                           n_inf_step0=int(zero[0].sum()), n_inf_total=int(zero.sum()),
                           b64=base64.b64encode(arr.tobytes()).decode())
            with open(os.path.join(odir, f'run_{key}.js'), 'w') as f:
                f.write('window.RUNS = window.RUNS || {};\nwindow.RUNS["%s"] = %s;\n' % (key, json.dumps(payload)))
            index.append(dict(key=key, N=N, den=den, M=payload['M'], steps=payload['steps'], n_inf_step0=payload['n_inf_step0'], n_inf_total=payload['n_inf_total']))
            print(key, 'M', payload['M'], 'inf at step0', payload['n_inf_step0'], 'inf total', payload['n_inf_total'], flush=True)
    json.dump(index, open(os.path.join(odir, 'index.json'), 'w'), indent=1)
    print('DATA DONE')
