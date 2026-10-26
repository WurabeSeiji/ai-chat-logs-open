#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""段階2 本走行の状態（results_main/hm_N{N}_den_{den}_states_500.npz）を 3D 表示用のデータに書き出す（読み出しのみ、2026-10-26）。

各歩・各関係の四元数 q = a + b i + c j + d k について：
  u = q/|q|（S³ 上の点）
  u' = p · u,  p = (1 + i + j + k)/2   （固定の回転。初期値は ±1, ±i, ±j, ±k に乗り、−1 が射影の極で無限遠に飛ぶため、
                                       極を初期値から外す。p は等長変換なので形は変わらない）
  (x, y, z) = (b', c', d') / (1 + a')   （−1 を極とする立体射影）
  半径 RCLIP を超えた点は、向きを保って半径 RCLIP に寄せる（clip フラグを残す）
  大きさ用に |q| / max|q|（その走行の全歩での最大）
群：初期状態で 1–i 面（終点 w 偶数）＝0、j–k 面（w 奇数）＝1。
出力：data/run_N{N}_den{den}.js（window.RUNS[key] = {N, den, M, steps, group:[...], b64: Int16 (x,y,z,amp) を base64}）
使い方: python make_3d_data_v1.py STATES_DIR STATE0_DIR OUT_DIR N1 N2 ..."""
import base64, json, os, sys
import numpy as np

RCLIP = 4.0
SCALE = 32767.0 / RCLIP


def qmul(p, q):
    a1, b1, c1, d1 = np.moveaxis(p, -1, 0); a2, b2, c2, d2 = np.moveaxis(q, -1, 0)
    return np.stack([a1*a2 - b1*b2 - c1*c2 - d1*d2, a1*b2 + b1*a2 + c1*d2 - d1*c2,
                     a1*c2 - b1*d2 + c1*a2 + d1*b2, a1*d2 + b1*c2 - c1*b2 + d1*a2], axis=-1)


P = np.array([0.5, 0.5, 0.5, 0.5])

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
            U = Q / amp[..., None]
            U2 = qmul(np.broadcast_to(P, U.shape), U)
            den_ = 1.0 + U2[..., 0]
            with np.errstate(divide='ignore', invalid='ignore'):
                X = U2[..., 1:] / den_[..., None]
            X[~np.isfinite(X)] = RCLIP
            r = np.linalg.norm(X, axis=2)
            clip = r > RCLIP
            X[clip] *= (RCLIP / r[clip])[:, None]
            A = amp / amp.max()
            arr = np.concatenate([np.round(X * SCALE), np.round(A[..., None] * 32767.0)], axis=2).astype('<i2')
            key = f'N{N}_den{den}'
            payload = dict(N=N, den=den, M=int(Q.shape[1]), steps=int(Q.shape[0] - 1), rclip=RCLIP, scale=SCALE,
                           group=group, nclip=int(clip.sum()), b64=base64.b64encode(arr.tobytes()).decode())
            with open(os.path.join(odir, f'run_{key}.js'), 'w') as f:
                f.write('window.RUNS = window.RUNS || {};\nwindow.RUNS["%s"] = %s;\n' % (key, json.dumps(payload)))
            index.append(dict(key=key, N=N, den=den, M=payload['M'], clipped_points=payload['nclip']))
            print(key, 'M', payload['M'], 'clipped', payload['nclip'], 'of', clip.size, flush=True)
    json.dump(index, open(os.path.join(odir, 'index.json'), 'w'), indent=1)
    print('DATA DONE')
