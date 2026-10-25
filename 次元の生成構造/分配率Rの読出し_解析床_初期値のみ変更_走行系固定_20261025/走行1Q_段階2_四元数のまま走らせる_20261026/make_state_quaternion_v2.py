#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""走行1Q 段階2（2026-10-26）の初期状態を作る。入力は段階1 の floors_init_quaternion/（q＝e^(iπv/N)·r·e^(jπw/N) と φ0）。

段階1 の復号（make_initial_quaternion_v1.decode、段階1 からのバイト同一コピー）で q の値だけから v, w, r₁ を読み、次の 2 種類の状態を書く。
  main   ：q = e^(iφ0) · i^v · r · j^w、 r = r₁·(−1)^⌊(v+w)/2⌋
            （j を i に潰すと e^(iφ0)·r₁·i^((v+w) mod 2) ＝ 走行1 の床 Z0 になる。これを検査する）
  control：q = Re Z0 + i·Im Z0（走行1 の床を 1–i 面にそのまま置く。j, k 成分 0）
辺表 ea, eb は段階1 で値から読んだもの（＝v, w）。
使い方: python make_state_quaternion_v2.py STAGE1_FLOOR_DIR OUT_MAIN OUT_CONTROL N1 N2 ..."""
import json
import math
import os
import sys

import numpy as np

from make_initial_quaternion_v1 import decode, qmul, exp_i, exp_j

if __name__ == '__main__':
    src, out_main, out_ctrl = sys.argv[1:4]
    Ns = [int(x) for x in sys.argv[4:]]
    os.makedirs(out_main, exist_ok=True); os.makedirs(out_ctrl, exist_ok=True)
    manifest = {}
    print(' N     M | 潰した main − Z0 | |Σq²|/H main  control | 成分エネルギー比 main (1,i,j,k)')
    for N in Ns:
        name = f'floor_90deg_all_N{N:05d}_analytic.npz'
        d = np.load(os.path.join(src, name))
        Q1 = np.asarray(d['Q'], dtype=np.float64); phi0 = float(d['phi0'])
        v, w, r1, z, _, _, _ = decode(Q1, N, phi0)
        r = r1 * (-1.0) ** ((v + w) // 2)
        rq = np.zeros((len(r), 4)); rq[:, 0] = r
        Qm = qmul(qmul(exp_i(phi0 + 0.5 * math.pi * v), rq), exp_j(0.5 * math.pi * w))
        collapsed = (Qm[:, 0] + 1j * Qm[:, 1]) + 1j * (Qm[:, 2] + 1j * Qm[:, 3])   # j→i：a + b i + c i + d i·i
        dcol = float(np.max(np.abs(collapsed - z)))
        Qc = np.zeros_like(Qm); Qc[:, 0] = z.real; Qc[:, 1] = z.imag

        def sqsum(Q):
            s = np.zeros(4)
            for q in Q:
                s += qmul(q, q)
            return float(np.linalg.norm(s) / np.sum(Q * Q))
        Em = np.sum(Qm * Qm, axis=0) / np.sum(Qm * Qm)
        for out, Q, tag in ((out_main, Qm, 'main'), (out_ctrl, Qc, 'control')):
            np.savez_compressed(os.path.join(out, name), Q0=Q, ea=v, eb=w, n=np.array(N), M=np.array(len(v)), phi0=np.array(phi0),
                                rule=np.array(f'RUN1Q stage2 {tag}: ' + ('q=e^(i phi0) i^v r j^w, r=r1(-1)^floor((v+w)/2)' if tag == 'main' else 'q=Re Z0 + i Im Z0 (run1 floor in 1-i plane)')))
        chk = dict(N=N, M=int(len(v)), collapse_main_minus_Z0=dcol, sqsum_main=sqsum(Qm), sqsum_control=sqsum(Qc),
                   energy_fraction_main_1ijk=[float(x) for x in Em])
        manifest[name] = chk
        print(f"{N:>2} {len(v):>5} | {dcol:.1e} | {chk['sqsum_main']:.2e} {chk['sqsum_control']:.2e} | " + ' '.join(f'{x:.3f}' for x in Em), flush=True)
    json.dump(manifest, open(os.path.join(out_main, 'manifest_state_quaternion.json'), 'w'), indent=2, ensure_ascii=False)
    print('QUATERNION STATES DONE')
