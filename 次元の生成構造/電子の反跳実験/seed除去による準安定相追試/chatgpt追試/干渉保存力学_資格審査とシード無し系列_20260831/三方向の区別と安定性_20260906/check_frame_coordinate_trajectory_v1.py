#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""三方向の区別と安定性 (4): 3枠座標での状態軌跡の読出し（読み出しのみ）。

問い（木原 2026-09-06）: 3方向が安定なら、3次元空間としての読出しが可能になるか。
最小の実演: 状態自身を枠座標で読む。
  X(t) = |⟨p, Z(t)⟩|²/‖Z‖²、Y(t) = |⟨q, Z(t)⟩|²/‖Z‖²、Zc(t) = 1−X−Y（面外成分。
  N=3 では厳密に核＝第三方向。N≥4 では核＋他平面の合算——ラベルは「面外」）。
検定（N=3 の予想）: 初期 (1/2, 1/2, 0) → 終状態 (1/3, 1/2, 1/6)（実験11-4の帰結）。
もし Y(t)=1/2 が全時間で成立するなら、インフレーションは「X が Z へ 1/6 を渡し、
Y は傍観者」という**一軸方向のエネルギー輸送**として読める——3次元読出しの実演。
入力: N3_N40_long10000_20260905 の npz（SHA台帳照合）。新規走行なし。
出力: check_frame_coordinate_trajectory_v1.{json,png}
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent
SWEEP = BASE.parent / "N3_N40_long10000_20260905"
NS = [3, 4, 5, 6, 10, 40]

ledger = {}
with open(SWEEP / "SHA256SUMS.txt") as fh:
    for line in fh:
        parts = line.split()
        if len(parts) == 2:
            ledger[parts[1]] = parts[0]

out = {}
fig, axs = plt.subplots(2, 3, figsize=(16, 8))
for ax, N in zip(axs.ravel(), NS):
    rel = f'results/hm_N{N}_den_{N}_states_10000.npz'
    h = hashlib.sha256(open(SWEEP / rel, 'rb').read()).hexdigest()
    assert ledger[rel] == h, f'INPUT GATE FAIL: {rel}'
    Z = np.asarray(np.load(SWEEP / rel)['Z'], dtype=np.complex128)
    rp, ip = Z[0].real, Z[0].imag
    p = rp / np.linalg.norm(rp)
    q = ip - (ip @ p) * p
    q = q / np.linalg.norm(q)
    H = np.sum(np.abs(Z) ** 2, axis=1)
    X = np.abs(Z @ p) ** 2 / H
    Y = np.abs(Z @ q) ** 2 / H
    Zc = 1.0 - X - Y
    res = {'start_XYZ': [float(X[0]), float(Y[0]), float(Zc[0])],
           'final_XYZ_mean_last1000': [float(X[-1000:].mean()), float(Y[-1000:].mean()),
                                       float(Zc[-1000:].mean())],
           'Y_deviation_from_half': {'max_all_time': float(np.abs(Y - 0.5).max()),
                                     'argmax_step': int(np.argmax(np.abs(Y - 0.5))),
                                     'final_mean_last1000': float(np.abs(Y[-1000:] - 0.5).mean())},
           'X_final_times_6': float(6 * X[-1000:].mean()),
           'Y_final_times_6': float(6 * Y[-1000:].mean()),
           'Z_final_times_6': float(6 * Zc[-1000:].mean())}
    out[str(N)] = res
    ax.plot(X, lw=0.8, label='X = |<p,Z>|^2/H')
    ax.plot(Y, lw=0.8, label='Y = |<q,Z>|^2/H')
    ax.plot(Zc, lw=0.8, label='Z = out-of-plane')
    for yv in (1 / 3, 1 / 2, 1 / 6):
        ax.axhline(yv, color='gray', ls=':', lw=0.5)
    ax.set_xscale('log')
    ax.set_title(f'N={N}')
    ax.grid(alpha=.3)
    ax.legend(fontsize=7)
    print(f"N={N}: 始点 XYZ=({X[0]:.6f},{Y[0]:.6f},{Zc[0]:.2e}) → 終点 "
          f"({X[-1000:].mean():.8f},{Y[-1000:].mean():.8f},{Zc[-1000:].mean():.8f}) | "
          f"6倍=({6*X[-1000:].mean():.5f},{6*Y[-1000:].mean():.5f},{6*Zc[-1000:].mean():.5f}) | "
          f"|Y−1/2| 全時間max={np.abs(Y-0.5).max():.3e}(step {np.argmax(np.abs(Y-0.5))})", flush=True)

with open(BASE / 'check_frame_coordinate_trajectory_v1.json', 'w') as fh:
    json.dump(out, fh, indent=2)
fig.suptitle('State trajectory in its own 3-frame coordinates (X, Y, Z fractions)', y=.995)
fig.tight_layout()
fig.savefig(BASE / 'check_frame_coordinate_trajectory_v1.png', dpi=140)
print('ALL DONE')
