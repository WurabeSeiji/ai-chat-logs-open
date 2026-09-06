#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""三方向の区別と安定性 (2): 元の複素平面 Π0 内で保存される方向の探索（読み出しのみ）。

問い（木原 2026-09-06）: 元の複素平面で保存される方向があれば、法線と合わせて
3軸が安定して読み出せるのではないか。

面内状態は複素2成分 v(t) = (p·Z(t), q·Z(t)) ∈ C²。面内の「方向」は CP¹ 点
（グローバル位相を除いた v̂ の向き）。面内回転 SO(2) の複素固有ベクトルは
円偏光対 ê± = (1, ±i)/√2 で、固有値は純位相——つまり **90°位相差構造（円偏光
の組み合わせ）そのものが、面内回転で唯一向きを変えない方向**。親 Z0 = Re+i·Im
（Re⊥Im・等振幅）はちょうど ê₊ に対応するはず——これを測る。
検定:
 (a) c₊(t) = |⟨ê₊, v̂(t)⟩|: 円偏光方向の純度（1=完全保存）。全歩。
 (b) CP¹ 連続ドリフト d_cont(t) = 1−|⟨v̂(t), v̂(t+1)⟩| と
     ストロボドリフト d_str(k) = 1−|⟨v̂(k·den), v̂((k+1)·den)⟩|（den=N）。
     実軸（X,Y 個別）の保存はストロボでのみ意味を持つ——ロック後に 0 になるか。
 (c) 文脈量: 面内エネルギー比 h1(t) = (|a|²+|b|²)/‖Z‖²（=1−f）。
入力: N3_N40_long10000_20260905 の npz（SHA台帳照合）。新規走行なし。
出力: check_inplane_direction_conservation_v1.{json,csv,png}
"""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent
SWEEP = BASE.parent / "N3_N40_long10000_20260905"
RL = BASE.parent / "N3_N40_rotation_lock_analysis_20260905"
NS = [3, 4, 5, 6, 10, 40]

ledger = {}
with open(SWEEP / "SHA256SUMS.txt") as fh:
    for line in fh:
        parts = line.split()
        if len(parts) == 2:
            ledger[parts[1]] = parts[0]

rl = json.load(open(RL / 'rotation_lock_analysis_v1.json'))
onset = {r[0]: int(r[2]) for r in rl['table']}

E_PLUS = np.array([1.0, 1.0j]) / np.sqrt(2.0)

out = {}
rows = []
fig, axs = plt.subplots(2, 3, figsize=(16, 8))
for ax, N in zip(axs.ravel(), NS):
    rel = f'results/hm_N{N}_den_{N}_states_10000.npz'
    h = hashlib.sha256(open(SWEEP / rel, 'rb').read()).hexdigest()
    assert ledger[rel] == h, f'INPUT GATE FAIL: {rel}'
    Z = np.asarray(np.load(SWEEP / rel)['Z'], dtype=np.complex128)
    T, M = Z.shape
    rp, ip = Z[0].real, Z[0].imag
    p = rp / np.linalg.norm(rp)
    q = ip - (ip @ p) * p
    q = q / np.linalg.norm(q)
    a = Z @ p
    b = Z @ q
    v = np.stack([a, b], axis=1)
    vn = np.linalg.norm(v, axis=1)
    vhat = v / vn[:, None]
    h1 = (np.abs(a) ** 2 + np.abs(b) ** 2) / np.sum(np.abs(Z) ** 2, axis=1)
    c_plus = np.abs(vhat @ E_PLUS.conj())
    d_cont = 1.0 - np.abs(np.einsum('ij,ij->i', vhat[:-1].conj(), vhat[1:]))
    den = N
    vs = vhat[::den]
    d_str = 1.0 - np.abs(np.einsum('ij,ij->i', vs[:-1].conj(), vs[1:]))
    sc = onset[N]
    ramp_end = max(sc - 10, 1)
    res = {'onset': sc,
           'c_plus': {'step0': float(c_plus[0]),
                      'ramp_min_(0..onset-10)': float(c_plus[:ramp_end].min()),
                      'global_min': float(c_plus.min()),
                      'argmin_step': int(np.argmin(c_plus)),
                      'final_mean_last1000': float(c_plus[-1000:].mean())},
           'cp1_drift_per_step': {'ramp_max': float(d_cont[:ramp_end].max()),
                                  'final_mean_last1000': float(d_cont[-1000:].mean())},
           'cp1_strobe_drift': {'first10_mean': float(d_str[:10].mean()),
                                'final_mean_last100': float(d_str[-100:].mean()),
                                'final_max_last100': float(d_str[-100:].max())},
           'inplane_fraction_h1': {'step0': float(h1[0]), 'final_mean_last1000': float(h1[-1000:].mean())}}
    out[str(N)] = res
    rows.append([N, sc, f"{c_plus[0]:.10f}", f"{res['c_plus']['ramp_min_(0..onset-10)']:.10f}",
                 f"{res['c_plus']['global_min']:.6f}", res['c_plus']['argmin_step'],
                 f"{res['c_plus']['final_mean_last1000']:.6f}",
                 f"{res['cp1_strobe_drift']['final_mean_last100']:.3e}",
                 f"{res['inplane_fraction_h1']['final_mean_last1000']:.6f}"])
    ax.plot(c_plus, lw=0.7, color='tab:green', label='c+ (circular purity)')
    ax.plot(np.arange(0, T, den)[:-1] if len(d_str) == len(vs) - 1 else np.arange(len(d_str)) * den,
            1.0 - d_str, lw=0.7, color='tab:purple', alpha=.7, label='1 - strobe CP1 drift')
    ax.axvline(sc, color='k', ls=':', lw=0.8)
    ax.set_title(f'N={N}')
    ax.set_xscale('log')
    ax.set_ylim(-0.02, 1.02)
    ax.grid(alpha=.3)
    ax.legend(fontsize=7)
    print(f"N={N}: c+ step0={c_plus[0]:.8f} 増幅期min={res['c_plus']['ramp_min_(0..onset-10)']:.8f} "
          f"全域min={c_plus.min():.4f}(step {np.argmin(c_plus)}) 終盤={res['c_plus']['final_mean_last1000']:.6f} | "
          f"ストロボCP¹ドリフト 終盤={res['cp1_strobe_drift']['final_mean_last100']:.2e} | "
          f"面内比 終盤={res['inplane_fraction_h1']['final_mean_last1000']:.4f}", flush=True)

HEADER = ['N', 'onset', 'c_plus_step0', 'c_plus_ramp_min', 'c_plus_global_min', 'argmin_step',
          'c_plus_final_mean', 'strobe_cp1_drift_final', 'h1_final_mean']
with open(BASE / 'check_inplane_direction_conservation_v1.csv', 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(HEADER); w.writerows(rows)
with open(BASE / 'check_inplane_direction_conservation_v1.json', 'w') as fh:
    json.dump(out, fh, indent=2)
fig.suptitle('In-plane direction conservation: circular purity c+ and strobe CP1 stability; dotted = onset', y=.995)
fig.tight_layout()
fig.savefig(BASE / 'check_inplane_direction_conservation_v1.png', dpi=140)
print('ALL DONE')
