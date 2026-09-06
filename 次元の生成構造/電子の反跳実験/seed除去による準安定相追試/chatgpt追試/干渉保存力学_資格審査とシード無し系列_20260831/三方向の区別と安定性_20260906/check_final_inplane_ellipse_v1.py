#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""三方向の区別と安定性 (3): 終状態の面内楕円の直接測定（読み出しのみ）。

検査2の c₊(終盤) が N=3,4 で同一値 0.994936 だった。円偏光純度の低下は
 (i) 振幅比 |b|/|a| が 1 からずれる（楕円化）か (ii) 位相差 arg(b/a) が 90° から
ずれるかの2経路がある。どちらかを直接測る。
候補仮説（検定対象）: 位相差 90° 維持のまま振幅比 r で楕円化した場合
 c₊ = (1+r)/√(2(1+r²)) であり、c₊=0.994936 ⇒ r=√(2/3) または √(3/2)。
これは面法線の傾き arccos√(2/3)=35.264°（検査1）と同じ √(2/3) 系の値。
入力: N3_N40_long10000_20260905 の npz（SHA台帳照合）。最終 1000 歩の平均で報告。
出力: check_final_inplane_ellipse_v1.json
"""
import hashlib
import json
from pathlib import Path

import numpy as np

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
for N in NS:
    rel = f'results/hm_N{N}_den_{N}_states_10000.npz'
    h = hashlib.sha256(open(SWEEP / rel, 'rb').read()).hexdigest()
    assert ledger[rel] == h, f'INPUT GATE FAIL: {rel}'
    Z = np.asarray(np.load(SWEEP / rel)['Z'], dtype=np.complex128)
    rp, ip = Z[0].real, Z[0].imag
    p = rp / np.linalg.norm(rp)
    q = ip - (ip @ p) * p
    q = q / np.linalg.norm(q)
    a = Z[-1000:] @ p
    b = Z[-1000:] @ q
    r_amp = np.abs(b) / np.abs(a)
    dphase = np.degrees(np.angle(b / a))
    r_mean = float(r_amp.mean())
    res = {'amp_ratio_b_over_a_mean': r_mean,
           'amp_ratio_std': float(r_amp.std()),
           'phase_diff_deg_mean': float(dphase.mean()),
           'phase_diff_deg_std': float(dphase.std()),
           'r_vs_sqrt2over3': float(r_mean / np.sqrt(2.0 / 3.0)),
           'r_vs_sqrt3over2': float(r_mean / np.sqrt(3.0 / 2.0)),
           'c_plus_predicted_from_r': float((1 + r_mean) / np.sqrt(2 * (1 + r_mean ** 2)))}
    out[str(N)] = res
    print(f"N={N}: |b/a|={r_mean:.8f}±{r_amp.std():.1e} 位相差={dphase.mean():.4f}°±{dphase.std():.1e} | "
          f"r/√(2/3)={res['r_vs_sqrt2over3']:.8f} r/√(3/2)={res['r_vs_sqrt3over2']:.8f} | "
          f"c₊予測={res['c_plus_predicted_from_r']:.8f}", flush=True)

with open(BASE / 'check_final_inplane_ellipse_v1.json', 'w') as fh:
    json.dump(out, fh, indent=2)
print('ALL DONE')
