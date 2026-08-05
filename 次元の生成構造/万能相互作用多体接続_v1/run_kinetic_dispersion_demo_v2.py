#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""運動学の実証 v2——計器修正版（反証と修正の記録対象）

v1 の計器欠陥: 重心を基本巻きモーメント Σ P(n) e^{2πin/N} で読んでいた。
奇数倍音のみのパケットは |W(n)|² が周期 N/2 を持つため、このモーメントは
恒等的にゼロ（実測: 数値ノイズの角度を読んでいた）。v1 の位置系列
[6.9, 10.0, 5.37, 8.0] は測定不能量のノイズである。

v2 の修正: 周期 N/2 の信号に対して well-defined な巻き数2のモーメント
    z₂ = Σ P(n) e^{2πi·2n/N},  位置 = (arg z₂)·N/(4π) mod N/2
を用いる。位置は mod N/2 で定義される（奇数倍音パケットの物理的な位置分解能）。

判定（実行前固定）:
    K1 v1計器の恒等零の確認: |z₁|/ΣP < 1e-3（ゼロベクトルのモーメント）。
    K2 並進: v2位置が予言 (x₀+0.05t) mod 8 と全記録時刻で ±0.2 セル以内。
    K3 形不変: PR（巻き数2定義の実効幅）が全時刻で初期値±10%。

使い方: python3 run_kinetic_dispersion_demo_v2.py
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec3 = importlib.util.spec_from_file_location("s3kin2", HERE / "run_genesis_v3_register_local_v1.py")
g3 = importlib.util.module_from_spec(spec3)
sys.modules[spec3.name] = g3
spec3.loader.exec_module(g3)
abl, V2 = g3.abl, g3.V2

N_GRAPH, NREG, DELTA = 12, 16, 3e-2
m = N_GRAPH * (N_GRAPH - 1) // 2
_, v, _, _, _, _, _, Z0c, wp0 = abl.build_init(N_GRAPH, False)
Z0c = Z0c / np.linalg.norm(Z0c)
seed_edge = g3.zero_closure_state(m, np.random.default_rng(98000))
odd_ks = [k for k in range(NREG) if k % 2 == 1]


class KineticEngine(V2):
    def __init__(self, n, C0, wp, omega1, **kw):
        super().__init__(n, C0, wp, **kw)
        self.disp = np.exp(1j * np.arange(self.nreg) * omega1)[None, :]

    def step(self):
        super().step()
        self.C = self.C * self.disp


def centroid_v2(C):
    """奇数倍音パケットの位置（mod N/2）: 巻き数2モーメント。v1恒等零の診断値つき。"""
    N = C.shape[1]
    ks = np.arange(N)
    Codd = C * ((ks % 2 == 1)[None, :])
    W = np.fft.ifft(Codd, axis=1) * N
    P = np.sum(np.abs(W) ** 2, axis=0)
    z1 = np.sum(P * np.exp(2j * np.pi * ks / N)) / P.sum()      # v1計器（恒等零のはず）
    z2 = np.sum(P * np.exp(2j * np.pi * 2 * ks / N)) / P.sum()  # v2計器
    x = float((np.angle(z2) * N / (4 * np.pi)) % (N / 2))
    pr = float(P.sum() ** 2 / np.sum(P ** 2) / 2)               # 周期N/2内の実効幅
    return x, pr, abs(z1)


prof = np.zeros(NREG, complex)
for k in odd_ks:
    prof[k] = 1.0 / np.sqrt(len(odd_ks))
C0 = np.zeros((m, NREG), complex)
C0[:, 2] = Z0c
for k in range(NREG):
    if abs(prof[k]) > 0:
        C0[:, k] += DELTA * prof[k] * seed_edge

omega1 = 2 * np.pi / NREG * 0.05
eng = KineticEngine(N_GRAPH, C0, wp0, omega1, vertex_on=True)
x0, pr0, z1_0 = centroid_v2(eng.C)
print(f"初期: 位置(mod 8)={x0:.3f} 実効幅PR={pr0:.2f}  v1計器の恒等零診断 |z1|={z1_0:.2e}")
rows = []
k1_ok, k2_ok, k3_ok = z1_0 < 1e-3, True, True
for t in range(400):
    eng.step()
    if (t + 1) % 50 == 0:
        x, pr, z1a = centroid_v2(eng.C)
        pred = (x0 + 0.05 * (t + 1)) % 8
        dev = min(abs(x - pred), 8 - abs(x - pred))
        rows.append({"t": t + 1, "x": round(x, 3), "pred": round(pred, 3),
                     "dev": round(dev, 3), "PR": round(pr, 2), "z1_abs": z1a})
        k1_ok &= z1a < 1e-3
        k2_ok &= dev <= 0.2
        k3_ok &= abs(pr - pr0) <= 0.1 * pr0
        print(f"  t={t+1:3d}: 位置={x:.3f} 予言={pred:.3f} 偏差={dev:.3f} PR={pr:.2f} |z1|={z1a:.1e}")
print(f"判定: K1(v1計器恒等零)={k1_ok}  K2(並進±0.2)={k2_ok}  K3(形不変±10%)={k3_ok}")
json.dump({"instrument_correction": "v1 fundamental-winding moment vanishes identically for odd-harmonic packets; v2 uses winding-2 moment, position defined mod N/2",
           "omega1_cells_per_step": 0.05, "x0": x0, "PR0": pr0,
           "rows": rows, "K1": bool(k1_ok), "K2": bool(k2_ok), "K3": bool(k3_ok),
           "all_pass": bool(k1_ok and k2_ok and k3_ok)},
          open(HERE / "kinetic_dispersion_demo_v2.json", "w"), indent=1)
print("saved")
