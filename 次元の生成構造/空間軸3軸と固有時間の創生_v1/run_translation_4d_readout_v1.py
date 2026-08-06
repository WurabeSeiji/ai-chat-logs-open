#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""4D読出し実験: xyz＋τ（R′時計）で読んだ並進の高密度サンプリング

先行論文の3Dデモ（run_genesis_3d_demo_v2.py・計器修正版）を高密度サンプリング
で再走行し、(x,y,z,τ) の4次元軌道を記録する。τ はエンジンのステップ＝R′時計の刻み。
サンプリングは3軸のシフト (0.05t, 0.03t, 0.02t) が同時に整数セルになる
t=100 の倍数で行う（分数セルシフト時は離散サンプリング効果で
パケットが見かけ上拡がり重心計が汚染されるため——1Dで実証済みのPR振動の3D版）。

対照検定（実行前固定）: t=100,200,300 の位置が保存済み
genesis_3d_demo_result_v2.json の値と一致すること。

使い方: python3 run_translation_4d_readout_v1.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
MB = HERE.parent / "万能相互作用多体接続_v1"

# 3Dデモv2の機構を再利用（v1原本が回収コード由来で cwd 相対のため、MBへ一時移動して読込）
import os
_cwd = os.getcwd()
os.chdir(MB)
try:
    spec = importlib.util.spec_from_file_location("g3dv2_4d", MB / "run_genesis_3d_demo_v2.py")
    g3d = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = g3d
    spec.loader.exec_module(g3d)
finally:
    os.chdir(_cwd)


def main() -> None:
    t0 = time.time()
    stored = json.loads((MB / "genesis_3d_demo_result_v2.json").read_text())

    # v2 main() と同一の初期化（シード同一）
    import importlib.util as iu
    spec_g3 = iu.spec_from_file_location("g3z_4d", MB / "run_genesis_v3_register_local_v1.py")
    g3 = iu.module_from_spec(spec_g3)
    sys.modules[spec_g3.name] = g3
    spec_g3.loader.exec_module(g3)
    abl, Engine3D, D = g3d.abl, g3d.Engine3D, g3d.D
    DELTA, VEL, N_GRAPH, m = g3d.DELTA, g3d.VEL, g3d.N_GRAPH, g3d.m
    odd_P3, centroid3_v2 = g3d.odd_P3, g3d.centroid3_v2

    seed_edge = g3.zero_closure_state(m, np.random.default_rng(98000))
    _, v, _, _, _, _, _, Z0c, wp0 = abl.build_init(N_GRAPH, False)
    Z0c = Z0c / np.linalg.norm(Z0c)
    C3_0 = np.zeros((m,) + D, complex)
    C3_0[:, 2, 0, 0] = Z0c
    prof = np.zeros(D)
    for k1 in range(D[0]):
        if k1 % 2 == 0:
            continue
        for k2 in range(D[1]):
            for k3 in range(D[2]):
                if (2 * k2) % D[1] == 0 and k2 != 0:
                    continue
                if (2 * k3) % D[2] == 0 and k3 != 0:
                    continue
                prof[k1, k2, k3] = 1.0
    prof /= np.linalg.norm(prof)
    C3_0 += DELTA * prof[None, :, :, :] * seed_edge[:, None, None, None]

    eng = Engine3D(N_GRAPH, C3_0, wp0, VEL, vertex_on=True)
    x0, pr0 = centroid3_v2(odd_P3(eng.C3()))
    mods = (4, 8, 8)
    rows = [{"t": 0, "x": [round(v_, 4) for v_ in x0], "PR": round(pr0, 2)}]
    T = 900
    EVERY = 100
    for t in range(T):
        eng.step()
        if (t + 1) % EVERY == 0:
            x, pr = centroid3_v2(odd_P3(eng.C3()))
            rows.append({"t": t + 1, "x": [round(v_, 4) for v_ in x], "PR": round(pr, 2)})
            print(f"  t={t+1:3d}: (x,y,z)=({x[0]:.3f},{x[1]:.3f},{x[2]:.3f}) PR={pr:.2f}")

    # 対照検定
    ok = True
    for srow in stored["rows"]:
        me = next(r for r in rows if r["t"] == srow["t"])
        dev = max(min(abs(me["x"][a] - srow["x"][a]),
                       mods[a] - abs(me["x"][a] - srow["x"][a])) for a in range(3))
        ok &= dev < 1e-3
    print(f"対照検定（保存済みv2の3点と一致）= {ok}")

    # 展開座標（予言に最も近い合同値を選ぶ・表示用）
    unwrapped = []
    for r in rows:
        u = []
        for a in range(3):
            pred = x0[a] + VEL[a] * r["t"]
            k = round((pred - r["x"][a]) / mods[a])
            u.append(round(r["x"][a] + k * mods[a], 4))
        unwrapped.append(u)

    out = {"VEL": list(VEL), "x0": [float(v_) for v_ in x0], "mods": list(mods),
           "EVERY": EVERY, "rows": rows, "unwrapped": unwrapped,
           "contrast_ok": bool(ok), "runtime_sec": time.time() - t0}
    (HERE / "translation_4d_readout_result_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
