#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""P7b: 時間方向の落下——局所時計レートの距離依存（重力時間遅れ対応）

P7の教訓: 静止塊はロックされ空間落下しない（運動=自前時計の特権）。
先に現れるのは時間方向の落下=局所時計レート変化。P7の時計計はseed依存で不良
→ 自己参照時計（領域状態の自己位相前進 arg⟨z_i(t+Δ), z_i(t)⟩）に修正。

測定（事前固定）: 二塊の分離 d0 ∈ {2,3,4} で、各領域の自己参照時計レートと
遠方参照領域（塊から最遠点）のレートの比 ρ_i(d0) = rate_i / rate_far。
判定 H_grav_t: ρ が d0 に系統依存（近いほど遅い/速いの単調則）→
時間方向落下（重力時間遅れ対応）の初検出。対称性検査: 同質量なら ρ_1=ρ_2。
使い方: python3 run_pre_p7b_clock_gravity_v1.py
"""
from __future__ import annotations
import importlib.util, json, os, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
MB = HERE.parent / "万能相互作用多体接続_v1"

def main():
    t0 = time.time()
    _cwd = os.getcwd(); os.chdir(MB)
    try:
        spec = importlib.util.spec_from_file_location("g3d_p7b", MB / "run_genesis_3d_demo_v2.py")
        g3d = importlib.util.module_from_spec(spec); sys.modules[spec.name] = g3d
        spec.loader.exec_module(g3d)
        spec2 = importlib.util.spec_from_file_location("g3z_p7b", MB / "run_genesis_v3_register_local_v1.py")
        g3 = importlib.util.module_from_spec(spec2); sys.modules[spec2.name] = g3
        spec2.loader.exec_module(g3)
    finally:
        os.chdir(_cwd)
    abl3, Engine3D, D = g3d.abl, g3d.Engine3D, g3d.D
    DELTA, N_GRAPH, m = g3d.DELTA, g3d.N_GRAPH, g3d.m
    _, v3, _, _, _, _, _, Z0c, wp0 = abl3.build_init(N_GRAPH, False)
    Z0c = Z0c / np.linalg.norm(Z0c)
    SIG = 0.8
    n2 = np.arange(D[1]); n3 = np.arange(D[2])
    def kprof(dy, dz):
        g2 = np.exp(-0.5 * ((np.minimum(n2, D[1] - n2)) / SIG) ** 2)
        g3_ = np.exp(-0.5 * ((np.minimum(n3, D[2] - n3)) / SIG) ** 2)
        Fk = np.fft.fftn(np.outer(g2, g3_))
        k2g, k3g = np.meshgrid(np.arange(D[1]), np.arange(D[2]), indexing="ij")
        Fk = Fk * np.exp(-2j * np.pi * (dy * k2g / D[1] + dz * k3g / D[2]))
        prof = np.zeros(D, complex)
        for k1 in range(D[0]):
            if k1 % 2 == 1:
                prof[k1] = Fk
        return prof / np.linalg.norm(prof)
    seeds = [g3.zero_closure_state(m, np.random.default_rng(98000 + i)) for i in range(2)]
    AMP = 5.0 * DELTA

    def region_w(cy, cz, hw=1.4):
        w = np.zeros(D)
        for a in range(D[1]):
            for b in range(D[2]):
                da = min(abs(a - cy), D[1] - abs(a - cy))
                db = min(abs(b - cz), D[2] - abs(b - cz))
                if da <= hw and db <= hw:
                    w[:, a, b] = 1.0
        return w

    out = {"rows": []}
    print(f"{'d0':>4} {'rate1':>10} {'rate2':>10} {'rate_far':>10} {'ρ1':>9} {'ρ2':>9}")
    for d0 in (2.0, 3.0, 4.0):
        C3_0 = np.zeros((m,) + D, complex); C3_0[:, 2, 0, 0] = Z0c
        for i, (dy, dz) in enumerate([(0.0, 0.0), (d0, 0.0)]):
            C3_0 = C3_0 + AMP * kprof(dy, dz)[None] * seeds[i][:, None, None, None]
        eng = Engine3D(N_GRAPH, C3_0, wp0, (0.0, 0.0, 0.0), vertex_on=True)
        regions = [region_w(0.0, 0.0), region_w(d0, 0.0),
                   region_w((d0 / 2) % D[1], 4.0)]   # 3つ目=遠方参照（垂直方向遠地）
        for _ in range(200):
            eng.step()
        zs_prev = None
        adv = [[], [], []]
        for t in range(800):
            eng.step()
            if t % 4 == 0:
                C = eng.C3()
                zs = [(C * w[None]).sum(axis=(1, 2, 3)) for w in regions]
                if zs_prev is not None:
                    for i in range(3):
                        adv[i].append(np.angle(np.vdot(zs_prev[i], zs[i])))
                zs_prev = zs
        rates = [float(np.mean(a)) for a in adv]
        r1, r2, rf = rates
        rho1 = r1 / rf if rf != 0 else float("nan")
        rho2 = r2 / rf if rf != 0 else float("nan")
        print(f"{d0:>4} {r1:>10.6f} {r2:>10.6f} {rf:>10.6f} {rho1:>9.5f} {rho2:>9.5f}")
        out["rows"].append({"d0": d0, "rates": rates, "rho1": rho1, "rho2": rho2})
    out["runtime_sec"] = time.time() - t0
    (HERE / "pre_p7b_clock_gravity_result_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s")

if __name__ == "__main__":
    main()
