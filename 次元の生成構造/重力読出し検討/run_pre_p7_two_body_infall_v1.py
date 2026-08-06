#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""重力読出し P7: 二塊の関係落下実験——万有引力は万能相互作用に入っているか

原理（事前記録・木原の枠組み）: 重力はM体でも究極2体（自分×残りの合成重心）。
落下の三方向分解: R方向=不可（反力=質量）／時間方向=局所時計レート変化／
空間方向=測地線落下。読むのは関係量2つ（相対測定=共通モード除去=非循環）:
  (a) 相対距離 d(t)（重心間・厳密辺長）——縮めば空間落下＝引力の直接実証
  (b) 領域時計レート比——時間方向落下（重力時間遅れ対応・γ橋）
力の項は一切書かない。万能頂点＋海だけで d が縮むかを見る。

対照: 単一パケット（相手なし）の位置安定性。
使い方: python3 run_pre_p7_two_body_infall_v1.py
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
        spec = importlib.util.spec_from_file_location("g3d_p7", MB / "run_genesis_3d_demo_v2.py")
        g3d = importlib.util.module_from_spec(spec); sys.modules[spec.name] = g3d
        spec.loader.exec_module(g3d)
        spec2 = importlib.util.spec_from_file_location("g3z_p7", MB / "run_genesis_v3_register_local_v1.py")
        g3 = importlib.util.module_from_spec(spec2); sys.modules[spec2.name] = g3
        spec2.loader.exec_module(g3)
    finally:
        os.chdir(_cwd)
    abl3, Engine3D, D = g3d.abl, g3d.Engine3D, g3d.D
    DELTA, N_GRAPH, m = g3d.DELTA, g3d.N_GRAPH, g3d.m
    odd_P3 = g3d.odd_P3

    _, v3, _, _, _, _, _, Z0c, wp0 = abl3.build_init(N_GRAPH, False)
    Z0c = Z0c / np.linalg.norm(Z0c)
    SIG = 0.8
    n2 = np.arange(D[1]); n3 = np.arange(D[2])
    def kprof(dy, dz):
        g2 = np.exp(-0.5 * ((np.minimum(n2, D[1] - n2)) / SIG) ** 2)
        g3_ = np.exp(-0.5 * ((np.minimum(n3, D[2] - n3)) / SIG) ** 2)
        bump = np.outer(g2, g3_)
        Fk = np.fft.fftn(bump)
        k2g, k3g = np.meshgrid(np.arange(D[1]), np.arange(D[2]), indexing="ij")
        Fk = Fk * np.exp(-2j * np.pi * (dy * k2g / D[1] + dz * k3g / D[2]))
        prof = np.zeros(D, complex)
        for k1 in range(D[0]):
            if k1 % 2 == 1:
                prof[k1] = Fk
        return prof / np.linalg.norm(prof)
    seeds = [g3.zero_closure_state(m, np.random.default_rng(98000 + i)) for i in range(2)]
    AMP = 5.0 * DELTA   # 重め

    def centroid_axis(P, w):
        Pw = P * w
        tot = Pw.sum()
        cy = np.angle(np.sum(Pw.sum(axis=(0, 2)) * np.exp(2j * np.pi * n2 / D[1]))) * D[1] / (2 * np.pi) % D[1]
        cz = np.angle(np.sum(Pw.sum(axis=(0, 1)) * np.exp(2j * np.pi * n3 / D[2]))) * D[2] / (2 * np.pi) % D[2]
        return float(cy), float(cz), float(tot)

    def region_w(cy, cz, hw=1.9):
        w = np.zeros(D)
        for a in range(D[1]):
            for b in range(D[2]):
                da = min(abs(a - cy), D[1] - abs(a - cy))
                db = min(abs(b - cz), D[2] - abs(b - cz))
                if da <= hw and db <= hw:
                    w[:, a, b] = 1.0
        return w

    def torus_d(p1, p2):
        dy = min(abs(p1[0] - p2[0]), D[1] - abs(p1[0] - p2[0]))
        dz = min(abs(p1[1] - p2[1]), D[2] - abs(p1[1] - p2[1]))
        return float(np.hypot(dy, dz))

    def run(offsets, label, T=3000, every=300):
        C3_0 = np.zeros((m,) + D, complex); C3_0[:, 2, 0, 0] = Z0c
        for i, (dy, dz) in enumerate(offsets):
            C3_0 = C3_0 + AMP * kprof(dy, dz)[None] * seeds[i][:, None, None, None]
        eng = Engine3D(N_GRAPH, C3_0, wp0, (0.0, 0.0, 0.0), vertex_on=True)
        rows = []
        # 領域時計: 各領域の状態位相前進（比のみ使用）
        phases = [[], []]
        for t in range(T + 1):
            if t % every == 0:
                P = odd_P3(eng.C3())
                cents = []
                for (dy, dz) in offsets:
                    w = region_w(dy, dz)
                    cy, cz, tot = centroid_axis(P, w)
                    cents.append((cy, cz, tot))
                if len(offsets) == 2:
                    d = torus_d(cents[0][:2], cents[1][:2])
                    rows.append({"t": t, "d": d,
                                  "c": [[round(c[0], 3), round(c[1], 3)] for c in cents],
                                  "P": [round(c[2], 4) for c in cents]})
                    print(f"  [{label}] t={t}: d={d:.4f} 重心={[(round(c[0],2),round(c[1],2)) for c in cents]}")
                else:
                    rows.append({"t": t, "c": [[round(cents[0][0], 3), round(cents[0][1], 3)]]})
                    print(f"  [{label}] t={t}: 重心=({cents[0][0]:.3f},{cents[0][1]:.3f})")
            if t < T:
                eng.step()
            # 領域時計位相（毎ステップ・後で比を取る）
            if len(offsets) == 2 and t % 10 == 0:
                C = eng.C3()
                for i, (dy, dz) in enumerate(offsets):
                    w = region_w(dy, dz)
                    zi = (C * w[None]).sum(axis=(1, 2, 3))
                    phases[i].append(np.angle(np.vdot(seeds[i], zi)))
        out = {"rows": rows}
        if len(offsets) == 2:
            ph0 = np.unwrap(np.array(phases[0])); ph1 = np.unwrap(np.array(phases[1]))
            r0 = np.polyfit(np.arange(len(ph0)), ph0, 1)[0]
            r1 = np.polyfit(np.arange(len(ph1)), ph1, 1)[0]
            out["clock_rates"] = [float(r0), float(r1)]
            out["clock_ratio"] = float(r0 / r1) if r1 != 0 else None
            print(f"  [{label}] 領域時計レート比 = {out['clock_ratio']:.6f}")
        return out

    print("== 対照: 単一塊 ==")
    ctrl = run([(0.0, 0.0)], "単一", T=3000, every=600)
    print("== 本試験: 二塊 d0=3 ==")
    two = run([(0.0, 0.0), (3.0, 0.0)], "二塊", T=3000, every=300)
    out = {"AMP_over_DELTA": 5.0, "control": ctrl, "two_body": two,
           "runtime_sec": time.time() - t0}
    (HERE / "pre_p7_two_body_infall_result_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s")

if __name__ == "__main__":
    main()
