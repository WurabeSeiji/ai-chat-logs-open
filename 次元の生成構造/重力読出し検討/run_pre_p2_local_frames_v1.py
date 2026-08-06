#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""重力読出し P2: 点ごとの局所枠とホロノミー（予備実験）

原理（事前記録）: Gauss三角測量は大域座標を要しない——各頂点の局所枠で
視線間の角を測れば内角和が出る。必要なのは点ごとの局所チャートと
枠間の遷移（接続）。枠をループ A→B→C→A で一周させた正味回転角＝
ホロノミー＝三角形を貫く曲率束（調査ノート§4・非循環系統）。

方法: P1-Bと同一の静止3パケット（レジスタ格子・(0,0)/(3,0)/(0,3)）。
各パケット領域 i の辺状態時系列 z_i(t)∈C^m から局所SVDで局所平面
（d1_i, d2_i）を構成。枠間の 2×2 重なり行列の回転部から遷移角 θ_ij を
取り、ホロノミー H = θ_AB + θ_BC + θ_CA を評価。

判定（記述的・事前固定）:
  各領域で局所平面が明確（sv比で確認）か／θ_ij の窓間安定性／
  H の値（平坦なら0。系統的非零なら曲率信号候補——ただし枠構成の
  人工物との判別は三角形サイズスケーリング（次段）まで保留）。

使い方: python3 run_pre_p2_local_frames_v1.py
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
        spec = importlib.util.spec_from_file_location("g3d_p2", MB / "run_genesis_3d_demo_v2.py")
        g3d = importlib.util.module_from_spec(spec); sys.modules[spec.name] = g3d
        spec.loader.exec_module(g3d)
        spec2 = importlib.util.spec_from_file_location("g3z_p2", MB / "run_genesis_v3_register_local_v1.py")
        g3 = importlib.util.module_from_spec(spec2); sys.modules[spec2.name] = g3
        spec2.loader.exec_module(g3)
    finally:
        os.chdir(_cwd)
    abl3, Engine3D, D = g3d.abl, g3d.Engine3D, g3d.D
    DELTA, N_GRAPH, m = g3d.DELTA, g3d.N_GRAPH, g3d.m

    _, v3, _, _, _, _, _, Z0c, wp0 = abl3.build_init(N_GRAPH, False)
    Z0c = Z0c / np.linalg.norm(Z0c)
    prof = np.zeros(D)
    for k1 in range(D[0]):
        if k1 % 2 == 0: continue
        for k2 in range(D[1]):
            for k3 in range(D[2]):
                if (2 * k2) % D[1] == 0 and k2 != 0: continue
                if (2 * k3) % D[2] == 0 and k3 != 0: continue
                prof[k1, k2, k3] = 1.0
    prof /= np.linalg.norm(prof)
    k2g, k3g = np.meshgrid(np.arange(D[1]), np.arange(D[2]), indexing="ij")
    def displaced_seed(rng_seed, dy, dz):
        se = g3.zero_closure_state(m, np.random.default_rng(rng_seed))
        grad = np.exp(2j * np.pi * (dy * k2g / D[1] + dz * k3g / D[2]))
        return prof[None, :, :, :] * se[:, None, None, None] * grad[None, None, :, :]
    import os as _os
    L = float(_os.environ.get("TRI_L", "3.0"))
    offsets = [(0.0, 0.0), (L, 0.0), (0.0, L)]
    C3_0 = np.zeros((m,) + D, complex)
    C3_0[:, 2, 0, 0] = Z0c
    for i, (dy, dz) in enumerate(offsets):
        C3_0 = C3_0 + DELTA * displaced_seed(98000 + i, dy, dz)
    eng = Engine3D(N_GRAPH, C3_0, wp0, (0.0, 0.0, 0.0), vertex_on=True)

    # 領域重み（P1と同一の窓）
    weights = []
    for (dy, dz) in offsets:
        w = np.zeros(D)
        for a in range(D[1]):
            for b in range(D[2]):
                da = min(abs(a - dy), D[1] - abs(a - dy))
                db = min(abs(b - dz), D[2] - abs(b - dz))
                if da <= 1.5 and db <= 1.5:
                    w[:, a, b] = 1.0
        weights.append(w)

    # 整定50 → T=300 サンプル毎3
    for _ in range(50):
        eng.step()
    series = [[], [], []]
    for t in range(300):
        eng.step()
        if t % 3 == 0:
            C = eng.C3()
            for i, w in enumerate(weights):
                zi = (C * w[None, :, :, :]).sum(axis=(1, 2, 3))
                series[i].append(zi.copy())
    frames = []
    for i in range(3):
        S = np.array(series[i])
        X = np.hstack([S.real, S.imag]); Xc = X - X.mean(axis=0)
        _, sv, Vt = np.linalg.svd(Xc, full_matrices=False)
        frames.append((Vt[0], Vt[1], sv[:4] / sv[0]))
        print(f"領域{i} 局所平面: sv比 = {[round(float(x),3) for x in sv[:4]/sv[0]]}")

    def trans_angle(fi, fj):
        O = np.array([[fi[0] @ fj[0], fi[0] @ fj[1]],
                      [fi[1] @ fj[0], fi[1] @ fj[1]]])
        theta = float(np.arctan2(O[0, 1] - O[1, 0], O[0, 0] + O[1, 1]))
        fid = float(np.linalg.norm(O))
        return theta, fid

    tAB, fAB = trans_angle(frames[0], frames[1])
    tBC, fBC = trans_angle(frames[1], frames[2])
    tCA, fCA = trans_angle(frames[2], frames[0])
    H = tAB + tBC + tCA
    Hn = (H + np.pi) % (2 * np.pi) - np.pi
    print(f"遷移角: θ_AB={np.degrees(tAB):+.4f}° θ_BC={np.degrees(tBC):+.4f}° "
          f"θ_CA={np.degrees(tCA):+.4f}°  （重なり強度 {fAB:.3f}/{fBC:.3f}/{fCA:.3f}）")
    print(f"ホロノミー H = {np.degrees(Hn):+.6f}°  （平坦なら0）")
    out = {"L": L, "sv_ratios": [[float(x) for x in f[2]] for f in frames],
           "theta_deg": {"AB": float(np.degrees(tAB)), "BC": float(np.degrees(tBC)),
                          "CA": float(np.degrees(tCA))},
           "overlap": [fAB, fBC, fCA],
           "holonomy_deg": float(np.degrees(Hn)),
           "runtime_sec": time.time() - t0}
    (HERE / ("pre_p2_local_frames_result_L%s_v1.json" % L)).write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s")

if __name__ == "__main__":
    main()
