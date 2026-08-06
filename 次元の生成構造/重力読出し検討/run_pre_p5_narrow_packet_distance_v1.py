#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""重力読出し P5: 狭パケット構成によるD1距離の再較正

P4の原因（パケット幅広・テール飽和）への対処: k2,k3方向を位置空間ガウス包絡
（幅σセル）のFTで構成した狭パケットに置換し、D1（位置領域相関距離）を再較正。

手順（事前固定）:
  0) 診断: 旧構成と新構成の差分場パケットの位置プロファイル（周辺分布・
     中心占有率）を比較。
  1) 較正: r=1..4 で G(r) 単調減衰か。
  2) 通過なら距離化→三角形（平坦対照・内角和180°）関門。
使い方: python3 run_pre_p5_narrow_packet_distance_v1.py
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
        spec = importlib.util.spec_from_file_location("g3d_p5", MB / "run_genesis_3d_demo_v2.py")
        g3d = importlib.util.module_from_spec(spec); sys.modules[spec.name] = g3d
        spec.loader.exec_module(g3d)
        spec2 = importlib.util.spec_from_file_location("g3z_p5", MB / "run_genesis_v3_register_local_v1.py")
        g3 = importlib.util.module_from_spec(spec2); sys.modules[spec2.name] = g3
        spec2.loader.exec_module(g3)
    finally:
        os.chdir(_cwd)
    abl3, Engine3D, D = g3d.abl, g3d.Engine3D, g3d.D
    DELTA, N_GRAPH, m = g3d.DELTA, g3d.N_GRAPH, g3d.m

    _, v3, _, _, _, _, _, Z0c, wp0 = abl3.build_init(N_GRAPH, False)
    Z0c = Z0c / np.linalg.norm(Z0c)
    seeds = [g3.zero_closure_state(m, np.random.default_rng(98000 + i)) for i in range(3)]

    # 狭パケット: 位置空間ガウス（幅σ）→ k空間プロファイル（k2,k3）。k1は奇数のみ。
    SIG = 0.8
    n2 = np.arange(D[1]); n3 = np.arange(D[2])
    def kprof(dy, dz):
        g2 = np.exp(-0.5 * ((np.minimum(n2, D[1] - n2)) / SIG) ** 2)
        g3_ = np.exp(-0.5 * ((np.minimum(n3, D[2] - n3)) / SIG) ** 2)
        bump = np.outer(np.roll(g2, 0), np.roll(g3_, 0))
        Fk = np.fft.fftn(bump)                       # (k2,k3)
        k2g, k3g = np.meshgrid(np.arange(D[1]), np.arange(D[2]), indexing="ij")
        Fk = Fk * np.exp(-2j * np.pi * (dy * k2g / D[1] + dz * k3g / D[2]))
        prof = np.zeros(D, complex)
        for k1 in range(D[0]):
            if k1 % 2 == 1:
                prof[k1] = Fk
        return prof / np.linalg.norm(prof)

    def dseed(i, dy, dz):
        return kprof(dy, dz)[None] * seeds[i][:, None, None, None]

    oddmask = (np.arange(D[0]) % 2 == 1)
    def run_system(offsets, T=150, diag=False):
        C3_0 = np.zeros((m,) + D, complex); C3_0[:, 2, 0, 0] = Z0c
        for i, (dy, dz) in enumerate(offsets):
            C3_0 = C3_0 + DELTA * dseed(i, dy, dz)
        eng = Engine3D(N_GRAPH, C3_0, wp0, (0.0, 0.0, 0.0), vertex_on=True)
        C3_r = np.zeros((m,) + D, complex); C3_r[:, 2, 0, 0] = Z0c
        ref = Engine3D(N_GRAPH, C3_r, wp0, (0.0, 0.0, 0.0), vertex_on=True)
        for _ in range(50):
            eng.step(); ref.step()
        Gsum = np.zeros((len(offsets), len(offsets))); Pw = np.zeros(len(offsets)); nS = 0
        diag_marg = None
        for t in range(T):
            eng.step(); ref.step()
            if t % 5 == 0:
                C = eng.C3() - ref.C3()
                psis = []
                for i in range(len(offsets)):
                    pk = np.tensordot(np.conj(seeds[i]), C, axes=(0, 0)) * oddmask[:, None, None]
                    psis.append(np.fft.ifftn(pk))
                if diag and diag_marg is None:
                    P0 = np.abs(psis[0]) ** 2
                    diag_marg = (P0.sum(axis=(0, 2)) / P0.sum(), P0.sum(axis=(0, 1)) / P0.sum())
                for i in range(len(offsets)):
                    Pw[i] += float(np.sum(np.abs(psis[i]) ** 2))
                    for j in range(len(offsets)):
                        if j > i:
                            Gsum[i, j] += abs(np.vdot(psis[i], psis[j]))
                nS += 1
        Pw /= nS; Gsum /= nS
        Gn = np.zeros_like(Gsum)
        for i in range(len(offsets)):
            for j in range(len(offsets)):
                if j > i:
                    Gn[i, j] = Gsum[i, j] / np.sqrt(Pw[i] * Pw[j])
        return (Gn, diag_marg) if diag else Gn

    # 0) 診断
    _, marg = run_system([(0.0, 0.0), (4.0, 0.0)], T=60, diag=True)
    print(f"[診断] 狭パケット位置周辺分布 n2: {[round(float(x),3) for x in marg[0]]}")
    print(f"        中心占有率 = {float(max(marg[0])):.3f}")

    # 1) 較正
    print("== 較正: G(r) ==")
    cal = []
    for r in (1.0, 1.5, 2.0, 2.5, 3.0, 4.0):
        Gn = run_system([(0.0, 0.0), (r, 0.0)])
        g = float(Gn[0, 1]); cal.append((r, g))
        print(f"  r={r}: G={g:.5f}")
    gs = [g for _, g in cal]
    monotone = all(gs[i] >= gs[i + 1] - 0.04 for i in range(len(gs) - 1))
    print(f"  単調減衰（許容0.04）= {monotone}")
    out = {"SIG": SIG, "center_occupancy": float(max(marg[0])), "calibration": cal,
           "monotone": monotone}

    if monotone:
        # 床補正＋単調化（cummin）した較正曲線で距離化
        rs = np.array([r for r, _ in cal]); gv = np.minimum.accumulate(np.array(gs))
        floor = gv[-1] * 0.98
        gcorr = (gv - floor) / (1 - floor)
        def dist_of(G):
            Gc = max((G - floor) / (1 - floor), 1e-4)
            return float(np.interp(Gc, gcorr[::-1], rs[::-1]))
        L = 2.5
        Gtri = run_system([(0.0, 0.0), (L, 0.0), (0.0, L)])
        dAB = dist_of(float(Gtri[0, 1])); dAC = dist_of(float(Gtri[0, 2])); dBC = dist_of(float(Gtri[1, 2]))
        print(f"== 三角形（L={L}）== d: AB={dAB:.4f} AC={dAC:.4f} BC={dBC:.4f} "
              f"（真値 {L}, {L}, {L*np.sqrt(2):.4f}）")
        import math
        tri_ok = dAB + dAC > dBC and dAB + dBC > dAC and dAC + dBC > dAB
        if tri_ok:
            A = math.acos(max(-1, min(1, (dAB**2 + dAC**2 - dBC**2) / (2 * dAB * dAC))))
            B = math.acos(max(-1, min(1, (dAB**2 + dBC**2 - dAC**2) / (2 * dAB * dBC))))
            Cc = math.acos(max(-1, min(1, (dAC**2 + dBC**2 - dAB**2) / (2 * dAC * dBC))))
            s = math.degrees(A + B + Cc)
            print(f"  三角不等式OK・内角和 = {s:.4f}°（平坦対照180°）")
            out["triangle"] = {"d": [dAB, dAC, dBC], "angle_sum_deg": s}
        else:
            print("  三角不等式 破れ")
            out["triangle"] = {"tri_inequality": False}
    out["runtime_sec"] = time.time() - t0
    (HERE / "pre_p5_narrow_packet_result_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s")

if __name__ == "__main__":
    main()
