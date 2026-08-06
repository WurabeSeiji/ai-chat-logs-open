#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""重力読出し P4: ゲージ・フリー距離汎関数（相関減衰距離）の構築と検分

原理（事前記録）: 方向の輸送（ゲージ）を捨て、距離だけの連鎖
  ゲージ不変スカラー → 距離 → 余弦定理 → 角度 → 内角和
で角度を得る（Regge の辺長のみ定式化と同じ回避）。

距離汎関数 D1（相関減衰）: 各パケット i の空間プロファイル
  ψ_i(cells,t) = ⟨seed_i(辺), C3(辺,cells)⟩
の正規化重なり G_ij = ⟨|⟨ψ_i,ψ_j⟩|⟩_t / √(P_i P_j)。モジュラスなので
局所・大域位相に不変＝ゲージ・フリー。

検分手順（事前固定）:
  1) 較正: 2パケットのレジスタ分離 r=1..4 で G(r) が単調減衰するか。
  2) 距離化: d(G) を較正曲線の逆写像（単調なら可能）で構成。
  3) 三角形（平坦対照）: 3パケット (0,0)/(L,0)/(0,L) の3つの G→d、
     三角不等式の成立と、余弦定理での内角和＝180°の検定。
  通れば主観側ゲージ不変距離連鎖が確立（G1）。通らなければ汎関数不採用の記録。

使い方: python3 run_pre_p4_gaugefree_distance_v1.py
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
        spec = importlib.util.spec_from_file_location("g3d_p4", MB / "run_genesis_3d_demo_v2.py")
        g3d = importlib.util.module_from_spec(spec); sys.modules[spec.name] = g3d
        spec.loader.exec_module(g3d)
        spec2 = importlib.util.spec_from_file_location("g3z_p4", MB / "run_genesis_v3_register_local_v1.py")
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
    seeds = [g3.zero_closure_state(m, np.random.default_rng(98000 + i)) for i in range(3)]
    def dseed(i, dy, dz):
        grad = np.exp(2j * np.pi * (dy * k2g / D[1] + dz * k3g / D[2]))
        return prof[None] * seeds[i][:, None, None, None] * grad[None, None, :, :]

    def run_system(offsets, T=150):
        C3_0 = np.zeros((m,) + D, complex); C3_0[:, 2, 0, 0] = Z0c
        for i, (dy, dz) in enumerate(offsets):
            C3_0 = C3_0 + DELTA * dseed(i, dy, dz)
        eng = Engine3D(N_GRAPH, C3_0, wp0, (0.0, 0.0, 0.0), vertex_on=True)
        # 参照系（パケットなし・同一初期化）——差分場で背景汚染を除去
        C3_r = np.zeros((m,) + D, complex); C3_r[:, 2, 0, 0] = Z0c
        ref = Engine3D(N_GRAPH, C3_r, wp0, (0.0, 0.0, 0.0), vertex_on=True)
        for _ in range(50):
            eng.step(); ref.step()
        Gsum = np.zeros((len(offsets), len(offsets))); Pw = np.zeros(len(offsets)); nS = 0
        for t in range(T):
            eng.step(); ref.step()
            if t % 5 == 0:
                C = eng.C3() - ref.C3()
                oddmask = (np.arange(D[0]) % 2 == 1)
                psis = []
                for i in range(len(offsets)):
                    psi_k = np.tensordot(np.conj(seeds[i]), C, axes=(0, 0))
                    psi_k = psi_k * oddmask[:, None, None]
                    psis.append(np.fft.ifftn(psi_k))   # 位置領域（odd_P3と同一規約）
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
        return Gn

    # 1) 較正: 分離 r 掃引（2パケット）
    print("== 較正: G(r) 単調性 ==")
    cal = []
    for r in (1.0, 1.5, 2.0, 2.5, 3.0, 4.0):
        Gn = run_system([(0.0, 0.0), (r, 0.0)])
        g = float(Gn[0, 1])
        cal.append((r, g))
        print(f"  r={r}: G={g:.5f}")
    gs = [g for _, g in cal]
    monotone = all(gs[i] >= gs[i + 1] - 1e-6 for i in range(len(gs) - 1))
    print(f"  単調減衰 = {monotone}")

    out = {"calibration": cal, "monotone": monotone}
    if monotone:
        # 2) 距離化: 較正曲線の線形補間逆写像
        rs = np.array([r for r, _ in cal]); gv = np.array(gs)
        def dist_of(G):
            return float(np.interp(G, gv[::-1], rs[::-1]))
        # 3) 三角形（平坦対照）
        print("== 三角形（平坦対照・L=2.5）==")
        L = 2.5
        Gtri = run_system([(0.0, 0.0), (L, 0.0), (0.0, L)])
        dAB = dist_of(float(Gtri[0, 1])); dAC = dist_of(float(Gtri[0, 2])); dBC = dist_of(float(Gtri[1, 2]))
        print(f"  G: AB={Gtri[0,1]:.5f} AC={Gtri[0,2]:.5f} BC={Gtri[1,2]:.5f}")
        print(f"  d: AB={dAB:.4f} AC={dAC:.4f} BC={dBC:.4f}  （真値: {L}, {L}, {L*np.sqrt(2):.4f}）")
        tri_ok = dAB + dAC > dBC and dAB + dBC > dAC and dAC + dBC > dAB
        import math
        if tri_ok:
            A = math.acos(max(-1, min(1, (dAB**2 + dAC**2 - dBC**2) / (2 * dAB * dAC))))
            B = math.acos(max(-1, min(1, (dAB**2 + dBC**2 - dAC**2) / (2 * dAB * dBC))))
            Cang = math.acos(max(-1, min(1, (dAC**2 + dBC**2 - dAB**2) / (2 * dAC * dBC))))
            s = math.degrees(A + B + Cang)
            print(f"  三角不等式 OK・内角和 = {s:.4f}°（平坦対照: 180°）")
            out["triangle"] = {"G": [float(Gtri[0,1]), float(Gtri[0,2]), float(Gtri[1,2])],
                                "d": [dAB, dAC, dBC], "angle_sum_deg": s}
        else:
            print("  三角不等式 破れ——距離汎関数D1は現形では不採用（記録）")
            out["triangle"] = {"tri_inequality": False}
    out["runtime_sec"] = time.time() - t0
    (HERE / "pre_p4_gaugefree_distance_result_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s")

if __name__ == "__main__":
    main()
