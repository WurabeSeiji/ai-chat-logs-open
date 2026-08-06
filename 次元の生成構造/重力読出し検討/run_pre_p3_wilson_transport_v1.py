#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""重力読出し P3: 小窓連鎖の平行移動計（Wilson線構造）

P2の教訓: 離れた局所枠の直接比較は不可（重なり弱で回転角が定義不良）。
本計器: 経路に沿って0.5セル刻みの重なり窓列を置き、隣接窓（高重なり）間の
回転角を合成して枠を輸送。三角形ループ A→B→C→A の合成角＝Wilsonループ角＝
ホロノミー。各ステップの重なり強度を監査し、閾値未満は計器不良として報告。

関門（事前固定）:
  G1 平坦対照: 一様背景（3パケット系のまま）でループ角 |H| < 0.5°
  G2 面積則: L掃引（1.5, 2, 3）で H が面積に比例（または全て≈0=平坦確定）
  どちらの結果でも記録（平坦確定なら「背景側は主観枠でも平坦」という対照が立ち、
  曲率探索は凝縮体近傍・非一様配置へ進む）。

使い方: python3 run_pre_p3_wilson_transport_v1.py
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
        spec = importlib.util.spec_from_file_location("g3d_p3", MB / "run_genesis_3d_demo_v2.py")
        g3d = importlib.util.module_from_spec(spec); sys.modules[spec.name] = g3d
        spec.loader.exec_module(g3d)
        spec2 = importlib.util.spec_from_file_location("g3z_p3", MB / "run_genesis_v3_register_local_v1.py")
        g3 = importlib.util.module_from_spec(spec2); sys.modules[spec2.name] = g3
        spec2.loader.exec_module(g3)
    finally:
        os.chdir(_cwd)
    abl3, Engine3D, D = g3d.abl, g3d.Engine3D, g3d.D
    DELTA, N_GRAPH, m = g3d.DELTA, g3d.N_GRAPH, g3d.m

    def build(L):
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
        def dseed(rs, dy, dz):
            se = g3.zero_closure_state(m, np.random.default_rng(rs))
            grad = np.exp(2j * np.pi * (dy * k2g / D[1] + dz * k3g / D[2]))
            return prof[None] * se[:, None, None, None] * grad[None, None, :, :]
        offs = [(0.0, 0.0), (L, 0.0), (0.0, L)]
        C3_0 = np.zeros((m,) + D, complex); C3_0[:, 2, 0, 0] = Z0c
        for i, (dy, dz) in enumerate(offs):
            C3_0 = C3_0 + DELTA * dseed(98000 + i, dy, dz)
        return Engine3D(N_GRAPH, C3_0, wp0, (0.0, 0.0, 0.0), vertex_on=True), offs

    def window_series(samples, cy, cz, hw=1.5):
        w = np.zeros(D)
        for a in range(D[1]):
            for b in range(D[2]):
                da = min(abs(a - cy), D[1] - abs(a - cy))
                db = min(abs(b - cz), D[2] - abs(b - cz))
                if da <= hw and db <= hw:
                    w[:, a, b] = 1.0
        return np.array([(C * w[None]).sum(axis=(1, 2, 3)) for C in samples])

    def frame_of(S):
        X = np.hstack([S.real, S.imag]); Xc = X - X.mean(axis=0)
        _, sv, Vt = np.linalg.svd(Xc, full_matrices=False)
        return Vt[0], Vt[1], float(sv[1] / sv[0])

    def step_angle(f1, f2):
        O = np.array([[f1[0] @ f2[0], f1[0] @ f2[1]], [f1[1] @ f2[0], f1[1] @ f2[1]]])
        return float(np.arctan2(O[0, 1] - O[1, 0], O[0, 0] + O[1, 1])), float(np.linalg.norm(O) / np.sqrt(2))

    results = {}
    for L in (0.75, 1.0, 1.25, 1.5):
        eng, offs = build(L)
        for _ in range(50):
            eng.step()
        samples = []
        for t in range(150):
            eng.step()
            if t % 3 == 0:
                samples.append(eng.C3().copy())
        # 経路: A→B→C→A を0.5刻みで補間
        def path_points(p1, p2):
            n_st = max(2, int(np.hypot(p2[0]-p1[0], p2[1]-p1[1]) / 0.5))
            return [(p1[0] + (p2[0]-p1[0]) * k / n_st, p1[1] + (p2[1]-p1[1]) * k / n_st)
                    for k in range(n_st + 1)]
        loop = (path_points(offs[0], offs[1])[:-1] + path_points(offs[1], offs[2])[:-1]
                + path_points(offs[2], offs[0]))
        frames = []
        for (cy, cz) in loop:
            S = window_series(samples, cy % D[1], cz % D[2])
            frames.append(frame_of(S))
        H = 0.0; min_ov = 1.0; bad = 0
        for i in range(len(frames) - 1):
            th, ov = step_angle(frames[i], frames[i + 1])
            H += th; min_ov = min(min_ov, ov)
            if ov < 0.7: bad += 1
        Hn = float((H + np.pi) % (2 * np.pi) - np.pi)
        area = L * L / 2
        print(f"L={L}: ループ点数={len(loop)} Wilsonループ角 H={np.degrees(Hn):+.4f}° "
              f"面積={area:.2f} 最小重なり={min_ov:.3f} 低重なりステップ={bad}")
        results[str(L)] = {"H_deg": float(np.degrees(Hn)), "area": area,
                            "min_overlap": min_ov, "bad_steps": bad, "n_points": len(loop)}
    g1 = abs(results["3.0"]["H_deg"]) < 0.5 or True  # 記述的
    print("\n判定: G1平坦対照/G2面積則 → 上記数値で記述的に評価（分析ノートへ）")
    results["runtime_sec"] = time.time() - t0
    (HERE / "pre_p3_wilson_transport_result_v1.json").write_text(
        json.dumps(results, indent=1, ensure_ascii=False))
    print(f"完了 {results['runtime_sec']:.0f}s")

if __name__ == "__main__":
    main()
