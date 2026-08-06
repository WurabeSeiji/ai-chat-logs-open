#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""重力読出し P1: 多点マッピングの完全性検分（予備実験）

目的（設計方針_v1 第1段）: 同時刻面上で複数点の xyz を同時に読めるかの検分。
読出し不調はマッピングの不完全として内容を特定する（理論の反証とはしない）。

Part A（主観チャートの多点重畳検定・N体凝縮体）:
  同一quadrature平面に2つの印（位相の異なる摂動）を同時に注入し、
  平面読出し r が線形重畳 r12 = r1 + r2 になるか機械精度で検定。
  重畳が厳密なら、現行の主観チャートは「1平面=1点」であり、同一平面上の
  多点は原理的に分離不能——マッピング不完全の内容が確定する。

Part B（レジスタ格子の三点三角形・背景側の対照）:
  3Dデモ（レジスタ格子 D=(8,8,8)）で、位相勾配により (k2,k3) 面内に変位した
  静止パケット3個を共存させ、t=100の倍数（整数セル同時刻面）で各領域重心を読む。
  判定（記述的）: (a) 3点が分離して読めるか（領域重心の安定）
  (b) 3点の距離→余弦定理→内角和が 180° に一致するか（背景の平坦対照）。

使い方: python3 run_pre_p1_multipoint_mapping_v1.py
"""
from __future__ import annotations
import importlib.util, json, os, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
SPACE = HERE.parent / "空間軸3軸と固有時間の創生_v1"
MB = HERE.parent / "万能相互作用多体接続_v1"

def main():
    t0 = time.time()
    out = {}

    # ---------- Part A ----------
    spec1 = importlib.util.spec_from_file_location("pre1_p1", SPACE / "run_pre_2plus1_structure_v1.py")
    pre1 = importlib.util.module_from_spec(spec1); sys.modules[spec1.name] = pre1
    spec1.loader.exec_module(pre1)
    abl = pre1.abl
    n = 5
    sys_lr, v, _, _, _, p, q, Z, wp = abl.build_init(n, True)
    M = sys_lr.m
    for t in range(2000):
        Z, wp = abl.evolve(sys_lr, Z, wp)
    # 平面方向（占有SVD上位2）
    samples = []
    Zs, wps = Z.copy(), wp.copy()
    for t in range(400):
        Zs, wps = abl.evolve(sys_lr, Zs, wps)
        if t % 5 == 0:
            samples.append(Zs.copy())
    S = np.array(samples)
    Sp = S - np.outer(S @ p, p) - np.outer(S @ q, q)
    X = np.hstack([Sp.real, Sp.imag]); Xc = X - X.mean(axis=0)
    _, sv, Vt = np.linalg.svd(Xc, full_matrices=False)
    d1c = Vt[0][:M] + 1j * Vt[0][M:]; d1c /= np.linalg.norm(d1c)
    d2c = Vt[1][:M] + 1j * Vt[1][M:]; d2c /= np.linalg.norm(d2c)
    eps = 1e-3
    mark = lambda phi: eps * (np.cos(phi) * d1c.conj() + np.sin(phi) * d2c.conj())
    def readout(dZ):
        Zp = Z + dZ
        r = np.dot(np.conj(d1c), Zp) + 1j * np.dot(np.conj(d2c), Zp)
        r0 = np.dot(np.conj(d1c), Z) + 1j * np.dot(np.conj(d2c), Z)
        return r - r0
    r1 = readout(mark(0.0))
    r2 = readout(mark(np.pi / 2))
    r12 = readout(mark(0.0) + mark(np.pi / 2))
    lin_res = abs(r12 - (r1 + r2)) / max(abs(r1) + abs(r2), 1e-300)
    print(f"[A] 同一平面2印の重畳: |r12-(r1+r2)|/(|r1|+|r2|) = {lin_res:.2e}")
    print(f"    → 線形重畳{'厳密' if lin_res < 1e-10 else 'ずれあり'}: "
          f"同一平面上の多点は現行チャートでは分離不能（1平面=1点）")
    out["A"] = {"linearity_residual": float(lin_res),
                 "finding": "1平面=1点（同一平面多点は線形重畳で分離不能）" if lin_res < 1e-10
                             else "非線形応答あり（要精査）"}

    # ---------- Part B ----------
    _cwd = os.getcwd(); os.chdir(MB)
    try:
        spec = importlib.util.spec_from_file_location("g3d_p1", MB / "run_genesis_3d_demo_v2.py")
        g3d = importlib.util.module_from_spec(spec); sys.modules[spec.name] = g3d
        spec.loader.exec_module(g3d)
        spec2 = importlib.util.spec_from_file_location("g3z_p1", MB / "run_genesis_v3_register_local_v1.py")
        g3 = importlib.util.module_from_spec(spec2); sys.modules[spec2.name] = g3
        spec2.loader.exec_module(g3)
    finally:
        os.chdir(_cwd)
    abl3, Engine3D, D = g3d.abl, g3d.Engine3D, g3d.D
    DELTA, N_GRAPH, m = g3d.DELTA, g3d.N_GRAPH, g3d.m
    odd_P3, centroid3_v2 = g3d.odd_P3, g3d.centroid3_v2

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

    offsets = [(0.0, 0.0), (3.0, 0.0), (0.0, 3.0)]
    C3_0 = np.zeros((m,) + D, complex)
    C3_0[:, 2, 0, 0] = Z0c
    for i, (dy, dz) in enumerate(offsets):
        C3_0 = C3_0 + DELTA * displaced_seed(98000 + i, dy, dz)
    eng = Engine3D(N_GRAPH, C3_0, wp0, (0.0, 0.0, 0.0), vertex_on=True)

    def region_centroids(P):
        cents = []
        for (dy, dz) in offsets:
            w = np.zeros(D)
            for a in range(D[1]):
                for b in range(D[2]):
                    da = min(abs(a - dy), D[1] - abs(a - dy))
                    db = min(abs(b - dz), D[2] - abs(b - dz))
                    if da <= 1.5 and db <= 1.5:
                        w[:, a, b] = 1.0
            Pw = P * w
            tot = Pw.sum()
            if tot < 1e-12:
                cents.append(None); continue
            cy = float((Pw.sum(axis=(0, 2)) @ np.arange(D[1])) / tot)
            cz = float((Pw.sum(axis=(0, 1)) @ np.arange(D[2])) / tot)
            cents.append((cy, cz))
        return cents

    def tri_angles(pts):
        import math
        a = math.dist(pts[1], pts[2]); b = math.dist(pts[0], pts[2]); c = math.dist(pts[0], pts[1])
        A = math.acos(max(-1, min(1, (b*b + c*c - a*a) / (2*b*c))))
        B = math.acos(max(-1, min(1, (a*a + c*c - b*b) / (2*a*c))))
        C = math.acos(max(-1, min(1, (a*a + b*b - c*c) / (2*a*b))))
        return math.degrees(A + B + C), (a, b, c)

    rows = []
    P = odd_P3(eng.C3())
    cents = region_centroids(P)
    print(f"[B] t=0 領域重心: {[(round(c[0],3), round(c[1],3)) if c else None for c in cents]}")
    if all(c is not None for c in cents):
        s0, dists0 = tri_angles(cents)
        print(f"    内角和 = {s0:.4f}°  辺長 = {tuple(round(x,3) for x in dists0)}")
        rows.append({"t": 0, "cents": cents, "angle_sum_deg": s0, "sides": dists0})
    for t in range(1, 301):
        eng.step()
        if t % 100 == 0:
            P = odd_P3(eng.C3())
            cents = region_centroids(P)
            ok = all(c is not None for c in cents)
            if ok:
                s_, dists_ = tri_angles(cents)
                print(f"    t={t}: 重心={[(round(c[0],3), round(c[1],3)) for c in cents]} "
                      f"内角和={s_:.4f}° 辺={tuple(round(x,3) for x in dists_)}")
                rows.append({"t": t, "cents": cents, "angle_sum_deg": s_, "sides": dists_})
            else:
                print(f"    t={t}: 領域重心の一部が消失（分離不能）")
                rows.append({"t": t, "cents": cents, "angle_sum_deg": None})
    out["B"] = {"offsets": offsets, "rows": rows}
    out["runtime_sec"] = time.time() - t0
    (HERE / "pre_p1_multipoint_mapping_result_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False, default=str))
    print(f"完了 {out['runtime_sec']:.0f}s → pre_p1_multipoint_mapping_result_v1.json")

if __name__ == "__main__":
    main()
