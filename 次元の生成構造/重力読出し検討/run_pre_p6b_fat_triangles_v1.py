#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""重力読出し P6b: 太い三角形による断面曲率の同定（面積計の較正と物理測定）

方向対3種（事前記録）:
  (a) 正則対 {u, iu}: 複素構造の平面。FS幾何の理論値: 小サイズ比→1（符号規約込み）
      ・正則断面曲率=最大。既知幾何による計器較正その1。
  (b) 全実対 {u1, u2}（⟨u2,iu1⟩=0）: 全実平面。理論値: 比→0（シンプレクティック
      面積が消える）＝Kähler角90°。較正その2。
  (c) 物理対 {d1, d2}（読出し平面）: 物理測定——読出し平面のKähler角と
      面積比例補正（断面曲率）。
方法: A=Z, B=N(Z+εu1), C=N(Z+εu2)。ε掃引で 比(ε) の切片（Kähler角）と
面積依存（曲率）を測る。辺∝ε・面積∝ε²（太い三角形の資格検査つき）。

使い方: python3 run_pre_p6b_fat_triangles_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time, math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
SPACE = HERE.parent / "空間軸3軸と固有時間の創生_v1"
spec1 = importlib.util.spec_from_file_location("pre1_p6b", SPACE / "run_pre_2plus1_structure_v1.py")
pre1 = importlib.util.module_from_spec(spec1); sys.modules[spec1.name] = pre1
spec1.loader.exec_module(pre1)
abl = pre1.abl

def main():
    t0 = time.time()
    n = 5
    sys_lr, v, _, _, _, p, q, Z, wp = abl.build_init(n, True)
    M = sys_lr.m
    for t in range(2000):
        Z, wp = abl.evolve(sys_lr, Z, wp)
    Z = Z / np.linalg.norm(Z)

    # 読出し平面方向（物理対）
    samples = []
    Zs, wps = Z.copy() * np.linalg.norm(Z), wp.copy()
    Zs = Z.copy(); 
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

    rng = np.random.default_rng(20260806)
    def gs(u):  # Zに直交化して正規化
        u = u - Z * np.vdot(Z, u)
        return u / np.linalg.norm(u)
    u = gs(rng.normal(size=M) + 1j * rng.normal(size=M))
    # 全実対: u1実ベース, u2は iu1 成分を除去
    u1r = gs(rng.normal(size=M) + 0j)
    u2r = rng.normal(size=M) + 0j
    u2r = u2r - Z * np.vdot(Z, u2r)
    u2r = u2r - u1r * np.vdot(u1r, u2r)
    u2r = u2r - (1j * u1r) * np.vdot(1j * u1r, u2r)
    u2r /= np.linalg.norm(u2r)

    pairs = {"(a)正則対{u,iu}": (u, 1j * u),
             "(b)全実対": (u1r, u2r),
             "(c)物理対{d1,d2}": (gs(d1c), gs(d2c - d1c * np.vdot(d1c, d2c)))}

    def tri(Za, Zb, Zc):
        oAB = np.vdot(Za, Zb); oBC = np.vdot(Zb, Zc); oCA = np.vdot(Zc, Za)
        dAB = math.acos(min(1, abs(oAB))); dBC = math.acos(min(1, abs(oBC)))
        dCA = math.acos(min(1, abs(oCA)))
        Phi = float(np.angle(oAB * oBC * oCA))
        s = 0.5 * (dAB + dBC + dCA)
        h2 = s * (s - dAB) * (s - dBC) * (s - dCA)
        return Phi, math.sqrt(max(h2, 0.0)), (dAB + dBC + dCA) / 3

    out = {"cases": {}}
    for name, (ua, ub) in pairs.items():
        print(f"== {name} ==")
        rows = []
        for eps in (0.005, 0.01, 0.02, 0.04, 0.08, 0.16):
            Zb = Z + eps * ua; Zb /= np.linalg.norm(Zb)
            Zc = Z + eps * ub; Zc /= np.linalg.norm(Zc)
            Phi, A_H, dm = tri(Z, Zb, Zc)
            ratio = Phi / (2 * A_H) if A_H > 1e-15 else float("nan")
            print(f"  ε={eps}: 辺~{dm:.4f} Φ={Phi:+.3e} A_H={A_H:.3e} 比={ratio:+.5f}")
            rows.append({"eps": eps, "side": dm, "Phi": Phi, "A_H": A_H, "ratio": ratio})
        # 切片と面積勾配（小εの3点で外挿・線形フィット vs A_H）
        As = np.array([r["A_H"] for r in rows[:4]]); Rs = np.array([r["ratio"] for r in rows[:4]])
        coef = np.polyfit(As, Rs, 1)
        print(f"  切片(Kähler角cos)={coef[1]:+.5f}  勾配(∝面積・曲率型)={coef[0]:+.3f}")
        out["cases"][name] = {"rows": rows, "intercept": float(coef[1]), "slope": float(coef[0])}
    out["runtime_sec"] = time.time() - t0
    (HERE / "pre_p6b_fat_triangles_result_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s")

if __name__ == "__main__":
    main()
