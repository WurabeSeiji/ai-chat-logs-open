#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""重力読出し P6: 面積経路による曲率読出し第一号（Bargmann面積 vs 辺長ヘロン面積）

原理（事前記録・論文0の処方の状態空間実装）:
  曲率の担い手は長さでなく面積（論文0: 辺長は厳密保存・面積が歪む）。
  状態はノルム保存で閉塞球面上に住む。3状態 A,B,C について
  - 辺長（測地・厳密）: d_ij = arccos |⟨Zi,Zj⟩| （ゲージ不変スカラー）
  - 面積の独立測定: Bargmann位相 Φ = arg(⟨A,B⟩⟨B,C⟩⟨C,A⟩)
    （基底・位相・パラメタ不変＝ゲージ不要・時計不要の面積計）
  - 平坦期待面積: ヘロン公式 A_H(d_AB, d_BC, d_CA)
  曲率読出し: 超過比 Φ/(2A_H)（規約: Bloch型では囲む面積=2×位相）と
  そのサイズ依存。定曲率なら Φ/2 ÷ A_H → 1 + K·A/3 型の一定曲率が出る。

判定（事前固定・記述的）:
  G1: Φ と A_H が三角形サイズ（Δt）とともに二乗スケールで縮む（面積の資格）
  G2: 比 Φ/(2A_H) がサイズによらず一定値へ（定曲率の読出し）→ K の初測定
  比が1なら平坦、≠1なら曲率検出。乱雑ならこの三状態選びが不適（記録）。

方法: N=5 凝縮体の時間発展から3スナップショット Z(t0), Z(t0+Δ), Z(t0+2Δ)
（球面上の小三角形）。Δ掃引でサイズスケーリング。

使い方: python3 run_pre_p6_area_curvature_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
SPACE = HERE.parent / "空間軸3軸と固有時間の創生_v1"
spec1 = importlib.util.spec_from_file_location("pre1_p6", SPACE / "run_pre_2plus1_structure_v1.py")
pre1 = importlib.util.module_from_spec(spec1); sys.modules[spec1.name] = pre1
spec1.loader.exec_module(pre1)
abl = pre1.abl

def main():
    t0 = time.time()
    n = 5
    sys_lr, v, _, _, _, p, q, Z, wp = abl.build_init(n, True)
    for t in range(2000):
        Z, wp = abl.evolve(sys_lr, Z, wp)
    # 軌道スナップショット列（十分長く保存）
    traj = [Z.copy()]
    for t in range(400):
        Z, wp = abl.evolve(sys_lr, Z, wp)
        traj.append(Z.copy())
    traj = [z / np.linalg.norm(z) for z in traj]

    import math
    def measure_triangle(i0, dlt):
        A_, B_, C_ = traj[i0], traj[i0 + dlt], traj[i0 + 2 * dlt]
        oAB = np.vdot(A_, B_); oBC = np.vdot(B_, C_); oCA = np.vdot(C_, A_)
        # 辺長（Fubini-Study測地）
        dAB = math.acos(min(1, abs(oAB))); dBC = math.acos(min(1, abs(oBC)))
        dCA = math.acos(min(1, abs(oCA)))
        # Bargmann位相（囲む面積の2倍が標準規約: Φ = −2×面積/…規約は比で吸収）
        Phi = float(np.angle(oAB * oBC * oCA))
        # ヘロン（平坦期待）
        s = 0.5 * (dAB + dBC + dCA)
        h2 = s * (s - dAB) * (s - dBC) * (s - dCA)
        A_H = math.sqrt(max(h2, 0.0))
        return dAB, dBC, dCA, Phi, A_H

    print(f"{'Δ':>4} {'辺長平均':>10} {'Φ(Bargmann)':>13} {'A_H(ヘロン)':>12} {'Φ/(2A_H)':>10}")
    rows = []
    for dlt in (2, 4, 8, 16, 32, 64, 128):
        vals = []
        for i0 in (0, 50, 100):
            if i0 + 2 * dlt >= len(traj): continue
            dAB, dBC, dCA, Phi, A_H = measure_triangle(i0, dlt)
            if A_H > 1e-14:
                vals.append((0.5 * (dAB + dBC) , Phi, A_H, Phi / (2 * A_H)))
        if not vals: continue
        m = np.mean(np.array(vals), axis=0)
        print(f"{dlt:>4} {m[0]:>10.5f} {m[1]:>13.3e} {m[2]:>12.3e} {m[3]:>10.4f}")
        rows.append({"delta": dlt, "side_mean": float(m[0]), "Phi": float(m[1]),
                      "A_H": float(m[2]), "ratio": float(m[3])})
    out = {"rows": rows, "runtime_sec": time.time() - t0}
    (HERE / "pre_p6_area_curvature_result_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s")

if __name__ == "__main__":
    main()
