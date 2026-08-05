#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GENESIS 3D デモ v2（確定版）——軸1の重心計を修正

経緯: v1（run_genesis_3d_demo_v1.py）の重心計は全軸で基本巻きモーメントを
使用していた。軸1（奇数倍音軸）では |W|² が周期 D₁/2 を持つためモーメントは
恒等的にゼロで、軸1の読み値はノイズだった（軸2・軸3は全時刻で予言と厳密一致）。
v2 は軸1のみ巻き数2モーメント（位置は mod D₁/2 = mod 4）に修正する。

判定（実行前固定）: 3軸とも位置が予言 (x₀+v·t) mod (4,8,8) と ±0.1 セル以内。

結果（2026-08-05）: t=100/200/300 の3時刻×3軸すべて偏差 0.0000、PR一定 3.39——ALL PASS。

使い方: python3 run_genesis_3d_demo_v2.py
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

# v1 のエンジン定義部を再利用（実行部の手前まで）
src = (HERE / "run_genesis_3d_demo_v1.py").read_text(encoding="utf-8")
head = src.split("# 初期状態:")[0]
ns = {"__name__": "g3d_v2"}
exec(head, ns)
abl, Engine3D, D, NF = ns["abl"], ns["Engine3D"], ns["D"], ns["NF"]
DELTA, VEL, N_GRAPH, m = ns["DELTA"], ns["VEL"], ns["N_GRAPH"], ns["m"]
odd_P3 = ns["odd_P3"]


def centroid3_v2(P):
    out = []
    for ax, Nax in enumerate(D):
        marg = P.sum(axis=tuple(a for a in range(3) if a != ax))
        if ax == 0:  # 奇数倍音軸: |W|²は周期 Nax/2 → 巻き数2モーメント、位置 mod Nax/2
            z = np.sum(marg * np.exp(2j * np.pi * 2 * np.arange(Nax) / Nax)) / marg.sum()
            out.append(float((np.angle(z) * Nax / (4 * np.pi)) % (Nax / 2)))
        else:
            z = np.sum(marg * np.exp(2j * np.pi * np.arange(Nax) / Nax)) / marg.sum()
            out.append(float((np.angle(z) * Nax / (2 * np.pi)) % Nax))
    pr = float(P.sum() ** 2 / np.sum(P ** 2))
    return out, pr


def main() -> None:
    spec_g3 = importlib.util.spec_from_file_location("g3z2", HERE / "run_genesis_v3_register_local_v1.py")
    g3 = importlib.util.module_from_spec(spec_g3)
    sys.modules[spec_g3.name] = g3
    spec_g3.loader.exec_module(g3)
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
    P0 = odd_P3(eng.C3())
    x0, pr0 = centroid3_v2(P0)
    print(f"初期: n⃗=({x0[0]:.3f},{x0[1]:.3f},{x0[2]:.3f}) PR={pr0:.2f}")
    rows = []
    ok = True
    mods = (4, 8, 8)
    for t in range(300):
        eng.step()
        if (t + 1) % 100 == 0:
            P = odd_P3(eng.C3())
            x, pr = centroid3_v2(P)
            pred = [(x0[a] + VEL[a] * (t + 1)) % mods[a] for a in range(3)]
            devs = [min(abs(x[a] - pred[a]), mods[a] - abs(x[a] - pred[a])) for a in range(3)]
            ok &= all(d <= 0.1 for d in devs)
            rows.append({"t": t + 1, "x": [round(v_, 3) for v_ in x],
                         "pred": [round(v_, 3) for v_ in pred],
                         "dev": [round(d, 4) for d in devs], "PR": round(pr, 2)})
            print(f"  t={t+1}: 実測=({x[0]:.3f},{x[1]:.3f},{x[2]:.3f}) "
                  f"予言=({pred[0]:.3f},{pred[1]:.3f},{pred[2]:.3f}) "
                  f"偏差={[f'{d:.4f}' for d in devs]} PR={pr:.2f}")
    print(f"判定: 3軸並進±0.1 = {ok}")
    json.dump({"note": "軸1の重心計を巻き数2モーメントに修正(位置mod4)。v1の軸1読み値は恒等零モーメントのノイズ",
               "x0": x0, "PR0": pr0, "VEL": list(VEL), "rows": rows, "all_pass": bool(ok)},
              open(HERE / "genesis_3d_demo_result_v2.json", "w"), indent=1)
    print("saved")


if __name__ == "__main__":
    main()
