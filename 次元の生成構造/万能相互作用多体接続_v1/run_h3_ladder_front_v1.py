#!/usr/bin/env python3
"""残穴③: E8梯子前線計器——上昇梯子の前線追跡と衝突干渉の特殊値検定

D2/D2′の計器不備（点火窓が低スライスに閉じる／負周波数梯子の誤検出）を
修正した第三世代計器。上昇梯子（kp=2, ks=1 から k=3,4,5,…）の前線
    F(t) = max{k ≤ Nreg//2 : P_k > 1e-12·P_tot}
を追跡する（Nreg//2 までに制限して負枝の別名を排除）。
上昇梯子と下降梯子（負枝）は k≈Nreg/2 で衝突する。偶Nregでは衝突点が
禁止Nyquistビンに一致し、奇Nregではビン間に落ちる——算術が幾何に
現れる最初の場所の候補。

判定（実行前固定）:
    L1 前線衝突時刻 t_coll（F が Nreg//2−1 に到達）: 特殊値 {124,137,144,248}
       の t_coll が ±1近傍平均の [0.5,2] 倍なら null。
    L2 衝突後干渉: 衝突後200stepの Nyquist近傍帯（Nreg//2±2）のパワー比
       （特殊値/近傍平均）が [0.5,2] なら null。
    E 探索: 偶奇系統（137 vs 136/138）の衝突挙動差を記録。

設定: v2エンジン、グラフN=5、δ=0.3、T=1500、
    Nreg ∈ {123,124,125,136,137,138,143,144,145,247,248,249}。

使い方: python3 run_h3_ladder_front_v1.py
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec3 = importlib.util.spec_from_file_location(
    "s3l", HERE / "run_stage3_sharedO_v2_and_hair_v1.py")
s3 = importlib.util.module_from_spec(spec3)
sys.modules[spec3.name] = s3
spec3.loader.exec_module(s3)
abl = s3.abl
gen3 = s3.gen3
V2 = s3.VertexEngineV2

SPECIALS = [124, 137, 144, 248]
REGS = sorted(set(sum([[s - 1, s, s + 1] for s in SPECIALS], [])))
T_LONG = 1500
DELTA = 0.3
POST = 200


def main() -> None:
    t0 = time.time()
    n, m = 5, 10
    _, v, _, _, _, _, _, Z0c, wp0 = abl.build_init(n, False)
    r2 = gen3.make_parent(n, seed=2)
    Csec = np.fft.fft(r2.relation_waves, axis=1) / n
    seed_state = Csec[:, 1] / np.linalg.norm(Csec[:, 1])
    print(f"H3 梯子前線 Nreg={REGS} δ={DELTA} T={T_LONG}")

    results = {}
    for reg in REGS:
        C0 = np.zeros((m, reg), complex)
        C0[:, 2] = Z0c / np.linalg.norm(Z0c)
        C0[:, 1] = DELTA * seed_state
        eng = V2(n, C0, wp0, vertex_on=True)
        half = reg // 2
        t_coll = None
        fronts = []
        ny_power = []
        ny_band = [k for k in range(max(1, half - 2), min(reg, half + 3))]
        for t in range(T_LONG):
            eng.step()
            P = np.sum(np.abs(eng.C) ** 2, axis=0)
            Pt = P.sum()
            occ = np.nonzero(P[:half + 1] > 1e-12 * Pt)[0]
            F = int(occ.max()) if len(occ) else 0
            fronts.append(F)
            if t_coll is None and F >= half - 1:
                t_coll = t + 1
            if t_coll is not None and len(ny_power) <= POST:
                ny_power.append(float(P[ny_band].sum() / Pt))
        results[reg] = {"t_coll": t_coll,
                         "front_at_100": fronts[99], "front_at_500": fronts[499],
                         "ny_band_mean": float(np.mean(ny_power)) if ny_power else None}
        print(f"  Nreg={reg:3d}: t_coll={t_coll} 前線(100/500)={fronts[99]}/{fronts[499]} "
              f"Nyquist帯={results[reg]['ny_band_mean']}", flush=True)

    verdicts = {}
    for s in SPECIALS:
        rs = results[s]
        nb = [results[s - 1], results[s + 1]]
        v_ = {}
        if rs["t_coll"] and all(x["t_coll"] for x in nb):
            v_["coll_ratio"] = float(rs["t_coll"] / np.mean([x["t_coll"] for x in nb]))
            v_["L1_null"] = bool(0.5 <= v_["coll_ratio"] <= 2.0)
        else:
            v_["coll_ratio"], v_["L1_null"] = None, None
        if rs["ny_band_mean"] and all(x["ny_band_mean"] for x in nb):
            v_["ny_ratio"] = float(rs["ny_band_mean"] / np.mean([x["ny_band_mean"] for x in nb]))
            v_["L2_null"] = bool(0.5 <= v_["ny_ratio"] <= 2.0)
        else:
            v_["ny_ratio"], v_["L2_null"] = None, None
        verdicts[s] = v_
        print(f"  Nreg={s}: 衝突比={v_['coll_ratio']} Nyquist比={v_['ny_ratio']} "
              f"→ L1={v_['L1_null']} L2={v_['L2_null']}")

    out = {"REGS": REGS, "SPECIALS": SPECIALS,
           "criteria": {"DELTA": DELTA, "T_LONG": T_LONG,
                         "L1": "collision-time ratio in [0.5,2]",
                         "L2": "post-collision Nyquist-band power ratio in [0.5,2]"},
           "results": {str(k): v_ for k, v_ in results.items()},
           "verdicts": {str(k): v_ for k, v_ in verdicts.items()},
           "runtime_sec": time.time() - t0}
    (HERE / "h3_ladder_front_result_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
