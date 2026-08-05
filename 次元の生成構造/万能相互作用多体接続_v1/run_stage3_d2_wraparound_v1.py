#!/usr/bin/env python3
"""段階3-D2′: E8特殊値の感度域実験——カスケード巻き付き再帰（v2エンジン）

背景: D2一次実験は点火窓がレジスタ規模に未到達で観測量が盲目（正直記録済み）。
    本実験は強駆動（δ=0.3）・長T（3000）でカスケード梯子（{1,2}から
    a+b−c 和則で整数梯子が昇る）をレジスタ一周（k≡0 mod Nreg の折返し）まで
    駆動し、巻き付き再帰＝自己干渉のレジームで特殊値選択性を検定する。

判定（実行前固定）:
    各特殊値 s ∈ {120,124,137,144,248} について ±1 近傍と比較:
    W1 巻き付き時刻 t_wrap（上位1割ビンの占有が 1e-10·全パワー を初超過）
       の比 ∈ [0.5,2] なら null。
    W2 巻き付き後のDC流入率（t_wrap後200stepの P₀ 増加率）の比 ∈ [0.5,2]
       なら null。巻き付き梯子は k≡0 (mod Nreg)＝DC を必ず通過するため、
       DC流入は再帰干渉の一次観測量。
    E1 探索: 137（奇）vs 136/138（偶）——梯子の折返しで k と k−Nreg の
       パリティが偶Nregでは保存・奇Nregでは反転する。奇偶系統差を記録。
    どちらかで比が帯域外なら選択性の証拠（線形相・一次レジスタ力学の
    全null系列に対する初の算術感度）。

使い方: python3 run_stage3_d2_wraparound_v1.py
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
    "s3w", HERE / "run_stage3_sharedO_v2_and_hair_v1.py")
s3 = importlib.util.module_from_spec(spec3)
sys.modules[spec3.name] = s3
spec3.loader.exec_module(s3)
abl = s3.abl
gen3 = s3.gen3
V2 = s3.VertexEngineV2

SPECIALS = [120, 124, 137, 144, 248]
T_LONG = 3000
DELTA = 0.3
POST = 200


def main() -> None:
    t0 = time.time()
    n, m = 5, 10
    _, v, _, _, _, _, _, Z0c, wp0 = abl.build_init(n, False)
    r2 = gen3.make_parent(n, seed=2)
    Csec = np.fft.fft(r2.relation_waves, axis=1) / n
    seed_state = Csec[:, 1] / np.linalg.norm(Csec[:, 1])
    all_reg = sorted(set(sum([[s - 1, s, s + 1] for s in SPECIALS], [])))
    print(f"D2′ 巻き付き再帰 Nreg={all_reg} δ={DELTA} T={T_LONG}")

    results = {}
    for reg in all_reg:
        C0 = np.zeros((m, reg), complex)
        C0[:, 2] = Z0c / np.linalg.norm(Z0c)
        C0[:, 1] = DELTA * seed_state
        eng = V2(n, C0, wp0, vertex_on=True)
        top = np.arange(reg) >= int(0.9 * reg)
        t_wrap = None
        P0_series = []
        ent_final = None
        for t in range(T_LONG):
            eng.step()
            P = np.sum(np.abs(eng.C) ** 2, axis=0)
            Ptot = P.sum()
            if t_wrap is None and P[top].sum() > 1e-10 * Ptot:
                t_wrap = t + 1
            if t_wrap is not None and len(P0_series) <= POST:
                P0_series.append(float(P[0]))
            if t == T_LONG - 1:
                p = P / Ptot
                ent_final = float(-np.sum(p[p > 0] * np.log(p[p > 0])))
        dc_rate = None
        if t_wrap is not None and len(P0_series) > 20:
            arr = np.array(P0_series)
            arr = arr[arr > 0]
            if len(arr) > 20:
                tt = np.arange(len(arr), dtype=float)
                A = np.vstack([tt, np.ones_like(tt)]).T
                coef, _, _, _ = np.linalg.lstsq(A, np.log(arr), rcond=None)
                dc_rate = float(coef[0])
        results[reg] = {"t_wrap": t_wrap, "dc_rate": dc_rate,
                         "entropy_final": ent_final,
                         "entropy_max": float(np.log(reg))}
        print(f"  Nreg={reg:3d}: t_wrap={t_wrap} DC率={dc_rate if dc_rate else '—'} "
              f"S_final={ent_final:.2f}/{np.log(reg):.2f}", flush=True)

    verdicts = {}
    for s in SPECIALS:
        r_s = results[s]
        nb = [results[s - 1], results[s + 1]]
        v_ = {}
        if r_s["t_wrap"] and all(x["t_wrap"] for x in nb):
            v_["wrap_ratio"] = float(r_s["t_wrap"] / np.mean([x["t_wrap"] for x in nb]))
            v_["W1_null"] = bool(0.5 <= v_["wrap_ratio"] <= 2.0)
        else:
            v_["wrap_ratio"] = None
            v_["W1_null"] = None
        if r_s["dc_rate"] and all(x["dc_rate"] for x in nb):
            v_["dc_ratio"] = float(r_s["dc_rate"] / np.mean([x["dc_rate"] for x in nb]))
            v_["W2_null"] = bool(0.5 <= v_["dc_ratio"] <= 2.0)
        else:
            v_["dc_ratio"] = None
            v_["W2_null"] = None
        verdicts[s] = v_
        print(f"  Nreg={s}: wrap比={v_['wrap_ratio']} DC比={v_['dc_ratio']} "
              f"→ W1={v_['W1_null']} W2={v_['W2_null']}")

    out = {"SPECIALS": SPECIALS, "criteria": {"DELTA": DELTA, "T_LONG": T_LONG,
            "W1": "wrap-time ratio in [0.5,2]", "W2": "post-wrap DC rate ratio in [0.5,2]"},
           "results": {str(k): v_ for k, v_ in results.items()},
           "verdicts": {str(k): v_ for k, v_ in verdicts.items()},
           "runtime_sec": time.time() - t0}
    (HERE / "stage3_d2_wraparound_result_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
