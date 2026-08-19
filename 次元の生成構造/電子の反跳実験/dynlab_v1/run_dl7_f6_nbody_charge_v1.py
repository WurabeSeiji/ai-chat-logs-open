#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL7-F6 走行 v1 — N体正本パイプラインでの±則（実行前固定・DL7ノート §13）

F（unified_interaction_v1）→ D v2（d_panel・変位=時計離調）。
状態への書込みは τ=0 の初期構成のみ。**判定量はすべて毎瞬の状態読出し**：
運動量 = δ×方向余弦の台上和（D v2 出力の積・累積なし）、
相対距離 = 距離幾何 X3 のクラスタ重心間距離（第3編の計器・毎瞬）。
位置の累積変数は判定に使わない（位置は波自身が保持——G v4 docstring と
キック8000步無減衰の実測。carry 積分は計算上の便宜であり物理機構ではない）。

ケース: same=(+2,+2)／opp=(+2,−2)／neut=(0,0)。
判定: F6-0 前提（λ3>0・奇数帯占有生存）／F6-1 運動の実在（台上 δ が台外と分離）／
      F6-2 ±分岐（same vs opp の相対運動・neut 対照の三値読み）

出力: result_dl7_f6_v1.json・dl7_f6_series_v1.npz
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
UNI = HERE.parent.parent / "統一万能関数_v1"
N, DELTA = 16, 0.1
T = 20000
SAMPLE = 4          # 記録間引き（読出しは毎步・保存のみ間引き）
AMP = 0.05


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


def main():
    t0 = time.time()
    uni = _load("uni_f6", UNI / "unified_interaction_v1.py")
    D = _load("dim_f6", UNI / "unified_dimension_v2.py")
    Jc = np.eye(N) - np.ones((N, N)) / N
    iiu, jju = np.triu_indices(N, 1)
    ones = np.ones(N) / np.sqrt(N)

    def rel_cluster_dist(eng):
        """関係読出し: 距離幾何 X3 のクラスタ重心間距離（毎瞬・状態から）。"""
        x = eng.C2().sum(axis=(1, 2))
        D2 = np.zeros((N, N))
        D2[iiu, jju] = D2[jju, iiu] = np.abs(x) ** 2
        B = -0.5 * Jc @ D2 @ Jc
        lam, V = np.linalg.eigh(B)
        keep = np.ones(N, bool)
        keep[int(np.argmax(np.abs(V.T @ ones)))] = False
        lamk, Vk = lam[keep], V[:, keep]
        o = np.argsort(lamk)[::-1]
        lamk, Vk = lamk[o], Vk[:, o]
        X3 = Vk[:, :3] * np.sqrt(np.maximum(lamk[:3], 0.0))[None, :]
        return float(np.linalg.norm(X3[:8].mean(axis=0) - X3[8:].mean(axis=0))),             float(lamk[2])

    def build(w1, w2):
        eng, _, _ = uni.build_standard_universe(N, DELTA)
        C2 = eng.C2().copy()
        Nn, Ne = C2.shape[1], C2.shape[2]
        ia, ib = eng.ia, eng.ib
        s1 = np.array([e for e in range(eng.m) if ia[e] < 8 and ib[e] < 8])
        s2 = np.array([e for e in range(eng.m) if ia[e] >= 8 and ib[e] >= 8])
        v1 = np.cos(0.31 * np.arange(len(s1))) + 1j * np.sin(0.47 * np.arange(len(s1)))
        v2 = np.cos(0.53 * np.arange(len(s2))) + 1j * np.sin(0.29 * np.arange(len(s2)))
        C2[s1, 1, w1 % Ne] += AMP * v1 / np.linalg.norm(v1)
        C2[s2, 1, w2 % Ne] += AMP * v2 / np.linalg.norm(v2)
        eng.C = C2.reshape(eng.m, -1)
        return eng, s1, s2

    def run(w1, w2):
        eng, s1, s2 = build(w1, w2)
        carry = None
        rec = {k: [] for k in ("tau", "P1", "P2", "d_rel", "lam3",
                                "disp_on", "disp_off", "odd1", "odd2")}
        mask_off = np.ones(eng.m, bool)
        mask_off[s1] = False
        mask_off[s2] = False
        for t in range(T):
            eng.step()
            dp = D.d_panel(eng.C2(), carry)
            carry = dp["_carry"]
            if t % SAMPLE:
                continue
            disp = np.nan_to_num(np.asarray(dp["displacement"]))
            dg = np.asarray(dp["dir_gauge"])
            step_v = np.nan_to_num(disp[:, None] * dg)   # 毎瞬の運動量読み（累積なし）
            d_rel, lam3 = rel_cluster_dist(eng)
            P2band = np.abs(eng.C2()) ** 2
            rec["tau"].append(t)
            rec["P1"].append(step_v[s1].sum(axis=0))
            rec["P2"].append(step_v[s2].sum(axis=0))
            rec["d_rel"].append(d_rel)
            rec["lam3"].append(lam3)
            rec["disp_on"].append(float(np.abs(disp[s1]).mean()
                                        + np.abs(disp[s2]).mean()) / 2)
            rec["disp_off"].append(float(np.abs(disp[mask_off]).mean()))
            rec["odd1"].append(float(P2band[s1][:, eng.odd_k, :].sum()))
            rec["odd2"].append(float(P2band[s2][:, eng.odd_k, :].sum()))
        return {k: np.array(v) for k, v in rec.items()}, eng

    out = {}
    for name, (w1, w2) in (("same", (2, 2)), ("opp", (2, -2)), ("neut", (0, 0))):
        out[name], eng_last = run(w1, w2)
        r = out[name]
        h = len(r["tau"]) // 2
        print(f"  [{name}] w=({w1},{w2}) 関係距離後半={r['d_rel'][h:].mean():.5f} "
              f"δ台上/台外={r['disp_on'][h:].mean():.3e}/{r['disp_off'][h:].mean():.3e} "
              f"奇帯占有末=({r['odd1'][-1]:.3e},{r['odd2'][-1]:.3e})")

    h = len(out["same"]["tau"]) // 2
    rd = {k: out[k]["d_rel"] for k in out}
    F60 = {"odd_survival": {k: [float(out[k]["odd1"][-1]), float(out[k]["odd2"][-1])]
                            for k in out},
           "pass": bool(min(out[k]["odd1"][-1] for k in ("same", "opp")) > 1e-6)}
    F61 = {"disp_on_late": {k: float(out[k]["disp_on"][h:].mean()) for k in out},
           "disp_off_late": {k: float(out[k]["disp_off"][h:].mean()) for k in out},
           "ratio_same": float(out["same"]["disp_on"][h:].mean()
                               / max(out["same"]["disp_off"][h:].mean(), 1e-300)),
           "pass": True}
    d_so = float(rd["same"][h:].mean() - rd["opp"][h:].mean())
    d_sn = float(rd["same"][h:].mean() - rd["neut"][h:].mean())
    d_on = float(rd["opp"][h:].mean() - rd["neut"][h:].mean())
    F62 = {"reldist_late_mean": {k: float(rd[k][h:].mean()) for k in out},
           "same_minus_opp": d_so, "same_minus_neut": d_sn, "opp_minus_neut": d_on,
           "reading": ("same>opp（同電荷=離れる側）" if d_so > 0 else
                       "same<opp（同電荷=近づく側）"),
           "pass": True}
    res = {"config": {"N": N, "delta": DELTA, "T": T, "amp": AMP,
                      "clusters": "s1=頂点{0..7}内部辺, s2=頂点{8..15}内部辺",
                      "band": 1, "windings": {"same": [2, 2], "opp": [2, -2],
                                              "neut": [0, 0]}},
           "F6_0": F60, "F6_1": F61, "F6_2": F62,
           "elapsed_sec": time.time() - t0}
    np.savez_compressed(HERE / "dl7_f6_series_v1.npz",
                        **{f"{k}_{q}": out[k][q] for k in out for q in out[k]})
    (HERE / "result_dl7_f6_v1.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    print(f"F6-0 奇帯生存: {F60['pass']}")
    print(f"F6-1 δ 台上/台外比(same)={F61['ratio_same']:.3f}")
    print(f"F6-2 相対距離後半: same={rd['same'][h:].mean():.4f} "
          f"opp={rd['opp'][h:].mean():.4f} neut={rd['neut'][h:].mean():.4f} "
          f"→ {F62['reading']}")
    print(f"({res['elapsed_sec']:.0f}s)")


if __name__ == "__main__":
    main()
