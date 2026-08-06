#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""予備実験 v5: 読出し整流説の検定——素電荷±1は「分母3の観測時計」の読みか

背景（v4/v4b）: 帯電種の不安定は海駆動の倍加・反転ウォーク（±2m, ±m）。
±1一意性は状態側の力学からは出なかった。残る候補=読出し側の整流。

鍵の代数（事前記録）: 倍加梯子 2^b·(±1) の折返し
  mod 3: 2≡−1, 4≡+1, 8≡−1, … → 全生成物が ±1 に折れる（|q|=1 普遍）
  mod 4: 2^b≡0 (b≥2) → 電荷が消えて読める（非物理的）
  mod 6: 2^b≡{2,4}≡±2 → ±2 と読める
  mod 5: 2^b≡{2,4,3,1} → 4値に散る
C論文の実測（タング幅: 分母3が幅≥0.92・幅比>461で支配）は、観測時計が
事実上分母3であることを与える。よって予言:

  H_rect: J=3 折返し類 {0, +1, −1}（m mod 3）での類重みは、ウォークの間
  ほぼ保存される（時間変動が生分布より桁で小さい）。対照 J=4,5,6 では
  保存が破れるか電荷が消える。成立すれば:
  素電荷の一意性 = 「支配的観測時計（分母3）が倍加ウォークを±1に読む」。

方法: v4b の保存済み窓系列（S1: +1+海 / S2: +2+海）と v3 原構成の系列から、
巻き数重み w(m,t) を J=3,4,5,6 で折返し、類重みの時間軌道と保存度
（変動係数 CV = std/mean）を生分布の支配重み CV と比較する。
新しい力学走行は不要（保存データのみ・決定論）。

使い方: python3 run_pre_readout_rectification_v5.py
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def fold_class(m: int, J: int) -> int:
    """m mod J を対称表示（−J/2, J/2] で返す。"""
    r = m % J
    if r > J // 2:
        r -= J
    return r


def analyze_series(series, label):
    ms = sorted(int(k) for k in series[0].keys())
    W = np.array([[s[str(m)] if str(m) in s else s[m] for m in ms] for s in series])
    # 生分布: 初期支配 m の重みの軌道
    w0 = W[0]
    m_dom = ms[int(np.argmax(w0))]
    raw = W[:, ms.index(m_dom)]
    raw_cv = float(np.std(raw) / max(np.mean(raw), 1e-300))
    out = {"label": label, "m_dom": m_dom, "raw_cv": raw_cv, "folds": {}}
    print(f"\n[{label}] 生支配 m={m_dom:+d}: CV={raw_cv:.3f} "
          f"(初{raw[0]:.3f}→終{raw[-1]:.3f})")
    for J in (3, 4, 5, 6):
        classes = {}
        for i, m in enumerate(ms):
            c = fold_class(m, J)
            classes.setdefault(c, []).append(i)
        traj = {c: W[:, idx].sum(axis=1) for c, idx in classes.items()}
        # 荷電類（c≠0）の保存度
        cvs = {}
        for c, tr in sorted(traj.items()):
            mean = float(np.mean(tr))
            cv = float(np.std(tr) / max(mean, 1e-300)) if mean > 1e-9 else None
            cvs[c] = {"mean": mean, "cv": cv,
                      "first": float(tr[0]), "last": float(tr[-1])}
        charged = [c for c in cvs if c != 0 and cvs[c]["mean"] > 1e-6]
        cv_ch = [cvs[c]["cv"] for c in charged if cvs[c]["cv"] is not None]
        cv_max = max(cv_ch) if cv_ch else None
        # |q|=1 集中度（最終窓、荷電類のみ）
        tot_ch_last = sum(abs(cvs[c]["last"]) for c in cvs if c != 0)
        q1_last = sum(cvs[c]["last"] for c in cvs if abs(c) == 1)
        conc = float(q1_last / tot_ch_last) if tot_ch_last > 1e-12 else None
        out["folds"][J] = {"classes": {str(c): v for c, v in cvs.items()},
                            "cv_max_charged": cv_max, "q1_concentration_last": conc}
        cv_txt = f"{cv_max:.3f}" if cv_max is not None else "n/a"
        conc_txt = f"{conc:.3f}" if conc is not None else "n/a"
        cls_txt = ", ".join(
            f"{c:+d}:{cvs[c]['first']:.3f}→{cvs[c]['last']:.3f}"
            for c in sorted(cvs) if cvs[c]["mean"] > 1e-6)
        print(f"  J={J}: 荷電類CV最大={cv_txt}  |q|=1集中度(終)={conc_txt}")
        print(f"        類軌道: {cls_txt}")
    return out


def main() -> None:
    t0 = time.time()
    out = {"analyses": []}
    v4b = json.loads((HERE / "pre_fixedpoint_pm1_result_v4b.json").read_text())
    v3 = json.loads((HERE / "pre_charged_stability_result_v3.json").read_text())

    for name in ("S1_+1+海25%", "S2_+2+海25%"):
        out["analyses"].append(analyze_series(v4b["cases"][name]["series"], name))
    # v3: q_series 形式（+1,+3,0,−1のみ）→ 完全分布ではないので参考扱い
    qs = v3["q_series"]
    series_v3 = [{"1": qs["+1"][i], "3": qs["+3"][i], "0": qs["0"][i],
                   "-1": qs["-1"][i]} for i in range(len(qs["+1"]))]
    out["analyses"].append(analyze_series(series_v3, "v3orig(部分帳簿・参考)"))

    # 判定: S1/S2 で J=3 の荷電類CVが生CVより小さく、他Jより|q|=1集中が高いか
    verdicts = {}
    for a in out["analyses"][:2]:
        j3 = a["folds"][3]
        better_than_raw = (j3["cv_max_charged"] is not None
                           and j3["cv_max_charged"] < a["raw_cv"])
        conc_rank = []
        for J in (3, 4, 5, 6):
            c = a["folds"][J]["q1_concentration_last"]
            conc_rank.append((J, -1.0 if c is None else c))
        best_J = max(conc_rank, key=lambda x: x[1])[0]
        verdicts[a["label"]] = {"J3_cv_better_than_raw": bool(better_than_raw),
                                 "best_conc_J": int(best_J)}
        print(f"\n[{a['label']}] J3荷電類CV<生CV: {better_than_raw}  "
              f"|q|=1集中が最大のJ: {best_J}")
    out["verdicts"] = verdicts
    out["runtime_sec"] = time.time() - t0
    (HERE / "pre_readout_rectification_result_v5.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"\n完了 {out['runtime_sec']:.1f}s → pre_readout_rectification_result_v5.json")


if __name__ == "__main__":
    main()
