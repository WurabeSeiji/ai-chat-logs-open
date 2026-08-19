#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL7 解析 v1 — S3 偏差の成長解析（保存済み dl7_series_v1.npz のみ使用）

S3（毛ゲージ不変性・周辺化軌道の一致）は T=1500 で 3.8e-6 と判定しきい 1e-9 を超えた。
仮説: 定理5.1 は厳密だが、巻きラベルの表現差に由来する丸め（1e-16級）が
非線形力学で指数的に増幅される（stage4 の T=200 実測は 3.3e-13）。
検査: ケース間最大偏差の時系列 dev(t) を測り、T=200 時点の値と成長率を報告する。
出力: result_dl7_s3_growth_v1.json
"""
import json
import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
D = np.load(HERE / "dl7_series_v1.npz")
names = ["same_pp", "same_mm", "opp_pm", "opp_mp", "neut"]
M = np.stack([D[f"{k}_marg"] for k in names])
dev = np.max(np.abs(M - M[0]), axis=0)
i200, i500, i1000 = 199, 499, 999
res = {
    "dev_at_200": float(dev[i200]), "dev_at_500": float(dev[i500]),
    "dev_at_1000": float(dev[i1000]), "dev_at_end": float(dev[-1]),
    "doubling_steps_200_to_1000": float(800 / np.log2(max(dev[i1000], 1e-300)
                                                      / max(dev[i200], 1e-300))),
    "verdict_roundoff_amplification": bool(dev[i200] < 1e-11),
}
(HERE / "result_dl7_s3_growth_v1.json").write_text(
    json.dumps(res, indent=1, ensure_ascii=False))
print(res)
