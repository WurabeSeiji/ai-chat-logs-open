#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""予備実験 v3系（コム測定3世代）——計器の失敗2件と確定版を1ファイルに記録

v3  位相勾配計: 欠陥計器。復調複素座標の unwrap 位相の傾きで周波数を測ったが、
    振幅がゼロ近傍を通ると位相スリップで汚染される。記録のため保存（結果は不採用）。
v3b 符号つきFFTピーク計（窓2000step）: 「0.1%精度の等間隔コム」を出したが、
    間隔Δは全NでFFTビン幅 2π/ns と厳密一致＝ビン量子化アーティファクト。
    「分解能の限界が構造をでっち上げる」実例として記録（図 fig_s5a）。
v3c 高分解能確定版（窓10倍＋Hanning＋8倍ゼロパディング）: コムは消滅し、
    全占有平面が単一固有時計 ω≈π/72/step で剛体回転（図 fig_s3a, fig_s5b）。

使い方: python3 run_pre_comb_v3_v3b_v3c.py [v3|v3b|v3c]（無指定は v3c のみ）
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec1 = importlib.util.spec_from_file_location("pre1_comb", HERE / "run_pre_2plus1_structure_v1.py")
pre1 = importlib.util.module_from_spec(spec1)
sys.modules[spec1.name] = pre1
spec1.loader.exec_module(pre1)
abl = pre1.abl


def collect(n, t_end, win, ev):
    sys_lr, v, _, _, _, p, q, Z, wp = abl.build_init(n, True)
    samples = []
    for t in range(t_end):
        Z, wp = abl.evolve(sys_lr, Z, wp)
        if win[0] <= t < win[1] and t % ev == 0:
            samples.append(Z.copy())
    S = np.array(samples)
    ns = S.shape[0]
    Sp = np.array([z - p * (p @ z) - q * (q @ z) for z in S])
    phi = np.unwrap(np.angle((S @ p) + 1j * (S @ q)))
    om_clock = float(np.polyfit(np.arange(ns), phi, 1)[0])
    X = np.hstack([Sp.real, Sp.imag])
    Xc = X - X.mean(axis=0)
    _, sv, Vt = np.linalg.svd(Xc, full_matrices=False)
    return Xc, sv, Vt, om_clock, ns


def run_v3():
    """位相勾配計（欠陥・記録用）。"""
    out = {}
    for n in (5, 6, 8):
        Xc, sv, Vt, om, ns = collect(n, 4000, (2000, 4000), 5)
        rows = []
        for k in range(0, 8, 2):
            c = (Xc @ Vt[k]) + 1j * (Xc @ Vt[k + 1])
            ph = np.unwrap(np.angle(c))
            w_signed = float(np.polyfit(np.arange(ns), ph, 1)[0])
            rows.append({"pair": k // 2 + 1, "sv_rel": float(sv[k] / sv[0]),
                          "w_over_clock_signed_conventional": w_signed / om,
                          "abs_ratio": abs(w_signed / om)})
        ratios = [r["abs_ratio"] for r in rows]
        gaps = [ratios[i] - ratios[i + 1] for i in range(len(ratios) - 1)]
        print(f"N={n}: ω_clock={om:+.4f}  |ω_k/ω_clock|={[round(x,4) for x in ratios]}")
        out[str(n)] = {"om_clock": om, "planes": rows, "gaps_over_clock": gaps}
    json.dump(out, open(HERE / "pre_signed_comb_result_v3.json", "w"),
              ensure_ascii=False, indent=1)
    print("saved v3")


def run_v3b():
    """符号つきFFTピーク計（窓2000step——ビン量子化アーティファクトの記録用）。"""
    out = {}
    for n in (5, 6, 8, 12):
        Xc, sv, Vt, om, ns = collect(n, 4000, (2000, 4000), 5)
        rows = []
        for k in range(0, 12, 2):
            if k + 1 >= Vt.shape[0]:
                break
            c = (Xc @ Vt[k]) + 1j * (Xc @ Vt[k + 1])
            F = np.abs(np.fft.fft(c - c.mean()))
            pk = int(np.argmax(F))
            f = pk / ns if pk <= ns // 2 else (pk - ns) / ns
            w = 2 * np.pi * f
            rows.append({"pair": k // 2 + 1, "sv_rel": float(sv[k] / sv[0]),
                          "abs_ratio": abs(w) / abs(om),
                          "line_purity": float(F[pk] / F.sum())})
        ratios = [r["abs_ratio"] for r in rows]
        gaps = [round(ratios[i] - ratios[i + 1], 4) for i in range(len(ratios) - 1)]
        print(f"N={n}: |ω_k|/ω_clock={[round(x,4) for x in ratios]}  間隔={gaps}")
        out[str(n)] = {"om_clock": om, "planes": rows, "gaps": gaps}
    json.dump(out, open(HERE / "pre_signed_comb_result_v3b.json", "w"),
              ensure_ascii=False, indent=1)
    print("saved v3b")


def run_v3c():
    """高分解能確定版（窓10倍＋Hanning＋8倍パディング）。"""
    PAD = 8
    out = {}
    for n in (5, 6, 8, 12):
        Xc, sv, Vt, om, ns = collect(n, 24000, (4000, 24000), 5)
        bin_w = 2 * np.pi / (ns * PAD)
        rows = []
        for k in range(0, 8, 2):
            c = (Xc @ Vt[k]) + 1j * (Xc @ Vt[k + 1])
            w_han = np.hanning(ns)
            F = np.abs(np.fft.fft((c - c.mean()) * w_han, n=ns * PAD))
            pk = int(np.argmax(F))
            f = pk / (ns * PAD) if pk <= ns * PAD // 2 else (pk - ns * PAD) / (ns * PAD)
            rows.append({"pair": k // 2 + 1, "sv_rel": float(sv[k] / sv[0]),
                          "abs_ratio": abs(2 * np.pi * f) / abs(om)})
        ratios = [r["abs_ratio"] for r in rows]
        gaps = [round(ratios[i] - ratios[i + 1], 5) for i in range(len(ratios) - 1)]
        print(f"N={n}: ω_clock={om:+.5f}  |ω_k|/ω_clock={[round(x,5) for x in ratios]}")
        out[str(n)] = {"om_clock": om, "planes": rows, "gaps": gaps,
                        "bin_over_clock": bin_w / abs(om)}
    json.dump(out, open(HERE / "pre_comb_highres_result_v3c.json", "w"),
              ensure_ascii=False, indent=1)
    print("saved v3c")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "v3c"
    t0 = time.time()
    {"v3": run_v3, "v3b": run_v3b, "v3c": run_v3c}[which]()
    print(f"({time.time()-t0:.0f}s)")
