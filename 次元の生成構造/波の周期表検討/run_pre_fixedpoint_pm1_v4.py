#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""予備実験 v4: ±1=和則不動点仮説の寿命検定

背景（v3）: 帯電種の崩壊チャネルは和則ウォーク m*=2m_B−m_s による高巻き数への
漏れ（+1→+3 をライブ観測）。仮説: ±1 はウォークの不動点（m_B=m_s=±1→m*=±1）
だから素電荷が ±1 に整流される。

ただし注意（事前記録）: 純粋 m 状態は任意の m で自己複製する
（m_B=m_s=m → m*=m）。よって不動点論だけでは ±1 の一意性は出ない。
一意性が本物なら追加機構（海 m=0 との頂点: 2m·0系→±2m, −m の倍加ウォーク、
φ(n)≤2^b 選択則との接続）が要る。本実験はそこを測り分ける。

方法: v3 と同一の帯電構成に η巻き数射影を施した4ケースを同条件で J=4000 走行:
  A: 純粋 m=+1（不動点・素電荷候補）
  B: 純粋 m=+2（不動点だが非±1——A と同寿命なら±1一意性は不動点論の外）
  C: 混合 {+1,+3}（非不動点——速い漏れを予言）
  D: v3原構成（基準）
測定: 窓40×100 の巻き数分布 P(m)（フェルミオン的マスク上）、
  初期支配巻き数の保持率、指数寿命 τ、漏れ先トップ3。

判定（事前固定）:
  H_fp: τ_A, τ_B ≫ τ_C（純粋=自己複製が混合より長寿命）→ 不動点機構は実在
  H_pm1: τ_A ≫ τ_B（±1 が +2 より長寿命）→ ±1一意性は力学に内在
        τ_A ≈ τ_B → 一意性は不動点の外（倍加ウォーク・選択則へ）

使い方: python3 run_pre_fixedpoint_pm1_v4.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_fp4", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
v1, toy, base = ex.v1, ex.toy, ex.base

S = 8.0
J_TOT = 4000
J_WIN = 40


def main() -> None:
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    shape = (n, ne)
    ms = np.arange(ne); mm = np.where(ms <= ne // 2, ms, ms - ne)
    ks = np.arange(n); kk = np.where(ks <= n // 2, ks, ks - n)
    ferm_k = (np.abs(kk) % 2 == 0) & (np.abs(kk) >= 4)

    def single_winding(v):
        f = np.fft.fft(v.reshape(shape), axis=0, norm="ortho")
        f[0, :] = 0.0
        f[n // 2:, :] = 0.0
        return np.fft.ifft(f, axis=0, norm="ortho").reshape(v.shape)

    def base_state():
        a = single_winding(v1.make_bundle(sp, (30, 32, 34), "A", scale=1.0)) * S
        b = single_winding(v1.make_bundle(sp, (30, 32, 34), "B", scale=1.0)) * S
        a = a + single_winding(v1.make_bundle(sp, (21,), "A", scale=1.0)) * (0.2 * S)
        return a, b

    def project_eta(v, m_set):
        f = np.fft.fft(v.reshape(shape), axis=1, norm="ortho")
        keep = np.isin(mm, list(m_set))
        f[:, ~keep] = 0.0
        return np.fft.ifft(f, axis=1, norm="ortho").reshape(v.shape)

    a0, b0 = base_state()
    pow0 = float(np.sum(np.abs(a0) ** 2) + np.sum(np.abs(b0) ** 2))

    def prep(m_set):
        a = project_eta(a0, m_set)
        b = project_eta(b0, m_set)
        pw = float(np.sum(np.abs(a) ** 2) + np.sum(np.abs(b) ** 2))
        if pw <= 1e-12:
            return None, None
        sc = np.sqrt(pow0 / pw)
        return a * sc, b * sc

    cases = {}
    cases["A_pure+1"] = prep({1})
    cases["B_pure+2"] = prep({2})
    cases["C_mix+1+3"] = prep({1, 3})
    cases["D_v3orig"] = (a0.copy(), b0.copy())

    out = {"S": S, "J_TOT": J_TOT, "J_WIN": J_WIN, "cases": {}}
    for name, (a, b) in cases.items():
        if a is None:
            print(f"{name}: 射影後パワー零——スキップ")
            continue
        # 初期支配巻き数
        fa = np.fft.fft(np.fft.fft(a.reshape(shape), axis=0, norm="ortho"),
                        axis=1, norm="ortho")
        fb = np.fft.fft(np.fft.fft(b.reshape(shape), axis=0, norm="ortho"),
                        axis=1, norm="ortho")
        P0 = (np.abs(fa) ** 2 + np.abs(fb) ** 2) / 2
        occ0 = P0 > (P0.max() * 1e-8)
        F0 = ferm_k[:, None] & occ0
        w0 = {int(m_): float(np.sum(P0[F0 & (mm[None, :] == m_)]))
              for m_ in range(-6, 7)}
        m_dom = max(w0, key=w0.get)

        series = []
        Pwin = np.zeros(shape)
        wins = []
        for j in range(J_TOT):
            a, b, _ = ex.collision_step_exact(a, b, sp)
            fa = np.fft.fft(np.fft.fft(a.reshape(shape), axis=0, norm="ortho"),
                            axis=1, norm="ortho")
            fb = np.fft.fft(np.fft.fft(b.reshape(shape), axis=0, norm="ortho"),
                            axis=1, norm="ortho")
            Pwin += (np.abs(fa) ** 2 + np.abs(fb) ** 2) / 2
            if (j + 1) % J_WIN == 0:
                P = Pwin / J_WIN
                occ = P > (P.max() * 1e-8)
                F = ferm_k[:, None] & occ
                w = {int(m_): float(np.sum(P[F & (mm[None, :] == m_)]))
                     for m_ in range(-6, 7)}
                series.append(w)
                wins.append(j + 1)
                Pwin = np.zeros(shape)

        dom = np.array([s[m_dom] for s in series])
        tot = np.array([sum(s.values()) for s in series])
        frac = dom / np.maximum(tot, 1e-300)
        ratio = float(dom[-1] / max(dom[0], 1e-300))
        tau = None
        if (dom > 0).all() and dom[0] > 0:
            coef = np.polyfit(np.array(wins, float), np.log(dom), 1)
            tau = float(-1.0 / coef[0]) if coef[0] < 0 else float("inf")
        # 漏れ先: 最終窓で初期支配以外のトップ3
        last = dict(series[-1])
        last.pop(m_dom, None)
        leaks = sorted(last.items(), key=lambda kv: -kv[1])[:3]
        print(f"{name}: 支配m={m_dom:+d} 保持率(重み比)={ratio:.3f} "
              f"純度 {frac[0]:.3f}→{frac[-1]:.3f} τ={tau if tau else 'n/a'}")
        print(f"    漏れ先トップ3: {[(m_, round(p_,4)) for m_, p_ in leaks]}")
        out["cases"][name] = {"m_dom": m_dom, "w0": w0,
                               "windows": wins, "series": series,
                               "retention": ratio, "purity_first": float(frac[0]),
                               "purity_last": float(frac[-1]), "tau": tau,
                               "leaks_top3": [[m_, p_] for m_, p_ in leaks]}

    cA = out["cases"].get("A_pure+1"); cB = out["cases"].get("B_pure+2")
    cC = out["cases"].get("C_mix+1+3")
    if cA and cB and cC:
        tA = cA["tau"] or float("inf"); tB = cB["tau"] or float("inf")
        tC = cC["tau"] or float("inf")
        h_fp = bool(min(tA, tB) > 3 * tC)
        h_pm1 = bool(tA > 3 * tB)
        print(f"\nH_fp（純粋≫混合）= {h_fp}   H_pm1（+1≫+2）= {h_pm1}")
        out["H_fp"] = h_fp; out["H_pm1"] = h_pm1
    out["runtime_sec"] = time.time() - t0
    (HERE / "pre_fixedpoint_pm1_result_v4.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s → pre_fixedpoint_pm1_result_v4.json")


if __name__ == "__main__":
    main()
