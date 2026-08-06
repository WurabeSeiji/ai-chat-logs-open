#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v4b正式版: 巻き数シフトによる純粋+2/真正混合/海入りの寿命比較
（当初ヒアドキュメント実行だったものの正式スクリプト化。JSON再現一致で検証済み）
事前登録: 孤立純粋種は任意mで自己複製(τ=∞)／海入りでτ(+1) vs τ(+2)比較。
注意（後知見）: 本海構成はQ_wind保存クラス外（v7-v8でNyquist折返しと同定）。
S1/S2の定量は保留付き。純粋種・混合の厳密安定の結論は影響なし。"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_fp4bs", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
v1, toy, base = ex.v1, ex.toy, ex.base

S = 8.0; J_TOT = 4000; J_WIN = 40

def main():
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    shape = (n, ne)
    ms = np.arange(ne); mm = np.where(ms <= ne // 2, ms, ms - ne)
    ks = np.arange(n); kk = np.where(ks <= n // 2, ks, ks - n)
    ferm_k = (np.abs(kk) % 2 == 0) & (np.abs(kk) >= 4)
    eta = 2 * np.pi * np.arange(ne) / ne
    def single_winding(v):
        f = np.fft.fft(v.reshape(shape), axis=0, norm="ortho"); f[0, :] = 0.0; f[n // 2:, :] = 0.0
        return np.fft.ifft(f, axis=0, norm="ortho").reshape(v.shape)
    def project_eta(v, m_set):
        f = np.fft.fft(v.reshape(shape), axis=1, norm="ortho")
        keep = np.isin(mm, list(m_set)); f[:, ~keep] = 0.0
        return np.fft.ifft(f, axis=1, norm="ortho").reshape(v.shape)
    def shift_eta(v, dm):
        return (v.reshape(shape) * np.exp(1j * dm * eta)[None, :]).reshape(v.shape)
    a0 = single_winding(v1.make_bundle(sp, (30, 32, 34), "A", scale=1.0)) * S
    b0 = single_winding(v1.make_bundle(sp, (30, 32, 34), "B", scale=1.0)) * S
    a0 = a0 + single_winding(v1.make_bundle(sp, (21,), "A", scale=1.0)) * (0.2 * S)
    pow0 = float(np.sum(np.abs(a0) ** 2) + np.sum(np.abs(b0) ** 2))
    a1 = project_eta(a0, {1}); b1 = project_eta(b0, {1})
    pw = float(np.sum(np.abs(a1) ** 2) + np.sum(np.abs(b1) ** 2))
    sc = np.sqrt(pow0 / pw); a1 *= sc; b1 *= sc
    sea_a = v1.make_bundle(sp, (29, 31, 33), "A", scale=1.0) * S
    sea_b = v1.make_bundle(sp, (29, 31, 33), "B", scale=1.0) * S
    sea_a = project_eta(sea_a, {0}); sea_b = project_eta(sea_b, {0})
    pws = float(np.sum(np.abs(sea_a) ** 2) + np.sum(np.abs(sea_b) ** 2))
    scs = np.sqrt(0.25 * pow0 / pws); sea_a *= scs; sea_b *= scs
    cases = {
        "A2_pure+2孤立": (shift_eta(a1, 1), shift_eta(b1, 1)),
        "M13_混合+1+3": ((a1 + shift_eta(a1, 2)) / np.sqrt(2), (b1 + shift_eta(b1, 2)) / np.sqrt(2)),
        "S1_+1+海25%": (a1 + sea_a, b1 + sea_b),
        "S2_+2+海25%": (shift_eta(a1, 1) + sea_a, shift_eta(b1, 1) + sea_b),
    }
    out = {"J_TOT": J_TOT, "J_WIN": J_WIN, "sea_power_frac": 0.25, "cases": {}}
    for name, (a, b) in cases.items():
        fa = np.fft.fft(np.fft.fft(a.reshape(shape), axis=0, norm="ortho"), axis=1, norm="ortho")
        fb = np.fft.fft(np.fft.fft(b.reshape(shape), axis=0, norm="ortho"), axis=1, norm="ortho")
        P0 = (np.abs(fa) ** 2 + np.abs(fb) ** 2) / 2
        occ0 = P0 > (P0.max() * 1e-8); F0 = ferm_k[:, None] & occ0
        w0 = {int(m_): float(np.sum(P0[F0 & (mm[None, :] == m_)])) for m_ in range(-6, 7)}
        m_dom = max(w0, key=w0.get)
        series = []; wins = []; Pwin = np.zeros(shape)
        for j in range(J_TOT):
            a, b, _ = ex.collision_step_exact(a, b, sp)
            fa = np.fft.fft(np.fft.fft(a.reshape(shape), axis=0, norm="ortho"), axis=1, norm="ortho")
            fb = np.fft.fft(np.fft.fft(b.reshape(shape), axis=0, norm="ortho"), axis=1, norm="ortho")
            Pwin += (np.abs(fa) ** 2 + np.abs(fb) ** 2) / 2
            if (j + 1) % J_WIN == 0:
                P = Pwin / J_WIN; occ = P > (P.max() * 1e-8); F = ferm_k[:, None] & occ
                series.append({int(m_): float(np.sum(P[F & (mm[None, :] == m_)])) for m_ in range(-6, 7)})
                wins.append(j + 1); Pwin = np.zeros(shape)
        dom = np.array([s[m_dom] for s in series]); tot = np.array([sum(s.values()) for s in series])
        frac = dom / np.maximum(tot, 1e-300)
        ratio = float(dom[-1] / max(dom[0], 1e-300))
        tau = None
        if (dom > 0).all() and dom[0] > 0:
            coef = np.polyfit(np.array(wins, float), np.log(dom), 1)
            tau = float(-1.0 / coef[0]) if coef[0] < 0 else float("inf")
        last = dict(series[-1]); last.pop(m_dom, None)
        leaks = sorted(last.items(), key=lambda kv: -kv[1])[:3]
        print(f"{name}: 支配m={m_dom:+d} 保持={ratio:.3f} 純度{frac[0]:.3f}→{frac[-1]:.3f} "
              f"τ={tau if tau and tau != float('inf') else '∞'}")
        out["cases"][name] = {"m_dom": m_dom, "retention": ratio, "purity_first": float(frac[0]),
            "purity_last": float(frac[-1]), "tau": tau,
            "leaks_top3": [[m_, p_] for m_, p_ in leaks], "windows": wins, "series": series}
    # 再現一致検定（旧JSONと比較）
    old = HERE / "pre_fixedpoint_pm1_result_v4b.json"
    if old.exists():
        prev = json.loads(old.read_text())
        ok = all(abs(prev["cases"][k]["retention"] - out["cases"][k]["retention"]) < 1e-9
                 for k in out["cases"] if k in prev.get("cases", {}))
        print(f"再現一致（旧JSONと保持率一致）= {ok}")
        out["reproduction_match"] = ok
    (HERE / "pre_fixedpoint_pm1_result_v4b.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    out["runtime_sec"] = time.time() - t0
    print(f"完了 {out['runtime_sec']:.0f}s")

if __name__ == "__main__":
    main()
