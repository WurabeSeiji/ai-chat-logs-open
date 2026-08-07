#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""予備実験 v8: 保存破れの正体＝ηNyquist折返し（巡回保存則の実証）

理論（事前記録）: 衝突写像は各点で (a,b) を φ=2R·s (s=Im b̄a) 回す点毎回転。
連続極限で δQ_wind = 2R∫s∂_ηs = R∮∂_η(s²) = 0（全状態で保存のはず）。
離散η格子上では点毎積＝スペクトルの巡回畳み込み＝巻き数保存は mod ne。
よって整数 Q_wind の見かけの破れは、ウォークがNyquistに届いた折返しのはず。

判定（事前固定）:
  H_alias: S1 の端パワー比（|m|≥ne/2−4）が O(0.1) まで蓄積し、
  corr(ΔQ_wind, 端パワー) < −0.8。D は帯域内に留まり Q_wind 厳密保存。

使い方: python3 run_pre_aliasing_v8.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_v8s", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
v1, toy, base = ex.v1, ex.toy, ex.base

S = 8.0; J_TOT = 2000; EVERY = 50

def main() -> None:
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    shape = (n, ne)
    ms = np.arange(ne); mm = np.where(ms <= ne // 2, ms, ms - ne)
    edge = np.abs(mm) >= (ne // 2 - 4)
    def single_winding(v):
        f = np.fft.fft(v.reshape(shape), axis=0, norm="ortho")
        f[0, :] = 0.0; f[n // 2:, :] = 0.0
        return np.fft.ifft(f, axis=0, norm="ortho").reshape(v.shape)
    def project_eta(v, m_set):
        f = np.fft.fft(v.reshape(shape), axis=1, norm="ortho")
        keep = np.isin(mm, list(m_set)); f[:, ~keep] = 0.0
        return np.fft.ifft(f, axis=1, norm="ortho").reshape(v.shape)
    a0 = single_winding(v1.make_bundle(sp, (30, 32, 34), "A", scale=1.0)) * S
    b0 = single_winding(v1.make_bundle(sp, (30, 32, 34), "B", scale=1.0)) * S
    a0 = a0 + single_winding(v1.make_bundle(sp, (21,), "A", scale=1.0)) * (0.2 * S)
    pow0 = float(np.sum(np.abs(a0) ** 2) + np.sum(np.abs(b0) ** 2))
    a1 = project_eta(a0, {1}); b1 = project_eta(b0, {1})
    pw = float(np.sum(np.abs(a1) ** 2) + np.sum(np.abs(b1) ** 2))
    sc = np.sqrt(pow0 / pw); a1 *= sc; b1 *= sc
    sea_a = project_eta(v1.make_bundle(sp, (29, 31, 33), "A", scale=1.0) * S, {0})
    sea_b = project_eta(v1.make_bundle(sp, (29, 31, 33), "B", scale=1.0) * S, {0})
    pws = float(np.sum(np.abs(sea_a) ** 2) + np.sum(np.abs(sea_b) ** 2))
    scs = np.sqrt(0.25 * pow0 / pws); sea_a *= scs; sea_b *= scs
    out = {"ne": int(ne), "chi_n": int(n), "J_TOT": J_TOT, "EVERY": EVERY, "cases": {}}
    for name, (a, b) in {"D": (a0.copy(), b0.copy()),
                          "S1": (a1 + sea_a, b1 + sea_b)}.items():
        Q, E, js = [], [], []
        for j in range(J_TOT):
            a, b, _ = ex.collision_step_exact(a, b, sp)
            if j % EVERY == 0:
                fa = np.fft.fft(np.fft.fft(a.reshape(shape), axis=0, norm="ortho"), axis=1, norm="ortho")
                fb = np.fft.fft(np.fft.fft(b.reshape(shape), axis=0, norm="ortho"), axis=1, norm="ortho")
                Pm = np.sum((np.abs(fa) ** 2 + np.abs(fb) ** 2) / 2, axis=0)
                Q.append(float(np.sum(mm * Pm)))
                E.append(float(Pm[edge].sum() / Pm.sum()))
                js.append(j)
        Q = np.array(Q); E = np.array(E)
        dQ = Q - Q[0]
        r = float(np.corrcoef(dQ[1:], E[1:])[0, 1]) if E[1:].std() > 0 else None
        print(f"{name}: Q {Q[0]:+.4f}→{Q[-1]:+.4f}  端パワー比 {E[0]:.2e}→{E[-1]:.2e}  "
              f"corr(ΔQ,端)={r:+.3f}" if r is not None else f"{name}: corr n/a")
        out["cases"][name] = {"j": js, "Q_wind": Q.tolist(), "edge_frac": E.tolist(),
                               "corr_dQ_edge": r}
    c = out["cases"]
    h = bool(c["S1"]["edge_frac"][-1] > 0.1 and (c["S1"]["corr_dQ_edge"] or 0) < -0.8
             and abs(c["D"]["Q_wind"][-1] - c["D"]["Q_wind"][0]) < 1e-3)
    out["H_alias"] = h
    print(f"H_alias = {h}")
    out["runtime_sec"] = time.time() - t0
    (HERE / "pre_aliasing_result_v8.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s → pre_aliasing_result_v8.json")

if __name__ == "__main__":
    main()
