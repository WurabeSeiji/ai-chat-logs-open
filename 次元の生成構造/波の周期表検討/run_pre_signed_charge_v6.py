#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""予備実験 v6: 符号付き電荷計器——正味巻き数電荷の保存と mod3 整流の両立

背景: v5の|q|=1普遍性はパワーベース（符号なし）で、±1類が対生成的に両方
成長するため正味電荷の追跡ができなかった。本実験は符号付き汎関数を測る。

代数（事前記録）:
  和則 m*=2m_B−m_s は η運動量保存（m_B+m_B−m_s=m*）。したがって
  Q_wind = Σ_{k,m} m·P(k,m) は保存候補。mod は環準同型なので、保存するなら
  Q3 = Σ fold3(m)·P も「3単位の交換を除いて」保存する（電荷保存 mod 3）。

判定（事前固定）:
  H_cons: Q_wind の変動係数 CV < 0.02（4000衝突・全窓）→ 正味巻き数電荷は保存
  H_Q3: Q3 も同水準で保存 → 「読める電荷（mod3）」の保存則が成立
  どちらか破れたらそのまま記録（パワー重みは厳密保存量でない可能性を含む）。

使い方: python3 run_pre_signed_charge_v6.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_sc6", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
v1, toy, base = ex.v1, ex.toy, ex.base

S = 8.0; J_TOT = 4000; J_WIN = 40

def main() -> None:
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    shape = (n, ne)
    ms = np.arange(ne); mm = np.where(ms <= ne // 2, ms, ms - ne)
    fold3 = ((mm + 1) % 3) - 1  # m mod 3 → {−1,0,+1}
    eta = 2 * np.pi * np.arange(ne) / ne

    def single_winding(v):
        f = np.fft.fft(v.reshape(shape), axis=0, norm="ortho")
        f[0, :] = 0.0; f[n // 2:, :] = 0.0
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

    cases = {"S1_+1+海": (a1 + sea_a, b1 + sea_b),
             "S2_+2+海": (shift_eta(a1, 1) + sea_a, shift_eta(b1, 1) + sea_b),
             "D_v3orig": (a0.copy(), b0.copy())}
    out = {"J_TOT": J_TOT, "J_WIN": J_WIN, "cases": {}}
    for name, (a, b) in cases.items():
        Qw, Q3s, Pt = [], [], []
        Pwin = np.zeros(shape)
        for j in range(J_TOT):
            a, b, _ = ex.collision_step_exact(a, b, sp)
            fa = np.fft.fft(np.fft.fft(a.reshape(shape), axis=0, norm="ortho"), axis=1, norm="ortho")
            fb = np.fft.fft(np.fft.fft(b.reshape(shape), axis=0, norm="ortho"), axis=1, norm="ortho")
            Pwin += (np.abs(fa) ** 2 + np.abs(fb) ** 2) / 2
            if (j + 1) % J_WIN == 0:
                P = Pwin / J_WIN
                Pm = np.sum(P, axis=0)                    # 全k合算のη分布
                Qw.append(float(np.sum(mm * Pm)))
                Q3s.append(float(np.sum(fold3 * Pm)))
                Pt.append(float(np.sum(Pm)))
                Pwin = np.zeros(shape)
        Qw = np.array(Qw); Q3s = np.array(Q3s); Pt = np.array(Pt)
        cvQ = float(np.std(Qw) / max(abs(np.mean(Qw)), 1e-300))
        cvQ3 = float(np.std(Q3s) / max(abs(np.mean(Q3s)), 1e-300))
        cvP = float(np.std(Pt) / max(np.mean(Pt), 1e-300))
        print(f"{name}: Q_wind {Qw[0]:+.4f}→{Qw[-1]:+.4f} (CV={cvQ:.4f})  "
              f"Q3 {Q3s[0]:+.4f}→{Q3s[-1]:+.4f} (CV={cvQ3:.4f})  総P CV={cvP:.4f}")
        out["cases"][name] = {"Q_wind_first": float(Qw[0]), "Q_wind_last": float(Qw[-1]),
            "Q_wind_cv": cvQ, "Q3_first": float(Q3s[0]), "Q3_last": float(Q3s[-1]),
            "Q3_cv": cvQ3, "P_total_cv": cvP,
            "Q_wind_series": Qw.tolist(), "Q3_series": Q3s.tolist()}
        out["cases"][name]["H_cons"] = bool(cvQ < 0.02)
        out["cases"][name]["H_Q3"] = bool(cvQ3 < 0.02)
        print(f"    H_cons(Q_wind保存)={out['cases'][name]['H_cons']}  H_Q3(Q3保存)={out['cases'][name]['H_Q3']}")
    out["runtime_sec"] = time.time() - t0
    (HERE / "pre_signed_charge_result_v6.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s → pre_signed_charge_result_v6.json")

if __name__ == "__main__":
    main()
