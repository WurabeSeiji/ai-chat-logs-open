#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v15: 閉じ込め検定——mod3可読性（v5整流計の拡張・第0.4版予言1）

操作的定義（事前記録）: J=3観測者の空間はη円のZ₃商。種が商上で一価 ⟺ m≡0(mod3)。
  可読パワー分率 f_read = Σ_{m≡0} P(m) / ΣP(m)
  可読電荷 Q_read = Σ_{m≡0} (m/3)·P(m)  （素電荷単位＝3生巻き）
判定（事前固定）:
  H_conf: クォーク型（m=+2種+海）は f_read(0)<0.1 から出発し、ウォークで
  複合チャネル（m≡0）が育つ分だけ f_read が増加（=ハドロン化のみが可読）。
  電子型（m=+3種+海）は f_read(0)>0.6 で維持。
  孤立対照: 純粋m=+2孤立は f_read=0 のまま（自己複製・複合なし）。

使い方: python3 run_pre_confinement_v15.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_v15", UIM / "run_ignition_fate_exact_v3.py")
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
    eta = 2 * np.pi * np.arange(ne) / ne
    read_mask = (mm % 3 == 0)

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
    sea_a = project_eta(v1.make_bundle(sp, (29, 31, 33), "A", scale=1.0) * S, {0})
    sea_b = project_eta(v1.make_bundle(sp, (29, 31, 33), "B", scale=1.0) * S, {0})
    pws = float(np.sum(np.abs(sea_a) ** 2) + np.sum(np.abs(sea_b) ** 2))
    scs = np.sqrt(0.25 * pow0 / pws); sea_a *= scs; sea_b *= scs

    cases = {
        "quark型m=+2+海": (shift_eta(a1, 1) + sea_a, shift_eta(b1, 1) + sea_b),
        "electron型m=+3+海": (shift_eta(a1, 2) + sea_a, shift_eta(b1, 2) + sea_b),
        "quark型m=+2孤立": (shift_eta(a1, 1), shift_eta(b1, 1)),
    }
    out = {"J_TOT": J_TOT, "J_WIN": J_WIN, "cases": {}}
    for name, (a, b) in cases.items():
        fr, Qr, wins = [], [], []
        Pwin = np.zeros(shape)
        for j in range(J_TOT):
            a, b, _ = ex.collision_step_exact(a, b, sp)
            fa = np.fft.fft(np.fft.fft(a.reshape(shape), axis=0, norm="ortho"), axis=1, norm="ortho")
            fb = np.fft.fft(np.fft.fft(b.reshape(shape), axis=0, norm="ortho"), axis=1, norm="ortho")
            Pwin += (np.abs(fa) ** 2 + np.abs(fb) ** 2) / 2
            if (j + 1) % J_WIN == 0:
                P = Pwin / J_WIN
                Pm = P.sum(axis=0)
                # 海（m=0）を除いた可読分率（種側の可読性を見る）
                tot_ns = float(Pm[mm != 0].sum())
                read_ns = float(Pm[read_mask & (mm != 0)].sum())
                fr.append(read_ns / max(tot_ns, 1e-300))
                Qr.append(float(np.sum((mm[read_mask] / 3) * Pm[read_mask])))
                wins.append(j + 1)
                Pwin = np.zeros(shape)
        fr = np.array(fr); Qr = np.array(Qr)
        print(f"{name}: f_read {fr[0]:.4f}→{fr[-1]:.4f}  Q_read {Qr[0]:+.4f}→{Qr[-1]:+.4f}")
        out["cases"][name] = {"windows": wins, "f_read": fr.tolist(),
                               "Q_read": Qr.tolist(),
                               "f_first": float(fr[0]), "f_last": float(fr[-1])}
    c = out["cases"]
    h = bool(c["quark型m=+2+海"]["f_first"] < 0.1
             and c["quark型m=+2+海"]["f_last"] > c["quark型m=+2+海"]["f_first"] + 0.1
             and c["electron型m=+3+海"]["f_first"] > 0.6
             and c["quark型m=+2孤立"]["f_last"] < 0.01)
    print(f"H_conf（閉じ込め=mod3可読性）= {h}")
    out["H_conf"] = h
    out["runtime_sec"] = time.time() - t0
    (HERE / "pre_confinement_result_v15.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s → pre_confinement_result_v15.json")

if __name__ == "__main__":
    main()
