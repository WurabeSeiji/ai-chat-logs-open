#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v11: スピン量子数の力学的判別計器——回転→海と相互作用→観測量のθ回帰周期

原理（事前記録）: 運動学的回帰（v10 Part2）は全二重項でスピノル的（自明）。
力学的判別: 種の二重項だけを U(θ)=exp(-iθσ_y/2) で回し、海と J 衝突相互作用させ、
観測量 O(θ) を測る。頂点角 φ=2R·Im(b̄a) の交差項 Im(b̄_sea a_sp) は種の二重項に
線形なので、力学は θ に一次で感応しうる:
  O(θ) の cos(θ/2) 成分（4π回帰）＝スピノル的結合の重み w_half
  O(θ) の cos(θ)   成分（2π回帰）＝ベクトル的結合
  O(θ) の cos(2θ)  成分（π回帰） ＝双線形二次結合

判定（事前固定）:
  H_spin: χ偶帯種（フェルミオン分類）の w_half が χ奇帯種（ボゾン分類）より
  有意に大きい → χパリティ＝スピン統計の力学的接続を実証。
  対照: 種+海の全体を同時に回す global 回転で O(θ) がほぼ平坦であること
  （計器の健全性）。
  観測量: O1=種帯(m=+1)パワー、O2=フェルミオンマスク総パワー、O3=Q3。

使い方: python3 run_pre_spin_dynamical_v11.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_v11", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
v1, toy, base = ex.v1, ex.toy, ex.base

S = 8.0; J_INT = 100; NTH = 33

def main() -> None:
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    shape = (n, ne)
    ms = np.arange(ne); mm = np.where(ms <= ne // 2, ms, ms - ne)
    ks = np.arange(n); kk = np.where(ks <= n // 2, ks, ks - n)
    ferm_k = (np.abs(kk) % 2 == 0) & (np.abs(kk) >= 4)
    fold3 = ((mm + 1) % 3) - 1

    def single_winding(v):
        f = np.fft.fft(v.reshape(shape), axis=0, norm="ortho")
        f[0, :] = 0.0; f[n // 2:, :] = 0.0
        return np.fft.ifft(f, axis=0, norm="ortho").reshape(v.shape)
    def project_eta(v, m_set):
        f = np.fft.fft(v.reshape(shape), axis=1, norm="ortho")
        keep = np.isin(mm, list(m_set)); f[:, ~keep] = 0.0
        return np.fft.ifft(f, axis=1, norm="ortho").reshape(v.shape)

    def make_species(chi_bins):
        a = single_winding(v1.make_bundle(sp, chi_bins, "A", scale=1.0)) * S
        b = single_winding(v1.make_bundle(sp, chi_bins, "B", scale=1.0)) * S
        a = project_eta(a, {1}); b = project_eta(b, {1})
        pw = float(np.sum(np.abs(a) ** 2) + np.sum(np.abs(b) ** 2))
        sc = np.sqrt(64.0 / pw)
        return a * sc, b * sc

    sp_even = make_species((30, 32, 34))   # χ偶=フェルミオン分類帯
    sp_odd = make_species((29, 31, 33))    # χ奇=ボゾン分類帯
    sea_a = project_eta(v1.make_bundle(sp, (37, 39, 41), "A", scale=1.0) * S, {0})
    sea_b = project_eta(v1.make_bundle(sp, (37, 39, 41), "B", scale=1.0) * S, {0})
    pws = float(np.sum(np.abs(sea_a) ** 2) + np.sum(np.abs(sea_b) ** 2))
    scs = np.sqrt(16.0 / pws); sea_a *= scs; sea_b *= scs

    thetas = np.linspace(0, 4 * np.pi, NTH)

    def observables(a, b):
        fa = np.fft.fft(np.fft.fft(a.reshape(shape), axis=0, norm="ortho"), axis=1, norm="ortho")
        fb = np.fft.fft(np.fft.fft(b.reshape(shape), axis=0, norm="ortho"), axis=1, norm="ortho")
        P = (np.abs(fa) ** 2 + np.abs(fb) ** 2) / 2
        O1 = float(P[:, mm == 1].sum())
        O2 = float(P[ferm_k, :].sum())
        Pm = P.sum(axis=0)
        O3 = float(np.sum(fold3 * Pm))
        return O1, O2, O3

    def harmonics(O):
        O = np.array(O); Oc = O - O.mean()
        c_half = 2 * np.mean(Oc * np.cos(thetas / 2))
        c_one = 2 * np.mean(Oc * np.cos(thetas))
        c_two = 2 * np.mean(Oc * np.cos(2 * thetas))
        tot = abs(c_half) + abs(c_one) + abs(c_two) + 1e-300
        return c_half, c_one, c_two, abs(c_half) / tot

    out = {"J_INT": J_INT, "NTH": NTH, "cases": {}}
    for label, (spa, spb), rotate_all in (
            ("χ偶種(F分類)", sp_even, False),
            ("χ奇種(B分類)", sp_odd, False),
            ("対照:全体回転", sp_even, True)):
        curves = {"O1": [], "O2": [], "O3": []}
        for th in thetas:
            c, s_ = np.cos(th / 2), np.sin(th / 2)
            # 二重項回転: (a,b)→(c·a−s·b, s·a+c·b)
            asp = c * spa - s_ * spb
            bsp = s_ * spa + c * spb
            if rotate_all:
                a = c * (spa + sea_a) - s_ * (spb + sea_b)
                b = s_ * (spa + sea_a) + c * (spb + sea_b)
            else:
                a = asp + sea_a
                b = bsp + sea_b
            for _ in range(J_INT):
                a, b, _ = ex.collision_step_exact(a, b, sp)
            O1, O2, O3 = observables(a, b)
            curves["O1"].append(O1); curves["O2"].append(O2); curves["O3"].append(O3)
        res = {}
        for name, O in curves.items():
            ch, c1, c2, w = harmonics(O)
            rng = float(max(O) - min(O))
            res[name] = {"c_half": float(ch), "c_one": float(c1), "c_two": float(c2),
                          "w_half": float(w), "range": rng, "mean": float(np.mean(O)),
                          "curve": [float(x) for x in O]}
            print(f"{label} {name}: 振れ={rng:.4f}({rng/np.mean(O)*100:.1f}%) "
                  f"c_θ/2={ch:+.4f} c_θ={c1:+.4f} c_2θ={c2:+.4f} → w_half={w:.3f}")
        out["cases"][label] = res
    # H_spin 判定（O1で代表）
    try:
        wF = out["cases"]["χ偶種(F分類)"]["O1"]["w_half"]
        wB = out["cases"]["χ奇種(B分類)"]["O1"]["w_half"]
        ctrl = out["cases"]["対照:全体回転"]["O1"]["range"]
        sig = out["cases"]["χ偶種(F分類)"]["O1"]["range"]
        h = bool(wF > wB + 0.1 and sig > 3 * ctrl)
        print(f"\nH_spin: w_half F={wF:.3f} vs B={wB:.3f}  信号/対照振れ比={sig/max(ctrl,1e-300):.1f} → {h}")
        out["H_spin"] = {"w_half_F": wF, "w_half_B": wB,
                          "signal_vs_control": sig / max(ctrl, 1e-300), "pass": h}
    except KeyError:
        pass
    out["runtime_sec"] = time.time() - t0
    (HERE / "pre_spin_dynamical_result_v11.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s → pre_spin_dynamical_result_v11.json")

if __name__ == "__main__":
    main()
