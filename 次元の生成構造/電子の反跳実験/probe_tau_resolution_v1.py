#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""相対位相の τ 分解能——「原理的に無理」の中身を切り分ける

問い（木原）
------------
本系は AB 二体しかなく分解能が極めて低いので強く量子化されている。
位置位相も強く量子化されているはずで、相対位相 Δθ の τ 発展をもっと細かく
読みたいが原理的に無理ではないか。

これは少なくとも三つの別々の主張に分かれる。混ぜずに一つずつ測る。

  Q1 読出しの量子化: 位置位相 θ の **値** は離散か。
     → 格子 512 セルより細かい変位に読出しが応答するかを測る。
  Q2 τ の量子化: 時間刻みは細分できるか。
     → 力学は「並進」と「衝突」の二つでできている。並進は連続量 v で
       いくらでも小さくできる。衝突 collision_step_exact は 1 步 = 1 回で
       半回はない。**細分できるのは並進だけ**という予想を検定する。
  Q3 ズームの可否: ω₀ を下げれば同じ軌道を細かく標本化できるか。
     → できるなら「分解能の問題」。同じ軌道にならないなら「原理の問題」。

検定
----
  W1 読出しの連続性: shift を 1 セル幅の 1/100 刻みで掃引し、重心角が
     階段になるか直線になるかを見る。階段なら Q1 は真。
  W2 実測 Δθ の増分分布: CR4 の 4000 步から |dΔθ| の最小値と分布を見る。
     量子があるなら最小増分に床ができる。
  W3 ω₀ ズーム: ω₀ を 1/4・1/16 にし、T を 4倍・16倍にして走らせ、
     規格化時刻で Δθ(τ/T) が重なるかを見る。重なれば細分可能。
  W4 衝突のみ（v≡0）: 並進を止め、衝突だけで Δθ がどれだけ動くかを測る。
     これが「1 步あたり削れない最小の位相移動」＝時間量子の候補。

使い方: python3 probe_tau_resolution_v1.py
出力  : result_tau_resolution_v1.json
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
UNI = HERE.parent / "統一万能関数_v1"


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


_uni = _load("uni_tr", UNI / "unified_interaction_v1.py")
K = _load("kin_tr", UNI / "unified_kinetic_v1.py")
_cr0 = _load("cr0_tr", HERE / "run_cr0_control_no_theta_v2.py")
_cr1 = _load("cr1_tr", HERE / "run_cr1_kinetic_feedback_v1.py")

DEG_A, DEG_B = -30.0, +30.0
OMEGA0 = np.pi / 72.0
PACK_A, PACK_B = tuple(range(1, 18)), tuple(range(1, 4))


def make_ab(sp, slope, icept, tag):
    case = _uni.two_body_base.explicit_packet_case(
        mode=tag, packet_a=PACK_A, packet_b=PACK_B,
        packet_a_shift=_cr0.shift_for_deg(DEG_A, slope, icept),
        packet_b_shift=_cr0.shift_for_deg(DEG_B, slope, icept))
    a = _uni.two_body_base.make_case_state(sp, case, "A", hair_enabled=True)
    b = _uni.two_body_base.make_case_state(sp, case, "B", hair_enabled=True)
    return a, b


def dtheta(a, b, n_chi, n_eta):
    pa, _ = _cr0.circle_position(a, n_chi, n_eta)
    pb, _ = _cr0.circle_position(b, n_chi, n_eta)
    return float(np.angle(np.exp(1j * (pa - pb))))


def main() -> None:
    t0 = time.time()
    sp = _uni.two_body_base.build_source_params(
        _uni.two_body_base.Params(high_n=63, recursive_collision_count=200))
    n_chi, n_eta = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, resid = _cr0.calibrate_shift(sp, n_chi, n_eta)
    cell_deg = 360.0 / n_chi
    cell_shift = cell_deg / slope          # 1 セルぶんの shift 量
    print(f"格子 n_chi={n_chi} n_eta={n_eta}  1セル={cell_deg:.4f}°  "
          f"（shift 単位で {cell_shift:.6e}）")
    out = {}

    # ---- W1 読出しの連続性 -------------------------------------------
    print("\n【W1】読出しは量子化されているか（1セル幅を 1/100 刻みで掃引）")
    for span_cells in (2.0, 0.1, 0.01):
        xs = np.linspace(0.0, span_cells * cell_shift, 201)
        ys = []
        for s in xs:
            a, _ = _cr0.make_pair(sp, s, s)
            ys.append(_cr0.circle_position(a, n_chi, n_eta)[0])
        yd = np.degrees(np.unwrap(np.asarray(ys)))
        sl, ic = np.polyfit(xs, yd, 1)
        r = float(np.max(np.abs(yd - (sl * xs + ic))))
        step = float(np.min(np.abs(np.diff(yd))))
        rng = float(yd.max() - yd.min())
        print(f"  掃引 {span_cells:5.2f} セル ({span_cells*cell_deg:9.6f}°): "
              f"応答幅 {rng:11.8f}°  直線残差 {r:.3e}°  "
              f"最小増分 {step:.3e}°")
        out[f"W1_span{span_cells}"] = {"range_deg": rng, "lin_resid_deg": r,
                                       "min_step_deg": step, "slope": float(sl)}
    print("  → 直線残差が最小増分と同程度なら連続。桁で下回れば階段（量子化）。")

    # ---- W2 実測 Δθ の増分分布 ---------------------------------------
    print("\n【W2】CR4 実測 Δθ の増分分布（4000 步・1步刻み）")
    for tag in ("case17_3", "case3_3", "case17_17"):
        p = HERE / f"cr4_relative_{tag}_data_v1.json"
        if not p.exists():
            continue
        chi = np.degrees(np.array(json.loads(p.read_text(encoding="utf-8"))["chi"]))
        d = np.abs(np.diff(chi))
        nz = d[d > 0]
        print(f"  {tag:10} 最小|dΔθ| {nz.min():.3e}°  中央 {np.median(d):.4f}°  "
              f"最大 {d.max():8.4f}°  ゼロ増分 {int((d==0).sum())}/{len(d)}")
        out[f"W2_{tag}"] = {"min_nonzero": float(nz.min()),
                            "median": float(np.median(d)),
                            "max": float(d.max()),
                            "n_zero": int((d == 0).sum()), "n": int(len(d))}

    # ---- W3 ω₀ ズーム -------------------------------------------------
    print("\n【W3】ω₀ を下げれば同じ軌道を細かく標本化できるか")
    zoom = {}
    for fac, T in ((1, 1500), (4, 6000), (16, 24000)):
        a, b = make_ab(sp, slope, icept, f"tr_zoom{fac}")
        omega, v = OMEGA0 / fac, 0.0
        chi_h, r_h = [], []
        for _ in range(T):
            chi = dtheta(a, b, n_chi, n_eta)
            r_now = float(_cr1.toy.theta_from_ab(a, b, sp).reflection_rate)
            acc = -4.0 * np.sin(omega / 2.0) ** 2 * chi
            v += acc
            omega += (1.0 - r_now) * acc
            chi_h.append(chi); r_h.append(r_now)
            a = K.k_translate_flat(a, -v, n_chi, n_eta)
            a, b, _ = _uni.collision_step_exact(a, b, sp)
        chi_h = np.degrees(np.array(chi_h))
        zoom[fac] = chi_h
        # 支配周期
        y = chi_h - chi_h.mean()
        Y = np.abs(np.fft.rfft(y * np.hanning(len(y))))
        f = np.fft.rfftfreq(len(y)); Y[0] = 0
        per = float(1.0 / f[int(np.argmax(Y))])
        print(f"  ω₀/{fac:<2d} T={T:6d}: Δθ範囲 [{chi_h.min():+8.3f},"
              f"{chi_h.max():+8.3f}]°  支配周期 {per:9.3f} 步  "
              f"周期/fac = {per/fac:8.3f}  r平均 {np.mean(r_h):.6f}")
        out[f"W3_fac{fac}"] = {"T": T, "period": per, "period_over_fac": per / fac,
                               "range": [float(chi_h.min()), float(chi_h.max())],
                               "r_mean": float(np.mean(r_h))}
    # 規格化時刻で重なるか
    print("\n  規格化時刻 τ/T で比べた軌道の一致（重なればズーム可能）")
    base = zoom[1]
    tb = np.linspace(0, 1, len(base))
    for fac in (4, 16):
        z = zoom[fac]
        zi = np.interp(tb, np.linspace(0, 1, len(z)), z)
        rms = float(np.sqrt(np.mean((zi - base) ** 2)))
        amp = float(np.sqrt(np.mean(base ** 2)))
        print(f"    ω₀/{fac:<2d} vs ω₀: RMS差 {rms:8.4f}°  "
              f"（基準RMS {amp:.4f}°、比 {rms/amp:6.3f}）")
        out[f"W3_match_{fac}"] = {"rms_diff": rms, "base_rms": amp,
                                  "ratio": rms / amp}

    # ---- W4 衝突のみ（並進を止める）------------------------------------
    print("\n【W4】並進 v≡0——衝突だけで Δθ はどれだけ動くか")
    a, b = make_ab(sp, slope, icept, "tr_collonly")
    chi_h = []
    for _ in range(1500):
        chi_h.append(dtheta(a, b, n_chi, n_eta))
        a, b, _ = _uni.collision_step_exact(a, b, sp)
    chi_h = np.degrees(np.array(chi_h))
    d = np.abs(np.diff(chi_h))
    print(f"  Δθ 範囲 [{chi_h.min():+8.4f},{chi_h.max():+8.4f}]°  "
          f"総変位 {chi_h.max()-chi_h.min():8.4f}°")
    print(f"  1步あたり |dΔθ|: 最小 {d.min():.3e}°  中央 {np.median(d):.4e}°  "
          f"最大 {d.max():.4e}°")
    out["W4_collision_only"] = {
        "range": [float(chi_h.min()), float(chi_h.max())],
        "span": float(chi_h.max() - chi_h.min()),
        "d_min": float(d.min()), "d_med": float(np.median(d)),
        "d_max": float(d.max()),
        "chi_head": [round(x, 6) for x in chi_h[:20]]}

    # 比較: 並進ありの 1 步あたり位相移動（W3 の fac=16、最も遅い走行）
    d16 = np.abs(np.diff(zoom[16]))
    print(f"  参考: ω₀/16 走行の 1步あたり |dΔθ| 中央 {np.median(d16):.4e}°")
    print("  → 衝突のみの位相移動が ω₀ を下げても消えないなら、"
          "それが削れない時間量子。")
    out["W4_ref_omega16_dmed"] = float(np.median(d16))

    out["meta"] = {"experiment": "tau_resolution_probe_v1",
                   "date": time.strftime("%Y-%m-%d %H:%M:%S"),
                   "n_chi": n_chi, "n_eta": n_eta, "cell_deg": cell_deg,
                   "calib_slope": slope, "calib_resid": resid}
    p = HERE / "result_tau_resolution_v1.json"
    p.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"\n保存 {p.name}  ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
