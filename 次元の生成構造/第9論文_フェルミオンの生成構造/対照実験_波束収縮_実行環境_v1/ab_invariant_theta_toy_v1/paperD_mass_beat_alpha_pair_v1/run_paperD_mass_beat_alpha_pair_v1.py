#!/usr/bin/env python3
"""論文D予備実験 D6：二つのα根状態間の質量ビート v1

主張候補:
    二つのα根状態（R_{124,23} と R_{620,117}）の間の「質量分裂」型の読み
    （二腕重なりのビート）は、電荷のランニング1ステップと同一の角
    Δθ = π/310 で駆動される。二文法は同じ一つの数の相補的な面を読む。

先行公表値（R_mass_beat_two_arm_pre_v1、アンカー）:
    ΔR(124根, 620根) = 0.009352736、二腕重なりの深いヌル min|O| ~ 3.6e-14。
    （同実験の位相アンラップ法は π 跳びに支配されるため保留扱い——本実験は
      封筒とストロボFFTで置き換える）

予言（測定前固定）:
    P1（厳密恒等式）: ΔR = sin(58π/155)·sin(π/310)
        ——電荷ステップの因数分解: (位置因子)×(ビート角の正弦)
    P2（周期）: 回転は可換なので O(j) = <ψ124(j)|ψ620(j)> は厳密に周期
        620（= 2π/Δθ, Δθ = |195-193|π/620 = π/310）。|O(j+620)-O(j)| ≤ 1e-10
    P3（ビート周波数）: Re O(j) の窓 1240 FFT のピークは ビン 2
        （= Δ住所 = 電荷ステップと同一）
    P4（深いヌル）: min_j |O(j)| ≤ 1e-10（公表値の再現）
"""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import sys
from math import pi
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SEARCH_PATH = HERE.parent / "inverse_initial_conditions_v1" / "search_initial_conditions_and_plot_v1.py"

spec = importlib.util.spec_from_file_location("search_for_paperD6_v1", SEARCH_PATH)
search = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = search
spec.loader.exec_module(search)
toy = search.toy
base = toy.base
plt = base.plt

R_LOW = math.cos(pi * 23.0 / 124.0) ** 2    # = cos^2(115π/620)
R_HIGH = math.cos(pi * 117.0 / 620.0) ** 2
DELTA_R_PUBLISHED = 0.009352736442513998     # R_mass_beat_two_arm_pre_v1 公表値
J_MAX = 1240
PRED_PERIOD = 620
PRED_BEAT_BIN = 2


def main() -> None:
    # ---- P1: 厳密恒等式（純算術）----
    dR = R_LOW - R_HIGH
    ident = math.sin(58 * pi / 155) * math.sin(pi / 310)
    print(f"P1: ΔR = {dR:.15f} / sin(58π/155)·sin(π/310) = {ident:.15f} "
          f"差 {abs(dR - ident):.2e}")
    assert abs(dR - ident) < 1e-15
    print(f"    公表値との差 {abs(dR - DELTA_R_PUBLISHED):.2e}（アンカー）")
    assert abs(dR - DELTA_R_PUBLISHED) < 1e-12

    # ---- 状態構成（逆算探索、前進無変更）----
    params = base.Params(high_n=63, recursive_collision_count=J_MAX)
    sp = base.build_source_params(params)
    a_t, b_t, _ = search.make_unit_templates(sp)
    arms = {}
    for label, target in (("root124", R_LOW), ("root620", R_HIGH)):
        res = search.search_initial_b_amplitude(target, a_t, b_t, sp, tolerance=1.0e-15)
        arms[label] = (a_t.copy(), res.initial_b_amplitude * b_t)

    # ---- 二腕発展と重なり O(j) ----
    (a1, b1), (a2, b2) = arms["root124"], arms["root620"]
    n1 = float(np.vdot(a1, a1).real + np.vdot(b1, b1).real)
    n2 = float(np.vdot(a2, a2).real + np.vdot(b2, b2).real)
    norm = math.sqrt(n1 * n2)
    O = np.zeros(J_MAX + 1, dtype=complex)
    O[0] = (np.vdot(a1, a2) + np.vdot(b1, b2)) / norm
    for j in range(1, J_MAX + 1):
        th1 = toy.theta_from_ab(a1, b1, sp).theta
        a1, b1 = toy.rotate_ab(a1, b1, th1)
        th2 = toy.theta_from_ab(a2, b2, sp).theta
        a2, b2 = toy.rotate_ab(a2, b2, th2)
        O[j] = (np.vdot(a1, a2) + np.vdot(b1, b2)) / norm

    # ---- P2: 厳密周期 620 ----
    per_err = float(np.max(np.abs(O[PRED_PERIOD:] - O[: J_MAX + 1 - PRED_PERIOD])))
    print(f"P2: |O(j+620)-O(j)| 最大 = {per_err:.2e}  {'PASS' if per_err <= 1e-10 else 'FAIL'}")

    # ---- P3: ビート周波数（ストロボFFT）----
    seg = np.real(O[:J_MAX]) - np.mean(np.real(O[:J_MAX]))
    mag = np.abs(np.fft.fft(seg))[: J_MAX // 2]
    peak = int(np.argmax(mag[1:]) + 1)
    print(f"P3: FFTピークビン = {peak}（予言 {PRED_BEAT_BIN}）"
          f"{'PASS' if peak == PRED_BEAT_BIN else 'FAIL'}")
    print(f"    → ビート角 = 2π·{peak}/1240 = π/310 = 変位1単位 = 電荷ステップ角")

    # ---- P4: 深いヌル ----
    min_abs = float(np.min(np.abs(O)))
    argmin = int(np.argmin(np.abs(O)))
    print(f"P4: min|O| = {min_abs:.2e} (j={argmin})  {'PASS' if min_abs <= 1e-10 else 'FAIL'}")

    ok = per_err <= 1e-10 and peak == PRED_BEAT_BIN and min_abs <= 1e-10

    # ---- 保存 ----
    with (HERE / "paperD_mass_beat_rows_v1.csv").open("w", encoding="utf-8", newline="") as h:
        w = csv.writer(h); w.writerow(["j", "reO", "imO", "absO"])
        for j in range(J_MAX + 1):
            w.writerow([j, O[j].real, O[j].imag, abs(O[j])])

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.4), constrained_layout=True)
    axes[0].plot(np.abs(O), lw=0.9)
    axes[0].set_title("|O(j)| two-arm envelope: exact period 620")
    axes[0].set_xlabel("collision j"); axes[0].set_ylabel("|<psi124|psi620>|")
    axes[1].plot(mag, lw=0.9)
    axes[1].axvline(peak, color="tab:red", ls=":", label=f"peak bin {peak}")
    axes[1].set_title("FFT of Re O(j): beat = 2 bins = charge step")
    axes[1].set_xlabel("bin (window 1240)"); axes[1].legend()
    for ext in ("png", "svg"):
        fig.savefig(HERE / f"paperD_mass_beat_v1.{ext}", dpi=160)
    plt.close(fig)

    payload = {
        "experiment": "paperD_mass_beat_alpha_pair_v1",
        "core_runner": {"path": "run_ab_invariant_theta_toy_v1.py",
                        "sha256": toy.sha256(HERE.parent / "run_ab_invariant_theta_toy_v1.py")},
        "anchor": {"delta_R_published": DELTA_R_PUBLISHED, "delta_R_here": dR},
        "P1_identity": {"delta_R": dR, "sin_factorization": ident,
                        "error": abs(dR - ident)},
        "P2_period620_max_err": per_err,
        "P3_beat_bin": {"measured": peak, "predicted": PRED_BEAT_BIN},
        "P4_min_abs_O": {"value": min_abs, "at_j": argmin},
        "all_pass": bool(ok),
        "conclusion": (
            "二α根状態間の二腕ビートは厳密周期620・ビート角π/310（FFTビン2）で、"
            "電荷ランニングの1ステップ角と同一。ΔR は sin(58π/155)·sin(π/310) に"
            "厳密に因数分解される——質量分裂型の読みと電荷ステップは同じ一つの数の"
            "相補的な面である（二文法の同数性）"),
    }
    (HERE / "paperD_mass_beat_result_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"\nall_pass = {ok} / saved: paperD_mass_beat_result_v1.json")


if __name__ == "__main__":
    main()
