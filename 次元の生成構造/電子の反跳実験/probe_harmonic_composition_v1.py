#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""倍音構成が決めるもの——τ の刻み・位置・相互作用

問い（木原）: 相対位相 Δθ の τ 発展をもっと細かく読めないか。原理的に無理か。

結論: 読めない。ただし分解能でも計器でもなく、**τ の刻みが量子化されている**
      から。刻みは倍音構成が決め、その値は整数を数えた有理数しか取らない。

規則（実測で同定・全て 1e-16 で一致）
--------------------------------------
θ は AB 合成 χ スペクトルの |k|≥4 かつ偶数のパワー占有率 r = sin²θ で決まる。
搬送波 q_A=+1, q_B=−1 で倍音 n は bin n±1 に移るので、寄与は

    倍音 3        → 1/2   （bin 2 と 4 に分かれ、4 だけが該当）
    奇数倍音 ≥5   → 1     （両側波とも該当）
    倍音 1・偶数倍音 → 0

各体はパワー 1 に規格化されるので、体ごとに自分の本数で割って平均する:

    r = ( w(A)/N_A + w(B)/N_B ) / 2 ,   w(P) = Σ_{n∈P} 寄与(n)

Δθ の周期は π/θ 步（弾性回転の周期 2π/θ の半分。Δθ は A↔B 交換で偶）。

そこから出ること
----------------
  * r は整数を数えた有理数しか取らない → θ は離散 → τ の刻みは離散。
    刻みの間に値はない。連続に細分することはできない（＝量子化）。
  * 連続族 1〜M は M→∞ で r → 1/2、周期 → 4 步。倍音を増やしても
    細かくならない。偶数が増えるのと同じだけ奇数も増えるため。
  * 細かくするには偶数倍音が奇数倍音を圧倒する必要がある。
  * 位置と相互作用は独立な自由度。r を固定して位置だけ、位置を固定して
    r だけ動かせる（本スクリプト T4 の対照対で実証）。

使い方: python3 probe_harmonic_composition_v1.py
出力  : result_harmonic_composition_v1.json
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from fractions import Fraction
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


_uni = _load("uni_hc", UNI / "unified_interaction_v1.py")
_cr0 = _load("cr0_hc", HERE / "run_cr0_control_no_theta_v2.py")
toy = _load("cr1_hc", HERE / "run_cr1_kinetic_feedback_v1.py").toy

DEG_A, DEG_B = -30.0, +30.0
NEW_SERIES = (1, 2, 3, 4, 6, 8, 10, 12, 14, 16)


def contribution(n):
    """倍音 n が P_f に入れる量。"""
    if n == 3:
        return 0.5
    return 1.0 if (n % 2 == 1 and n >= 5) else 0.0


def r_predicted(pa, pb):
    w = lambda p: sum(contribution(n) for n in p)
    return (w(pa) / len(pa) + w(pb) / len(pb)) / 2.0


def make(sp, slope, icept, pa, pb, deg_b=DEG_B):
    case = _uni.two_body_base.explicit_packet_case(
        mode="hc", packet_a=tuple(pa), packet_b=tuple(pb),
        packet_a_shift=_cr0.shift_for_deg(DEG_A, slope, icept),
        packet_b_shift=_cr0.shift_for_deg(deg_b, slope, icept))
    return (_uni.two_body_base.make_case_state(sp, case, "A", hair_enabled=True),
            _uni.two_body_base.make_case_state(sp, case, "B", hair_enabled=True))


def amp_chi(psi, n_chi, n_eta):
    return np.sqrt(np.sum(np.abs(psi.reshape(n_chi, n_eta)) ** 2, axis=1))


def W15(psi, n_chi, n_eta):
    """裾の重み s=1.5 の幅［°］（probe_width_definitions_v2 で決めた計器）。"""
    d = amp_chi(psi, n_chi, n_eta) ** 1.5
    d = d / d.sum()
    return float(1.0 / np.sum(d ** 2)) * 360.0 / n_chi


def lobes(psi, n_chi, n_eta, thr=0.5):
    """振幅が最大の thr 倍以上のピーク（位置［°］, 高さ）。像が割れたかを見る。"""
    A = amp_chi(psi, n_chi, n_eta)
    A = A / A.max()
    raw = [(i * 360.0 / n_chi, A[i]) for i in range(n_chi)
           if A[i] >= thr and A[i] >= A[(i - 1) % n_chi] and A[i] >= A[(i + 1) % n_chi]]
    out = []
    for d, v in raw:
        if out and abs(d - out[-1][0]) < 3.0:
            if v > out[-1][1]:
                out[-1] = (d, v)
        else:
            out.append((d, v))
    return out


def period_of(sp, n_chi, n_eta, a, b, T=4000):
    """Δθ の支配周期を実測（衝突のみ・並進なし）。"""
    ch = []
    for _ in range(T):
        pa, _ = _cr0.circle_position(a, n_chi, n_eta)
        pb, _ = _cr0.circle_position(b, n_chi, n_eta)
        ch.append(float(np.degrees(np.angle(np.exp(1j * (pa - pb))))))
        a, b, _ = _uni.collision_step_exact(a, b, sp)
    ch = np.array(ch)
    y = ch - ch.mean()
    Y = np.abs(np.fft.rfft(y * np.hanning(len(y))))
    f = np.fft.rfftfreq(len(y)); Y[0] = 0
    i = int(np.argmax(Y))
    return (float(1 / f[i]), float(100 * Y[i] ** 2 / np.sum(Y ** 2)),
            float(np.median(np.abs(np.diff(ch)))), ch)


def main() -> None:
    t0 = time.time()
    sp = _uni.two_body_base.build_source_params(
        _uni.two_body_base.Params(high_n=63, recursive_collision_count=200))
    n_chi, n_eta = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = _cr0.calibrate_shift(sp, n_chi, n_eta)
    out = {}

    # ---- T1 規則の検証 ------------------------------------------------
    print("【T1】r の規則  r = (w(A)/N_A + w(B)/N_B)/2   w: 倍音3→1/2, 奇数≥5→1, 他→0")
    print(f"  {'A':>26} {'B':>14} {'r予測':>10} {'r実測':>11} {'差':>9} {'分数':>9}")
    T1 = []
    tests = [((1, 2, 3), (1, 2, 3)), ((1, 2, 3, 4), (1, 2, 3, 4)),
             (tuple(range(1, 6)), tuple(range(1, 6))),
             (tuple(range(1, 8)), tuple(range(1, 8))),
             (tuple(range(1, 18)), tuple(range(1, 18))),
             (tuple(range(1, 18)), (1, 2, 3)),
             ((3,), (3,)), ((3, 5), (3, 5)),
             ((1, 3, 5, 7, 9, 13, 17), (1, 3, 5, 7, 9, 13, 17)),
             (NEW_SERIES, NEW_SERIES), (NEW_SERIES, (1, 2, 3))]
    for pa, pb in tests:
        a, b = make(sp, slope, icept, pa, pb)
        r = float(toy.theta_from_ab(a, b, sp).reflection_rate)
        pr = r_predicted(pa, pb)
        # 連続 1..M のときだけ "1..M" と略す。飛び飛びの集合を略すと
        # 連続族と誤読されるので、その場合は末尾を省略した実列を出す。
        def lab_of(p, w):
            if tuple(p) == tuple(range(1, max(p) + 1)):
                return f"1..{max(p)}"
            s = str(tuple(p))
            return s if len(s) <= w else s[:w - 4] + "…)"
        sa, sb = lab_of(pa, 26), lab_of(pb, 14)
        print(f"  {sa:>26} {sb:>14} {pr:10.6f} {r:11.6f} {abs(pr-r):9.1e} "
              f"{str(Fraction(r).limit_denominator(500)):>9}")
        T1.append({"a": list(pa), "b": list(pb), "r_pred": pr, "r": r,
                   "frac": str(Fraction(r).limit_denominator(500))})
    out["T1_rule"] = T1

    # ---- T2 連続族の頭打ち ---------------------------------------------
    print("\n【T2】連続 1〜M（対称）は M→∞ で r→1/2・周期→4 步（細かくならない）")
    print(f"  {'M':>3} {'偶数':>5} {'奇数≥5':>6} {'r':>10} {'分数':>8} {'θ[°]':>8} {'π/θ[步]':>9}")
    T2 = []
    for M in [3, 4, 5, 6, 7, 9, 13, 17, 25, 33, 49]:
        p = tuple(range(1, M + 1))
        a, b = make(sp, slope, icept, p, p)
        ro = toy.theta_from_ab(a, b, sp)
        r, th = float(ro.reflection_rate), float(ro.theta)
        per = np.pi / th if th > 1e-14 else float("inf")
        ev = sum(1 for n in p if n % 2 == 0)
        od = sum(1 for n in p if n % 2 == 1 and n >= 5)
        print(f"  {M:3d} {ev:5d} {od:6d} {r:10.6f} "
              f"{str(Fraction(r).limit_denominator(500)):>8} "
              f"{np.degrees(th):8.3f} {per:9.4f}")
        T2.append({"M": M, "even": ev, "odd5": od, "r": r, "period": per})
    out["T2_contiguous"] = T2

    # ---- T3 提案系列 --------------------------------------------------
    print(f"\n【T3】提案系列 {NEW_SERIES}（奇数≥3 は倍音3 の1本のみ）")
    print(f"  {'構成':30} {'r':>10} {'分数':>8} {'θ[°]':>8} {'予測π/θ':>9} "
          f"{'実測周期':>9} {'占有%':>7} {'1步Δθ':>8} {'W1.5':>8} {'|z|':>8} {'像':>3}")
    T3 = []
    for lab, pa, pb in [("A=B=提案系列", NEW_SERIES, NEW_SERIES),
                        ("A=提案系列 B=(1,2,3)", NEW_SERIES, (1, 2, 3)),
                        ("現行 A=1..17 B=1..3", tuple(range(1, 18)), (1, 2, 3)),
                        ("A=B=1..17", tuple(range(1, 18)), tuple(range(1, 18))),
                        ("A=B=(1,2,3)", (1, 2, 3), (1, 2, 3))]:
        a, b = make(sp, slope, icept, pa, pb)
        ro = toy.theta_from_ab(a, b, sp)
        r, th = float(ro.reflection_rate), float(ro.theta)
        w = W15(a, n_chi, n_eta)
        z = _cr0.circle_position(a, n_chi, n_eta)[1]
        per, share, dmed, ch = period_of(sp, n_chi, n_eta, a, b)
        nl = len(lobes(a, n_chi, n_eta))
        print(f"  {lab:30} {r:10.6f} {str(Fraction(r).limit_denominator(500)):>8} "
              f"{np.degrees(th):8.3f} {np.pi/th:9.4f} {per:9.4f} {share:7.2f} "
              f"{dmed:8.3f} {w:8.3f} {z:8.6f} {nl:3d}")
        T3.append({"label": lab, "r": r, "theta_deg": float(np.degrees(th)),
                   "period_pred": float(np.pi / th), "period_meas": per,
                   "share": share, "dtheta_med": dmed, "W15": w, "z": z,
                   "n_lobes": nl,
                   "chi_head": [round(float(x), 3) for x in ch[:30]]})
    out["T3_series"] = T3
    print("\n  提案系列の Δθ 最初の30步［°］:")
    print("   " + " ".join(f"{x:+.1f}" for x in T3[0]["chi_head"]))

    # ---- T4 位置と相互作用は独立 ---------------------------------------
    print("\n【T4】位置と相互作用は独立な自由度（対照対）")
    print(f"  {'構成':22} {'r':>10} {'θ[°]':>8} {'π/θ':>8} {'W1.5[°]':>9} "
          f"{'|z|':>8} {'像':>3}")
    T4 = []
    for lab, p in [("(1,2,3,4)", (1, 2, 3, 4)), ("(3,2,4,6)", (3, 2, 4, 6)),
                   ("1..7", tuple(range(1, 8))),
                   ("(1,3,5,7,9,13,17)", (1, 3, 5, 7, 9, 13, 17)),
                   ("偶数のみ 2..16", tuple(range(2, 17, 2))),
                   (str(NEW_SERIES), NEW_SERIES)]:
        a, b = make(sp, slope, icept, p, p)
        ro = toy.theta_from_ab(a, b, sp)
        r, th = float(ro.reflection_rate), float(ro.theta)
        per = np.pi / th if th > 1e-14 else float("inf")
        lb = lobes(a, n_chi, n_eta)
        print(f"  {lab:22} {r:10.6f} {np.degrees(th):8.3f} {per:8.3f} "
              f"{W15(a, n_chi, n_eta):9.3f} "
              f"{_cr0.circle_position(a, n_chi, n_eta)[1]:8.6f} {len(lb):3d}")
        T4.append({"label": lab, "r": r, "period": per,
                   "W15": W15(a, n_chi, n_eta),
                   "z": _cr0.circle_position(a, n_chi, n_eta)[1],
                   "lobes": [[round(d, 2), round(float(v), 4)] for d, v in lb]})
    out["T4_independence"] = T4
    print("  ・(1,2,3,4) と (3,2,4,6): r・θ・周期が完全同一で W1.5 だけ違う")
    print("  ・1..7 と (1,3,5,7,9,13,17): W1.5 がほぼ同じで r が 2.2 倍違う")
    print("  ・偶数のみは像が 180° 間隔で 2 個・高さ比 1.000 → |z|=0（回折格子）")

    out["meta"] = {"experiment": "harmonic_composition_probe_v1",
                   "date": time.strftime("%Y-%m-%d %H:%M:%S"),
                   "n_chi": n_chi, "n_eta": n_eta,
                   "new_series": list(NEW_SERIES)}
    p_ = HERE / "result_harmonic_composition_v1.json"
    p_.write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float),
                  encoding="utf-8")
    print(f"\n保存 {p_.name}  ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
