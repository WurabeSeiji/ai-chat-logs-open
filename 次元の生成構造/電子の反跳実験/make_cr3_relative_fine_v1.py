#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CR3: 相対位相だけを 1步刻みで見る

CR1/CR2 の可視化は 20 步おきに間引いていた。ところが Δθ の周期は実測 5.5 步で、
**ナイキスト（周期 40 步）を大きく割っている**。20 步おきの標本化では
単振動がエイリアスして複雑な準周期に見える。

本データは **1步刻み**（間引きなし）で、τ 方向の共通回転を取り除いた
相対位相 Δθ とその速度 dΔθ/dτ を記録する。単振動なら位相空間 (Δθ, Δθ̇) は
閉じた楕円を描く。

構成は CR1 と同一（A: 倍音1〜17 / B: 倍音1〜3 / 初期 ∓30° / κ=1−r / ω₀=π/72）。

出力: cr3_relative_fine_data_v1.json
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


_uni = _load("uni_cr3", UNI / "unified_interaction_v1.py")
K = _load("kin_cr3", UNI / "unified_kinetic_v1.py")
_cr0 = _load("cr0_cr3", HERE / "run_cr0_control_no_theta_v2.py")
_cr1 = _load("cr1_cr3", HERE / "run_cr1_kinetic_feedback_v1.py")

PACKET_A = tuple(range(1, 18))
PACKET_B = tuple(range(1, 4))
DEG_A, DEG_B = -30.0, +30.0
T_STEPS = 4000            # 1步刻み・周期5.5步なので約727周期
OMEGA0 = np.pi / 72.0


def main() -> None:
    t0 = time.time()
    sp = _uni.two_body_base.build_source_params(
        _uni.two_body_base.Params(high_n=63, recursive_collision_count=200))
    n_chi, n_eta = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = _cr0.calibrate_shift(sp, n_chi, n_eta)
    case = _uni.two_body_base.explicit_packet_case(
        mode="cr3", packet_a=PACKET_A, packet_b=PACKET_B,
        packet_a_shift=_cr0.shift_for_deg(DEG_A, slope, icept),
        packet_b_shift=_cr0.shift_for_deg(DEG_B, slope, icept))
    a = _uni.two_body_base.make_case_state(sp, case, "A", hair_enabled=True)
    b = _uni.two_body_base.make_case_state(sp, case, "B", hair_enabled=True)

    omega, v = OMEGA0, 0.0
    chi_l, mid_l, pr_l, Rp_l, r_l, v_l, om_l, cre, cim = ([] for _ in range(9))

    for s in range(T_STEPS):
        pa, _ = _cr0.circle_position(a, n_chi, n_eta)
        pb, _ = _cr0.circle_position(b, n_chi, n_eta)
        chi = float(np.angle(np.exp(1j * (pa - pb))))
        mid = float(np.angle(np.exp(1j * pa) + np.exp(1j * pb)))  # 共通回転（重心方位）
        r_now = float(_cr1.toy.theta_from_ab(a, b, sp).reflection_rate)
        acc = -4.0 * np.sin(omega / 2.0) ** 2 * chi
        v += acc
        omega += (1.0 - r_now) * acc

        RpA, _, _, _ = _cr1.cone_components(a, n_chi, n_eta)
        Z = complex(np.sum(a * a) + np.sum(b * b))
        chi_l.append(round(chi, 8)); mid_l.append(round(mid, 6))
        pr_l.append(round(_cr0.participation_ratio(a, n_chi, n_eta), 2))
        Rp_l.append(float(f"{RpA:.6g}")); r_l.append(round(r_now, 6))
        v_l.append(float(f"{v:.6g}")); om_l.append(round(omega, 8))
        cre.append(float(f"{Z.real:.4g}")); cim.append(float(f"{Z.imag:.4g}"))

        a = K.k_translate_flat(a, -v, n_chi, n_eta)
        a, b, _ = _uni.collision_step_exact(a, b, sp)

    chi = np.array(chi_l)
    dchi = np.gradient(chi)                     # dΔθ/dτ（1步刻みなので差分＝速度）

    out = {
        "experiment": "cr3_relative_fine_v1",
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "config": {"packet_a": list(PACKET_A), "packet_b": list(PACKET_B),
                   "deg_a": DEG_A, "deg_b": DEG_B, "T": T_STEPS,
                   "every": 1, "omega0": OMEGA0, "kappa": "1-r"},
        "chi": chi_l, "dchi": [float(f"{x:.6g}") for x in dchi],
        "mid": mid_l, "pr": pr_l, "Rp": Rp_l, "r": r_l,
        "v": v_l, "omega": om_l, "cloRe": cre, "cloIm": cim,
    }
    p = HERE / "cr3_relative_fine_data_v1.json"
    p.write_text(json.dumps(out, ensure_ascii=False, separators=(",", ":")),
                 encoding="utf-8")

    # 単振動かどうかの一次診断
    y = chi - chi.mean()
    Y = np.abs(np.fft.rfft(y * np.hanning(len(y))))
    f = np.fft.rfftfreq(len(y)); Y[0] = 0
    i = int(np.argmax(Y))
    share = 100 * Y[i] ** 2 / np.sum(Y ** 2)
    print(f"保存: {p.name}  {p.stat().st_size/1024:.0f}KB  ({time.time()-t0:.1f}s)")
    print(f"  Δθ 範囲 = [{np.degrees(chi).min():+.3f}, {np.degrees(chi).max():+.3f}]°")
    print(f"  第1ピーク 周期={1/f[i]:.4f} 步   パワー占有={share:.2f}%")
    print(f"  （20步おきの標本化では周期 5.5 步がナイキスト 40 步を割り、"
          f"エイリアスして複雑に見えていた）")
    print(f"  共通回転 mid の総回転 = "
          f"{np.unwrap(np.array(mid_l))[-1]/(2*np.pi)-np.unwrap(np.array(mid_l))[0]/(2*np.pi):+.3f} 周")


if __name__ == "__main__":
    main()
