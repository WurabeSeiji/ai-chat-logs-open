#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CR2: 倍音本数だけを変えた対照 —— {1,2} と {1,2,3}

倍音は位置を決めるのではなく、**位置という概念を成立させている**。
単一倍音 M=1 では |ψ|² が周波数 2k しか持たず、第1円周モーメントが恒等的に
ゼロになって位置が定義できない。M を増やすと位置の確度が

    |z| = 1 − 1/M

で立ち上がる（実測: M=2→0.500000、M=3→0.666667、M=17→0.941176）。

本実験は A・B ともに同じ倍音集合を与え、**本数だけを 2 と 3 で変える**。
他は CR1 と完全に同一（初期位置 ∓30°・κ=1−r・ω₀=π/72・T=40000）。
位置の確度 0.5 と 0.667 の差だけが動力学に何をするかを見る。

出力: cr2_harmonic_compare_data_v1.json
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


_uni = _load("uni_cr2", UNI / "unified_interaction_v1.py")
K = _load("kin_cr2", UNI / "unified_kinetic_v1.py")
_cr0 = _load("cr0_cr2", HERE / "run_cr0_control_no_theta_v2.py")
_cr1 = _load("cr1_cr2", HERE / "run_cr1_kinetic_feedback_v1.py")

CASES = [(1, 2), (1, 2, 3)]
DEG_A, DEG_B = -30.0, +30.0
T_STEPS = 40000
EVERY = 20
PROF_EVERY = 400
PROF_BINS = 128
OMEGA0 = np.pi / 72.0


def run_case(packet):
    sp = _uni.two_body_base.build_source_params(
        _uni.two_body_base.Params(high_n=63, recursive_collision_count=200))
    n_chi, n_eta = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = _cr0.calibrate_shift(sp, n_chi, n_eta)

    case = _uni.two_body_base.explicit_packet_case(
        mode="cr2", packet_a=packet, packet_b=packet,
        packet_a_shift=_cr0.shift_for_deg(DEG_A, slope, icept),
        packet_b_shift=_cr0.shift_for_deg(DEG_B, slope, icept))
    a = _uni.two_body_base.make_case_state(sp, case, "A", hair_enabled=True)
    b = _uni.two_body_base.make_case_state(sp, case, "B", hair_enabled=True)

    _, z0a = _cr0.circle_position(a, n_chi, n_eta)
    _, z0b = _cr0.circle_position(b, n_chi, n_eta)

    omega, v = OMEGA0, 0.0
    frames, profiles = [], []

    for s in range(T_STEPS):
        pa, za = _cr0.circle_position(a, n_chi, n_eta)
        pb, zb = _cr0.circle_position(b, n_chi, n_eta)
        chi = float(np.angle(np.exp(1j * (pa - pb))))
        r_now = float(_cr1.toy.theta_from_ab(a, b, sp).reflection_rate)
        acc = -4.0 * np.sin(omega / 2.0) ** 2 * chi
        v += acc
        omega += (1.0 - r_now) * acc          # κ = 1−r

        if s % EVERY == 0:
            RpA, m2A, _, _ = _cr1.cone_components(a, n_chi, n_eta)
            RpB, m2B, _, _ = _cr1.cone_components(b, n_chi, n_eta)
            Z = complex(np.sum(a * a) + np.sum(b * b))
            frames.append([
                s, round(pa, 6), round(pb, 6),
                round(_cr0.participation_ratio(a, n_chi, n_eta), 3),
                round(_cr0.participation_ratio(b, n_chi, n_eta), 3),
                float(f"{RpA:.6g}"), float(f"{RpB:.6g}"),
                round(chi, 6), round(2.0 * np.sin(abs(chi) / 2.0), 6),
                float(f"{Z.real:.4g}"), float(f"{Z.imag:.4g}"),
                round(omega, 8), float(f"{v:.6g}"), round(r_now, 6),
                round(za, 5), round(zb, 5)])

        if s % PROF_EVERY == 0:
            def prof(psi):
                P = np.sum(np.abs(psi.reshape(n_chi, n_eta)) ** 2, axis=1)
                P = P.reshape(PROF_BINS, -1).sum(axis=1)
                return [float(f"{x:.4g}") for x in (P / P.max())]
            profiles.append({"t": s, "A": prof(a), "B": prof(b)})

        a = K.k_translate_flat(a, -v, n_chi, n_eta)
        a, b, _ = _uni.collision_step_exact(a, b, sp)

    return {"packet": list(packet), "M": len(packet),
            "z0": [z0a, z0b], "frames": frames, "profiles": profiles}


def spectrum(x, top=6):
    y = np.asarray(x, float) - np.mean(x)
    n = len(y)
    Y = np.abs(np.fft.rfft(y * np.hanning(n)))
    f = np.fft.rfftfreq(n)
    Y[0] = 0.0
    tot = float(np.sum(Y ** 2))
    idx = np.argsort(Y)[::-1]
    pk = [i for i in idx if 0 < i < len(Y) - 1 and Y[i] > Y[i - 1] and Y[i] > Y[i + 1]][:top]
    return [{"period": float(1 / f[i]), "rel": float(Y[i] / Y[pk[0]]),
             "share": float(100 * Y[i] ** 2 / tot)} for i in pk]


def main() -> None:
    t0 = time.time()
    runs = []
    for pk in CASES:
        print(f"走行中: 倍音 {pk} …", end="", flush=True)
        r = run_case(pk)
        F = np.array(r["frames"], float)
        chi = F[:, 7]
        r["spec_chi"] = spectrum(chi)
        r["spec_Rp"] = spectrum(F[:, 5])
        runs.append(r)
        print(f" 完了  |z|={r['z0'][0]:.6f}  "
              f"Δθ範囲=[{np.degrees(chi).min():+.2f},{np.degrees(chi).max():+.2f}]°")

    out = {
        "experiment": "cr2_harmonic_compare_v1",
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "config": {"deg_a": DEG_A, "deg_b": DEG_B, "T": T_STEPS, "every": EVERY,
                   "prof_every": PROF_EVERY, "prof_bins": PROF_BINS,
                   "omega0": OMEGA0, "kappa": "1-r"},
        "columns": ["t", "posA", "posB", "prA", "prB", "RpA", "RpB", "chi",
                    "chord", "cloRe", "cloIm", "omega", "v", "r", "zA", "zB"],
        "runs": runs,
    }
    p = HERE / "cr2_harmonic_compare_data_v1.json"
    p.write_text(json.dumps(out, ensure_ascii=False, separators=(",", ":")),
                 encoding="utf-8")
    print(f"\n保存: {p.name}  {p.stat().st_size/1024:.0f}KB  ({time.time()-t0:.1f}s)")

    print("\n=== 比較 ===")
    print(f"  {'':22} {'M=2 {1,2}':>18} {'M=3 {1,2,3}':>18}")
    for lbl, fn in (
        ("位置の確度 |z|", lambda r: r["z0"][0]),
        ("PR 初期[セル]", lambda r: r["frames"][0][3]),
        ("PR 平均[セル]", lambda r: float(np.mean([f[3] for f in r["frames"]]))),
        ("Δθ 最小[°]", lambda r: float(np.degrees(min(f[7] for f in r["frames"])))),
        ("Δθ 最大[°]", lambda r: float(np.degrees(max(f[7] for f in r["frames"])))),
        ("R′² 変動[%]", lambda r: float(100 * np.std([f[5] for f in r["frames"]])
                                        / np.mean([f[5] for f in r["frames"]]))),
        ("|Σz²| 最大", lambda r: float(max(abs(complex(f[9], f[10]))
                                          for f in r["frames"]))),
        ("反射率 r 平均", lambda r: float(np.mean([f[13] for f in r["frames"]]))),
        ("Δθ 第1周期[步]", lambda r: r["spec_chi"][0]["period"]),
        ("第1ピーク占有[%]", lambda r: r["spec_chi"][0]["share"]),
    ):
        print(f"  {lbl:22} {fn(runs[0]):>18.6g} {fn(runs[1]):>18.6g}")


if __name__ == "__main__":
    main()
