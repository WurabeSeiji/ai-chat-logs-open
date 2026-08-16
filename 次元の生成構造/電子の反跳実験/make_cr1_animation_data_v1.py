#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CR1 の τ 発展を動画化するためのデータ生成（κ = 1−r 確定版）

出力: cr1_animation_data_v1.json
  frames[]: t, posA, posB（円周角[rad]）, prA, prB（参加率＝波束の広がり[セル]）,
            RpA, RpB（R′²）, chi（Δθ[rad]）, chord（弦 R″/R）,
            cloRe, cloIm（閉包 Σz² の実部・虚部）, omega, v, r
  profile[]: 一定間隔で χ 上のパワー分布そのもの（波形表示用・間引き）
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


_uni = _load("uni_anim", UNI / "unified_interaction_v1.py")
K = _load("kin_anim", UNI / "unified_kinetic_v1.py")
_cr0 = _load("cr0_anim", HERE / "run_cr0_control_no_theta_v2.py")
_cr1 = _load("cr1_anim", HERE / "run_cr1_kinetic_feedback_v1.py")

T_STEPS = 40000
EVERY = 20            # フレーム間引き（2000 フレーム）
PROF_EVERY = 400      # 波形そのものを保存する間隔（100 枚）
PROF_BINS = 128       # χ 512 → 128 に縮約


def main() -> None:
    t0 = time.time()
    sp = _uni.two_body_base.build_source_params(
        _uni.two_body_base.Params(high_n=63, recursive_collision_count=200))
    n_chi, n_eta = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = _cr0.calibrate_shift(sp, n_chi, n_eta)
    a, b = _cr0.make_pair(sp, _cr0.shift_for_deg(-30.0, slope, icept),
                          _cr0.shift_for_deg(+30.0, slope, icept))

    omega = np.pi / 72.0
    v = 0.0
    frames, profiles = [], []

    for s in range(T_STEPS):
        pa, _ = _cr0.circle_position(a, n_chi, n_eta)
        pb, _ = _cr0.circle_position(b, n_chi, n_eta)
        chi = float(np.angle(np.exp(1j * (pa - pb))))
        r_now = float(_cr1.toy.theta_from_ab(a, b, sp).reflection_rate)
        kappa = 1.0 - r_now                      # κ = 1−r（透過率・確定）
        acc = -4.0 * np.sin(omega / 2.0) ** 2 * chi
        v += acc
        omega += kappa * acc

        if s % EVERY == 0:
            RpA, _, _, _ = _cr1.cone_components(a, n_chi, n_eta)
            RpB, _, _, _ = _cr1.cone_components(b, n_chi, n_eta)
            Z = complex(np.sum(a * a) + np.sum(b * b))
            frames.append([
                s, round(pa, 6), round(pb, 6),
                round(_cr0.participation_ratio(a, n_chi, n_eta), 3),
                round(_cr0.participation_ratio(b, n_chi, n_eta), 3),
                float(f"{RpA:.6g}"), float(f"{RpB:.6g}"),
                round(chi, 6), round(2.0 * np.sin(abs(chi) / 2.0), 6),
                float(f"{Z.real:.4g}"), float(f"{Z.imag:.4g}"),
                round(omega, 8), float(f"{v:.6g}"), round(r_now, 6)])

        if s % PROF_EVERY == 0:
            def prof(psi):
                P = np.sum(np.abs(psi.reshape(n_chi, n_eta)) ** 2, axis=1)
                P = P.reshape(PROF_BINS, -1).sum(axis=1)
                return [float(f"{x:.4g}") for x in (P / P.max())]
            profiles.append({"t": s, "A": prof(a), "B": prof(b)})

        a = K.k_translate_flat(a, -v, n_chi, n_eta)
        a, b, _ = _uni.collision_step_exact(a, b, sp)

    out = {
        "experiment": "cr1_kinetic_feedback (kappa = 1-r)",
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "config": {"packet_a": "1..17", "packet_b": "1..3",
                   "deg_a": -30.0, "deg_b": 30.0, "T": T_STEPS,
                   "every": EVERY, "prof_every": PROF_EVERY,
                   "prof_bins": PROF_BINS,
                   "omega0": float(np.pi / 72.0), "kappa": "1-r",
                   "n_chi": n_chi, "n_eta": n_eta},
        "columns": ["t", "posA", "posB", "prA", "prB", "RpA", "RpB",
                    "chi", "chord", "cloRe", "cloIm", "omega", "v", "r"],
        "frames": frames,
        "profiles": profiles,
    }
    p = HERE / "cr1_animation_data_v1.json"
    p.write_text(json.dumps(out, ensure_ascii=False, separators=(",", ":")),
                 encoding="utf-8")
    print(f"frames={len(frames)}  profiles={len(profiles)}  "
          f"{p.stat().st_size/1024:.0f}KB  ({time.time()-t0:.1f}s)")
    F = np.array([f[1:] for f in frames], float)
    print(f"  posA[rad] {F[:,0].min():+.3f}..{F[:,0].max():+.3f}   "
          f"posB {F[:,1].min():+.3f}..{F[:,1].max():+.3f}")
    print(f"  PR   A {F[:,2].min():.1f}..{F[:,2].max():.1f}   "
          f"B {F[:,3].min():.1f}..{F[:,3].max():.1f}")
    print(f"  R′²  A {F[:,4].min():.3e}..{F[:,4].max():.3e}")
    print(f"  Δθ[°] {np.degrees(F[:,6]).min():+.2f}..{np.degrees(F[:,6]).max():+.2f}")


if __name__ == "__main__":
    main()
