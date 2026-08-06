#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""予備実験 v5: v4 と同一シード・同一手順の再走行＋時系列保存（図用）

v4 は要約値のみ保存していた。本スクリプトは同一の走行から
quadrature 時系列（plane1/plane2/axis の復調複素座標）を保存し、
v4 要約値との一致を対照検定してから書き出す。

使い方: python3 run_pre_kinematics_series_v5.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec4 = importlib.util.spec_from_file_location("pre4s", HERE / "run_pre_kinematics_quadrature_v4.py")
pre4 = importlib.util.module_from_spec(spec4)
sys.modules[spec4.name] = pre4
spec4.loader.exec_module(pre4)
abl = pre4.abl
N, T_SETTLE, T_MEAS, EV, KICK_EPS = pre4.N, pre4.T_SETTLE, pre4.T_MEAS, pre4.EV, pre4.KICK_EPS


def main() -> None:
    t0 = time.time()
    stored = json.loads((HERE / "pre_kinematics_quadrature_result_v4.json").read_text())

    S0, p, q = pre4.run_with_kick(None)
    ns = S0.shape[0]
    Sp0 = np.array([z - p * (p @ z) - q * (q @ z) for z in S0])
    X = np.hstack([Sp0.real, Sp0.imag])
    Xc = X - X.mean(axis=0)
    _, sv, Vt = np.linalg.svd(Xc, full_matrices=False)
    M = S0.shape[1]

    def cdir(k):
        u = Vt[k]
        d = u[:M] + 1j * u[M:]
        return d / np.linalg.norm(d)
    dirs = {"plane1": cdir(0), "plane2": cdir(2), "plane3": cdir(4)}

    def series(S):
        Sp = np.array([z - p * (p @ z) - q * (q @ z) for z in S])
        phi = np.unwrap(np.angle((S @ p) + 1j * (S @ q)))
        om = float(np.polyfit(np.arange(len(S)), phi, 1)[0])
        r = {k: (Sp @ np.conj(d)) * np.exp(-1j * phi) for k, d in dirs.items()}
        return om, r

    om0, r0 = series(S0)
    S2, _, _ = pre4.run_with_kick(dirs["plane1"])
    om2, r2 = series(S2)
    S3, _, _ = pre4.run_with_kick(dirs["plane2"])
    om3, r3 = series(S3)

    # 対照検定: v4 要約値との一致
    def drift(r):
        ang = np.unwrap(np.angle(r))
        return float(np.polyfit(np.arange(len(r)), ang, 1)[0])
    checks = {
        "om_clock": (om0, stored["K1_base"]["om_clock"]),
        "base_plane2_drift": (drift(r0["plane2"]) / om0,
                                stored["K1_base"]["plane2"]["angle_drift_per_sample"]
                                / stored["K1_base"]["om_clock"]),
        "kick_plane2_drift": (drift(r3["plane2"]) / om3,
                                stored["K3_plane2_kick"]["plane2"]["angle_drift_per_sample"]
                                / stored["K3_plane2_kick"]["om_clock"]),
    }
    ok = all(abs(a - b) < 1e-6 for a, b in checks.values())
    for k, (a, b) in checks.items():
        print(f"  対照 {k}: 再走行 {a:+.6f} vs v4 {b:+.6f}")
    print(f"  対照検定一致 = {ok}")

    np.savez_compressed(
        HERE / "pre_kinematics_series_v5.npz",
        base_plane1=r0["plane1"], base_plane2=r0["plane2"], base_plane3=r0["plane3"],
        kick1_plane1=r2["plane1"], kick2_plane2=r3["plane2"],
        om_clock=np.array([om0, om2, om3]), sv_rel=(sv / sv[0])[:12])
    json.dump({"contrast_ok": bool(ok), "ns": ns,
               "checks": {k: [a, b] for k, (a, b) in checks.items()},
               "runtime_sec": time.time() - t0},
              open(HERE / "pre_kinematics_series_v5.json", "w"), ensure_ascii=False, indent=1)
    print(f"saved ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
