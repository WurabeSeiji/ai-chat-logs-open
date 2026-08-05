#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""予備実験 v4: quadrature座標での運動学（M3'）——静止＝時計同期、運動＝離調

計器の導出（設計、実行前に固定）:
    読出し位置 = 固有時計で復調した各方向成分 r_k(t) = ⟨d_k, Z⊥(t)⟩ e^{-iφ(t)}。
    静止 = r が定数（時計に同期した成分）。
    運動 = r の偏角が等速で回る（時計に対する離調 δω がそのまま速度）。
ソフト期待（探索・記述）:
    K1 無摂動基線: r(t) の偏角ドリフトは小さい（基線静止）。
    K2 平面方向の静的キック: 変位が保持される（減衰時定数を記録）——静止系の存在。
    K3 別固有対方向のキック: 復調位置がその対の離調レートで等速ドリフト——運動＝離調。
        整合検定: ドリフトレート ≈ 基線で実測した当該対の周波数 − 時計。

N=8（空間が存在する側・計算軽量）を使用。バリオン領域（N=4,5）は対象外。

使い方: python3 run_pre_kinematics_quadrature_v4.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec1 = importlib.util.spec_from_file_location("pre1v4", HERE / "run_pre_2plus1_structure_v1.py")
pre1 = importlib.util.module_from_spec(spec1)
sys.modules[spec1.name] = pre1
spec1.loader.exec_module(pre1)
abl = pre1.abl

N = 8
T_SETTLE = 6000          # 準安定への整定
T_MEAS = 8000            # 測定区間
EV = 5
KICK_EPS = 1e-3          # キック振幅（状態ノルム比）


def run_with_kick(kick_dir=None):
    """整定→（キック）→測定。測定区間の Z サンプルと親平面射影を返す。"""
    sys_lr, v, _, _, _, p, q, Z, wp = abl.build_init(N, True)
    for _ in range(T_SETTLE):
        Z, wp = abl.evolve(sys_lr, Z, wp)
    if kick_dir is not None:
        Z = Z + KICK_EPS * np.linalg.norm(Z) * kick_dir
    samples = []
    for t in range(T_MEAS):
        Z, wp = abl.evolve(sys_lr, Z, wp)
        if t % EV == 0:
            samples.append(Z.copy())
    return np.array(samples), p, q


def analyze(S, p, q, dirs):
    ns = S.shape[0]
    Sp = np.array([z - p * (p @ z) - q * (q @ z) for z in S])
    phi = np.unwrap(np.angle((S @ p) + 1j * (S @ q)))
    om_clock = float(np.polyfit(np.arange(ns), phi, 1)[0])
    res = {"om_clock": om_clock}
    for name, d in dirs.items():
        c = Sp @ np.conj(d)
        r = c * np.exp(-1j * phi)
        ang = np.unwrap(np.angle(r))
        drift = float(np.polyfit(np.arange(ns), ang, 1)[0])       # rad/サンプル
        res[name] = {"mean_abs": float(np.mean(np.abs(r))),
                     "std_abs": float(np.std(np.abs(r))),
                     "angle_drift_per_sample": drift,
                     "drift_over_clock": drift / om_clock}
    return res


def main() -> None:
    t0 = time.time()
    # ---- 基線走行と方向の抽出 ----
    S0, p, q = run_with_kick(None)
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
    d_plane1 = cdir(0)      # 最上位平面
    d_plane2 = cdir(2)      # 第2固有対
    d_plane3 = cdir(4)
    dirs = {"plane1": d_plane1, "plane2": d_plane2, "plane3": d_plane3}

    base = analyze(S0, p, q, dirs)
    print(f"[K1] 基線（N={N}, 測定{T_MEAS}step）: ω_clock={base['om_clock']:+.5f}/サンプル")
    for k in dirs:
        b = base[k]
        print(f"     {k}: |r|={b['mean_abs']:.4f}±{b['std_abs']:.4f} "
              f"偏角ドリフト/時計={b['drift_over_clock']:+.5f}")

    # ---- K2 平面1方向の静的キック ----
    S2, _, _ = run_with_kick(d_plane1)
    k2 = analyze(S2, p, q, dirs)
    dr2 = k2["plane1"]
    disp2 = dr2["mean_abs"] - base["plane1"]["mean_abs"]
    print(f"[K2] 平面1キック: |r|={dr2['mean_abs']:.4f}（基線比 {disp2:+.4f}） "
          f"ドリフト/時計={dr2['drift_over_clock']:+.5f}")

    # ---- K3 平面2方向のキック（運動＝離調の検定） ----
    S3, _, _ = run_with_kick(d_plane2)
    k3 = analyze(S3, p, q, dirs)
    dr3 = k3["plane2"]
    pred = base["plane2"]["drift_over_clock"]
    meas = dr3["drift_over_clock"]
    print(f"[K3] 平面2キック: |r|={dr3['mean_abs']:.4f}（基線 {base['plane2']['mean_abs']:.4f}） ")
    print(f"     復調ドリフト/時計: キック後 {meas:+.5f} vs 基線の同方向 {pred:+.5f}"
          f"（運動＝離調なら一致）")

    out = {"N": N, "T_SETTLE": T_SETTLE, "T_MEAS": T_MEAS, "KICK_EPS": KICK_EPS,
           "sv_rel_top8": (sv / sv[0])[:8].tolist(),
           "K1_base": base, "K2_plane1_kick": k2, "K3_plane2_kick": k3,
           "runtime_sec": time.time() - t0}
    (HERE / "pre_kinematics_quadrature_result_v4.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
