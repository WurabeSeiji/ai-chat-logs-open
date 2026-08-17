#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""1 步の中を開く——衝突は細分できるか

根拠
----
collision_step_exact は積分器ではなく **閉形式の厳密解** である。正本
（run_ignition_fate_exact_v3.py 冒頭）に明記されている:

    頂点流 da/dτ = −2s·b, db/dτ = +2s·a,  s = Im(b̄a)
    は s を厳密に保存する（ds/dτ = 0）。ゆえに流れは各格子点で角速度 2s の
    厳密な回転であり、閉形式 a' = cosφ·a − sinφ·b, φ = 2R·Im(b̄a) で解ける。

つまり 1 步は「時刻 τ=1 における流れの厳密解」であって、原子的な事象ではない。
生成子 s が流れに沿って保存するなら、**τ=t（0<t<1）の厳密解も同じ式で
書ける**（回転角を φ·t にするだけ）。弾性部 rotate_ab(θ) も同型の回転。

したがって τ の細分は原理的に可能なはずである。ただし 1 步は
  Rot_φ ∘ Rot_θ
の合成であり、θ と φ を毎回読み直すと、K 分割の合成は 1 步と一致するとは
限らない（分解の順序誤差）。一致しないとしても K→∞ で収束するなら、
**連続極限が存在する＝τ は本当に細分できる**ということになる。

検定（実行前に固定）
--------------------
  S1 生成子の保存: Im(b̄a) が Rot_φ で不変か（機械精度）。
  S2 θ の保存: θ が Rot_θ で不変か、Rot_φ で不変か。
  S3 群性: t=1/K を K 回で t=1 を再現するか。K を倍にしたときの差の減り方。
     収束すれば連続極限が存在する。
  S4 步の中の Δθ: K=64 で 1 步の中を覗く。Δθ が步の中で何度動くか、
     巻いているか（|Δθ| の総変化 > 180° なら 1步刻みでは巻き数が見えない）。
  S5 別解像度の軌道: K=1 と K=64 で同じ τ 区間を走らせ、支配周期を比べる。
     5.376 步が K=64 でも残れば実在、消えればエイリアス。

使い方: python3 probe_substep_v1.py
出力  : result_substep_v1.json
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


_uni = _load("uni_ss", UNI / "unified_interaction_v1.py")
_cr0 = _load("cr0_ss", HERE / "run_cr0_control_no_theta_v2.py")
toy = _uni.two_body_base if False else None
_ex_toy = _load("cr1_ss", HERE / "run_cr1_kinetic_feedback_v1.py").toy

DEG_A, DEG_B = -30.0, +30.0
PACK_A, PACK_B = tuple(range(1, 18)), tuple(range(1, 4))


def collision_step_frac(a, b, sp, t):
    """collision_step_exact の時刻 t 版（t=1 で正本と一致）。

    正本 collision_step_exact と同じ順序・同じ量を使い、回転角だけ t 倍する。
    弾性部 θ→θ·t、非弾性部 φ→φ·t。φ は正本どおり **θ 回転後** の a,b から作る。
    """
    ro = _ex_toy.theta_from_ab(a, b, sp)
    a, b = _ex_toy.rotate_ab(a, b, ro.theta * t)
    r = float(ro.reflection_rate)
    if r > 0.0:
        phi = 2.0 * r * np.imag(np.conj(b) * a) * t
        c, s_ = np.cos(phi), np.sin(phi)
        a, b = c * a - s_ * b, s_ * a + c * b
    return a, b, ro


def make_ab(sp, slope, icept, tag):
    case = _uni.two_body_base.explicit_packet_case(
        mode=tag, packet_a=PACK_A, packet_b=PACK_B,
        packet_a_shift=_cr0.shift_for_deg(DEG_A, slope, icept),
        packet_b_shift=_cr0.shift_for_deg(DEG_B, slope, icept))
    a = _uni.two_body_base.make_case_state(sp, case, "A", hair_enabled=True)
    b = _uni.two_body_base.make_case_state(sp, case, "B", hair_enabled=True)
    return a, b


def dtheta_deg(a, b, n_chi, n_eta):
    pa, _ = _cr0.circle_position(a, n_chi, n_eta)
    pb, _ = _cr0.circle_position(b, n_chi, n_eta)
    return float(np.degrees(np.angle(np.exp(1j * (pa - pb)))))


def spec_top(x, top=3):
    y = np.asarray(x, float); y = y - y.mean()
    Y = np.abs(np.fft.rfft(y * np.hanning(len(y)))); f = np.fft.rfftfreq(len(y))
    Y[0] = 0
    tot = float(np.sum(Y ** 2))
    idx = np.argsort(Y)[::-1]
    pk = [i for i in idx if 0 < i < len(Y) - 1 and Y[i] > Y[i - 1]
          and Y[i] > Y[i + 1]][:top]
    return [(float(1 / f[i]), float(100 * Y[i] ** 2 / tot)) for i in pk]


def main() -> None:
    t0 = time.time()
    sp = _uni.two_body_base.build_source_params(
        _uni.two_body_base.Params(high_n=63, recursive_collision_count=200))
    n_chi, n_eta = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = _cr0.calibrate_shift(sp, n_chi, n_eta)
    out = {}

    # ---- S1/S2 生成子と θ の保存 --------------------------------------
    print("【S1/S2】生成子 Im(b̄a) と θ の保存")
    a, b = make_ab(sp, slope, icept, "ss_gen")
    ro = _ex_toy.theta_from_ab(a, b, sp)
    a1, b1 = _ex_toy.rotate_ab(a, b, ro.theta)
    th_after_theta = _ex_toy.theta_from_ab(a1, b1, sp).theta
    s_before = np.imag(np.conj(b1) * a1)
    r = float(ro.reflection_rate)
    phi = 2.0 * r * s_before
    c, s_ = np.cos(phi), np.sin(phi)
    a2, b2 = c * a1 - s_ * b1, s_ * a1 + c * b1
    s_after = np.imag(np.conj(b2) * a2)
    th_after_phi = _ex_toy.theta_from_ab(a2, b2, sp).theta
    ds = float(np.max(np.abs(s_after - s_before)))
    print(f"  Rot_φ による Im(b̄a) の最大変化: {ds:.3e}  "
          f"（0 なら生成子は厳密保存＝流れは厳密回転）")
    print(f"  θ: 初期 {ro.theta:.12f}  Rot_θ 後 {th_after_theta:.12f}  "
          f"変化 {abs(th_after_theta-ro.theta):.3e}")
    print(f"  θ: Rot_φ 後 {th_after_phi:.12f}  "
          f"変化 {abs(th_after_phi-ro.theta):.3e}")
    out["S1_gen_drift"] = ds
    out["S2_theta"] = {"theta0": ro.theta, "after_rot_theta": th_after_theta,
                       "after_rot_phi": th_after_phi}

    # ---- S3 群性・収束 ------------------------------------------------
    print("\n【S3】t=1/K を K 回 vs 正本 1 步（同一の初期状態から）")
    a0, b0 = make_ab(sp, slope, icept, "ss_group")
    aF, bF, _ = _uni.collision_step_exact(a0, b0, sp)
    prev = None
    rows = []
    for K in (1, 2, 4, 8, 16, 32, 64, 128):
        a, b = a0.copy(), b0.copy()
        for _ in range(K):
            a, b, _ = collision_step_frac(a, b, sp, 1.0 / K)
        dif = float(np.max(np.abs(a - aF)) + np.max(np.abs(b - bF)))
        conv = "" if prev is None else f"  前段との差 {abs(prev-dif)/max(dif,1e-300):7.3f}"
        # K 分割どうしの差（収束の見方）
        rows.append((K, dif))
        print(f"  K={K:4d}: 正本との最大差 {dif:.6e}{conv}")
        prev = dif
    # K どうしの差（連続極限の存在）
    print("  K を倍にしたときの解どうしの差（連続極限があれば単調に減る）")
    states = {}
    for K in (8, 16, 32, 64, 128, 256):
        a, b = a0.copy(), b0.copy()
        for _ in range(K):
            a, b, _ = collision_step_frac(a, b, sp, 1.0 / K)
        states[K] = (a, b)
    Ks = [8, 16, 32, 64, 128]
    for K in Ks:
        d = float(np.max(np.abs(states[K][0] - states[2 * K][0])))
        print(f"   ||ψ_K − ψ_2K||∞  K={K:4d}: {d:.6e}")
        out[f"S3_conv_K{K}"] = d
    out["S3_vs_exact"] = {str(k): v for k, v in rows}

    # ---- S4 步の中の Δθ ------------------------------------------------
    print("\n【S4】1 步の中の Δθ（K=64 で覗く・並進なし）")
    a, b = make_ab(sp, slope, icept, "ss_inside")
    K = 64
    inside = []
    for step in range(6):
        seq = [dtheta_deg(a, b, n_chi, n_eta)]
        for _ in range(K):
            a, b, _ = collision_step_frac(a, b, sp, 1.0 / K)
            seq.append(dtheta_deg(a, b, n_chi, n_eta))
        u = np.degrees(np.unwrap(np.radians(np.array(seq))))
        total = float(u[-1] - u[0])
        pathlen = float(np.sum(np.abs(np.diff(u))))
        inside.append({"step": step, "start": seq[0], "end": seq[-1],
                       "net_deg": total, "path_deg": pathlen})
        print(f"  步{step}: Δθ {seq[0]:+8.3f}° → {seq[-1]:+8.3f}°  "
              f"正味 {total:+9.3f}°  経路長 {pathlen:9.3f}°  "
              f"{'巻きあり' if abs(total) > 180 else ''}")
    out["S4_inside"] = inside

    # ---- S5 K=1 と K=64 の軌道比較 ------------------------------------
    print("\n【S5】同じ τ 区間 [0,300] を K=1 と K=64 で走らせた軌道")
    res = {}
    for K in (1, 8, 64):
        a, b = make_ab(sp, slope, icept, f"ss_traj{K}")
        ser = []
        for _ in range(300 * K):
            ser.append(dtheta_deg(a, b, n_chi, n_eta))
            a, b, _ = collision_step_frac(a, b, sp, 1.0 / K)
        ser = np.array(ser)
        # τ 単位の周期に直す
        sp_ = [(p / K, s) for p, s in spec_top(ser)]
        u = np.degrees(np.unwrap(np.radians(ser)))
        print(f"  K={K:3d}: Δθ範囲[{ser.min():+8.3f},{ser.max():+8.3f}]°  "
              f"正味回転 {(u[-1]-u[0])/360.0:+9.4f} 周  "
              f"支配周期[τ] " + "  ".join(f"{p:7.4f}({s:5.2f}%)" for p, s in sp_))
        res[K] = {"range": [float(ser.min()), float(ser.max())],
                  "turns": float((u[-1] - u[0]) / 360.0),
                  "spec_tau": sp_}
        # K=1 の標本点だけ取り出して比較
        if K > 1:
            sub = ser[::K]
            base = res[1]["_ser"]
            n = min(len(sub), len(base))
            dif = float(np.max(np.abs(np.degrees(np.angle(
                np.exp(1j * np.radians(sub[:n] - base[:n])))))))
            print(f"        K=1 の標本点との最大差 {dif:.6e}°")
            res[K]["max_diff_at_integer_tau"] = dif
        else:
            res[1]["_ser"] = ser
    res[1].pop("_ser", None)
    out["S5"] = res

    out["meta"] = {"experiment": "substep_probe_v1",
                   "date": time.strftime("%Y-%m-%d %H:%M:%S"),
                   "n_chi": n_chi, "n_eta": n_eta}
    p = HERE / "result_substep_v1.json"
    p.write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float),
                 encoding="utf-8")
    print(f"\n保存 {p.name}  ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
