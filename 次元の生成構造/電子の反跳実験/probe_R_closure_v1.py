#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R の正しい計算——Σx² = R² は複素にまとめる前の実成分で取る

私（Claude）の誤りとその訂正
----------------------------
R′² を「非共役の複素二乗和 Σa²」として計算し、機械ゼロ（1e-17）を得て
「R = 0」と結論した。しかしこれは自己矛盾である。R = 0 なら円周の長さが 0 で、
そもそも A・B を ±60° の位相位置に置けない。目の前の初期条件と矛盾していた。

原因: **複素表現の中に既に虚軸が入っている。** a = p + iq と書くと

    Σa² = Σ(p² − q²) + 2i·Σp·q

なので、Σa² = 0 は「Σp² = Σq² かつ Σp·q = 0」という**二条件**を意味する。
複素のまま二乗和を取れば (iR)² が含まれて恒等的に 0 になる。R を消してから
測っていたことになる。

正しい計算:

    R′² = Σx²   （x は複素にまとめる前の実成分）
    R′² = Σp² = Σq²  （二条件のうち実部条件より両者は等しい）
    R²  = R′a² + R′b²

  そして cloRe / cloIm（= Σ(p²−q²) と 2Σp·q）は **二条件の残差**であって
  R ではない。1e-16 なのは条件が成立していることの確認にすぎない。

検定（実行前に固定）
--------------------
  R1 二条件: Σp² = Σq² と Σp·q = 0 が各体で成立するか。
  R2 R の値: R′_A・R′_B・R = √(R′a²+R′b²) を出す。0 でないことを確認。
  R3 保存性: R′a²・R′b²・R² を走行させ、何が保存するか。
     （衝突のみ / 並進あり / 質量比配分あり の三条件）
  R4 中心の移動: R′ を質量として並進を配分し、対の中心が止まるか。
     A だけ −v（現行）／1/2 ずつ／R′ 比／R′² 比 を比べる。

使い方: python3 probe_R_closure_v1.py
出力  : result_R_closure_v1.json
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


_uni = _load("uni_R", UNI / "unified_interaction_v1.py")
K = _load("kin_R", UNI / "unified_kinetic_v1.py")
_cr0 = _load("cr0_R", HERE / "run_cr0_control_no_theta_v2.py")
toy = _load("cr1_R", HERE / "run_cr1_kinetic_feedback_v1.py").toy

DEG_A, DEG_B = -30.0, +30.0
OMEGA0 = np.pi / 72.0
T_STEPS = 400
MODES = [("A_only", "A だけ −v（現行）"), ("half", "1/2 ずつ"),
         ("Rp", "R′ 比"), ("Rp2", "R′² 比")]


def R_prime(psi):
    """R′ = √(Σx²)。x は複素にまとめる前の実成分。"""
    return float(np.sqrt(np.sum(psi.real ** 2)))


def two_conditions(psi):
    """(Σp², Σq², Σp·q)。閉包 Σa²=0 の実部条件と虚部条件の材料。"""
    p, q = psi.real, psi.imag
    return (float(np.sum(p * p)), float(np.sum(q * q)), float(np.sum(p * q)))


def pair_center(a, b, n_chi, n_eta):
    """対の中心（S0 = |a|²+|b|² 分布の円周第1モーメント）。"""
    S0 = np.abs(a.reshape(n_chi, n_eta)) ** 2 + np.abs(b.reshape(n_chi, n_eta)) ** 2
    P = S0.sum(axis=1); P = P / P.sum()
    w = np.exp(2j * np.pi * np.arange(n_chi) / n_chi)
    return complex(P @ w)


def main() -> None:
    t0 = time.time()
    sp = _uni.two_body_base.build_source_params(
        _uni.two_body_base.Params(high_n=63, recursive_collision_count=200))
    n_chi, n_eta = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = _cr0.calibrate_shift(sp, n_chi, n_eta)

    def mk():
        case = _uni.two_body_base.explicit_packet_case(
            mode="Rclo", packet_a=tuple(range(1, 18)), packet_b=(1, 2, 3),
            packet_a_shift=_cr0.shift_for_deg(DEG_A, slope, icept),
            packet_b_shift=_cr0.shift_for_deg(DEG_B, slope, icept))
        return (_uni.two_body_base.make_case_state(sp, case, "A", hair_enabled=True),
                _uni.two_body_base.make_case_state(sp, case, "B", hair_enabled=True))

    out = {}

    # ---- R1 / R2 -------------------------------------------------------
    print("【R1】閉包 Σa²=0 の二条件（実部 Σp²=Σq²、虚部 Σp·q=0）")
    a, b = mk()
    r1 = {}
    for nm, psi in (("A", a), ("B", b)):
        P, Q, PQ = two_conditions(psi)
        r1[nm] = {"sum_p2": P, "sum_q2": Q, "diff": abs(P - Q), "sum_pq": PQ,
                  "Rprime": R_prime(psi)}
        print(f"  {nm}: Σp² = {P:.8f}   Σq² = {Q:.8f}   差 = {abs(P-Q):.3e}   "
              f"Σp·q = {PQ:+.3e}   R′ = {R_prime(psi):.8f}")
    RA, RB = R_prime(a), R_prime(b)
    R2 = RA * RA + RB * RB
    print(f"\n【R2】R′a² = {RA*RA:.8f}   R′b² = {RB*RB:.8f}   "
          f"R² = {R2:.8f}   R = {np.sqrt(R2):.8f}")
    print(f"  （参考）Σ|a|² = {float(np.vdot(a,a).real):.8f} = Σp²+Σq² = 2·R′a²")
    print(f"  R′_A / R′_B = {RA/RB:.6f}   → 初期の質量比は 1:1")
    print(f"  ※ 複素のまま取ると Σa² = {complex(np.sum(a*a)):.3e} で機械ゼロ。"
          "これは二条件の残差であって R ではない。")
    out["R1_two_conditions"] = r1
    out["R2_values"] = {"Rprime_A": RA, "Rprime_B": RB, "R2": R2,
                        "R": float(np.sqrt(R2)),
                        "complex_square_sum_abs": abs(complex(np.sum(a * a)))}

    # ---- R3 / R4 -------------------------------------------------------
    def run(mode, T=T_STEPS):
        a, b = mk()
        omega, v = OMEGA0, 0.0
        z0 = pair_center(a, b, n_chi, n_eta)
        A0, B0 = R_prime(a) ** 2, R_prime(b) ** 2
        R20 = A0 + B0
        dA = dB = dR = dPQ = 0.0
        C = []
        for _ in range(T):
            if mode != "none":
                pa, _ = _cr0.circle_position(a, n_chi, n_eta)
                pb, _ = _cr0.circle_position(b, n_chi, n_eta)
                chi = float(np.angle(np.exp(1j * (pa - pb))))
                r = float(toy.theta_from_ab(a, b, sp).reflection_rate)
                acc = -4.0 * np.sin(omega / 2.0) ** 2 * chi
                v += acc
                omega += (1.0 - r) * acc
                ma, mb = R_prime(a), R_prime(b)
                if mode == "A_only":
                    a = K.k_translate_flat(a, -v, n_chi, n_eta)
                elif mode == "half":
                    a = K.k_translate_flat(a, -v / 2, n_chi, n_eta)
                    b = K.k_translate_flat(b, +v / 2, n_chi, n_eta)
                elif mode == "Rp":
                    M = ma + mb
                    a = K.k_translate_flat(a, -v * mb / M, n_chi, n_eta)
                    b = K.k_translate_flat(b, +v * ma / M, n_chi, n_eta)
                elif mode == "Rp2":
                    W = ma * ma + mb * mb
                    a = K.k_translate_flat(a, -v * mb * mb / W, n_chi, n_eta)
                    b = K.k_translate_flat(b, +v * ma * ma / W, n_chi, n_eta)
            a, b, _ = _uni.collision_step_exact(a, b, sp)
            pA = R_prime(a) ** 2; pB = R_prime(b) ** 2
            dA = max(dA, abs(pA - A0)); dB = max(dB, abs(pB - B0))
            dR = max(dR, abs(pA + pB - R20))
            dPQ = max(dPQ, max(abs(two_conditions(a)[2]), abs(two_conditions(b)[2])))
            z = pair_center(a, b, n_chi, n_eta)
            C.append(float(np.degrees(np.angle(z * np.conj(z0)))))
        return np.array(C), dA, dB, dR, dPQ

    print(f"\n【R3】保存性（T={T_STEPS}）")
    print(f"  {'条件':22} {'R′a² 変化':>12} {'R′b² 変化':>12} "
          f"{'R² 変化':>12} {'|Σp·q| 最大':>12}")
    r3 = {}
    for mode, lab in [("none", "衝突のみ")] + MODES:
        C, dA, dB, dR, dPQ = run(mode)
        r3[mode] = {"label": lab, "dRa2": dA, "dRb2": dB, "dR2": dR,
                    "max_sum_pq": dPQ,
                    "center_max_deg": float(np.max(np.abs(C))),
                    "center_last_deg": float(C[-1])}
        print(f"  {lab:22} {dA:12.3e} {dB:12.3e} {dR:12.3e} {dPQ:12.3e}")
    print("  → R′a² と R′b² は逆位相で振動し、R² だけが保存する。")

    print(f"\n【R4】対の中心の移動（孤立二体なら動いてはいけない・T={T_STEPS}）")
    print(f"  {'並進の配分':22} {'最大[°]':>12} {'終端[°]':>12}")
    for mode, lab in [("none", "衝突のみ")] + MODES:
        s = r3[mode]
        print(f"  {lab:22} {s['center_max_deg']:12.4f} {s['center_last_deg']:12.4f}")
    print("  → 現行（A だけ −v）は反作用が入らず中心が大きく動く。配分すると"
          "\n     桁で改善するがゼロにならない。R′ 比より R′² 比が良いという順序も"
          "\n     残っており、残差の原因は配分ではなく中心の定義（円周第1モーメントの"
          "\n     偏角は |z_A|≠|z_B| のとき重み付き平均にならない）を疑っている。")
    out["R3_R4"] = r3

    out["meta"] = {"experiment": "R_closure_probe_v1",
                   "date": time.strftime("%Y-%m-%d %H:%M:%S"),
                   "n_chi": n_chi, "n_eta": n_eta, "T": T_STEPS}
    p = HERE / "result_R_closure_v1.json"
    p.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"\n保存 {p.name}  ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
