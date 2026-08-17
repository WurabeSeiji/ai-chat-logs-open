#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""無名性——A・B の名前づけはゲージか。何が不変量か。

問い（木原）: AB が区別できると考えるのは無理ではないか。無名性の原則で、
明確に区別できる保存される毛が無い限り分からない。

ゲージ変換の定義
----------------
衝突 collision_step_exact は (a,b) の実直交回転でできている。したがって
**初期状態を任意の角 α で回すこと**が名前の付け替えに当たる:

    a → cosα·a − sinα·b ,  b → sinα·a + cosα·b        (G_α)

α=90° は A↔B の入れ替え。中間の α は「半分だけ A」という状態を作る。

G_α のもとでの各量（各格子点ごと、pointwise）
----------------------------------------------
  S0 = |a|²+|b|²      不変
  S2 = 2·Im(b̄a)       不変
  S3 = |a|²−|b|²  ┐
  S1 = 2·Re(b̄a)   ┘   二成分が **角 2α で回る**（不変ではない）
  a²+b²               不変（非共役の閉包）

  → 不変量は S0・S2・S1²+S3²・a²+b²。S3 と S1 は個別にはゲージ依存だが、
    Ψ = S3 + i·S1 と置くと Ψ → Ψ·e^{2iα} なので、**Ψ の位相差**は不変。

検定（実行前に固定）
--------------------
  G1 静止状態: α を振って各候補量を測り、α に依存しないものを同定する。
  G2 衝突のみの発展: 力学が G_α と可換なら、不変量は走行後も不変のはず。
  G3 CR4 本走行（並進あり）: 並進 k_translate_flat は **a にだけ**掛かる。
     これは名前づけに依存する操作なので、G_α と可換でないはず。
     可換でなければ **力学そのものが無名性を破っている**ことになる。

  判定は「α=0 との最大差」。1e-12 以下を不変とみなす。

使い方: python3 probe_gauge_invariants_v1.py
出力  : result_gauge_invariants_v1.json
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


_uni = _load("uni_g", UNI / "unified_interaction_v1.py")
K = _load("kin_g", UNI / "unified_kinetic_v1.py")
_cr0 = _load("cr0_g", HERE / "run_cr0_control_no_theta_v2.py")
_cr1 = _load("cr1_g", HERE / "run_cr1_kinetic_feedback_v1.py")
toy = _cr1.toy

DEG_A, DEG_B = -30.0, +30.0
ALPHAS = [0.0, 1.0, 5.0, 15.0, 45.0, 90.0]
OMEGA0 = np.pi / 72.0


def observables(a, b, n_chi, n_eta):
    """候補量。ゲージ依存のものも混ぜて並べる（依存することを示すため）。"""
    A = a.reshape(n_chi, n_eta)
    B = b.reshape(n_chi, n_eta)
    S0 = np.abs(A) ** 2 + np.abs(B) ** 2
    S3 = np.abs(A) ** 2 - np.abs(B) ** 2
    S1 = 2.0 * np.real(np.conj(B) * A)
    S2 = 2.0 * np.imag(np.conj(B) * A)

    P0 = S0.sum(axis=1); P0 = P0 / P0.sum()          # 対の位置分布（不変候補）
    w = np.exp(2j * np.pi * np.arange(n_chi) / n_chi)
    z0 = complex(P0 @ w)

    # Ψ = S3 + i·S1 は G_α で e^{2iα} 倍される。位相差だけが不変。
    Psi = (S3 + 1j * S1).sum(axis=1)
    Z0 = complex(Psi.sum())
    Z1 = complex(Psi @ w)

    pa, za = _cr0.circle_position(a, n_chi, n_eta)
    pb, zb = _cr0.circle_position(b, n_chi, n_eta)

    return {
        # --- 不変候補 ---
        "r": float(toy.theta_from_ab(a, b, sp_g).reflection_rate),
        "sumS0": float(S0.sum()),
        "sumS2": float(S2.sum()),
        "sumAbsS2": float(np.abs(S2).sum()),
        "sumHypot": float(np.hypot(S1, S3).sum()),
        "closureAbs": float(abs(complex(np.sum(a * a) + np.sum(b * b)))),
        "pair_pos_deg": float(np.degrees(np.angle(z0))),
        "pair_z": float(abs(z0)),
        "PsiPhaseDiff_deg": float(np.degrees(
            np.angle(Z1 * np.conj(Z0)))),          # Ψ の位相差（不変候補）
        "PsiRatio": float(abs(Z1) / max(abs(Z0), 1e-300)),
        # --- ゲージ依存とわかっているもの（対照） ---
        "dtheta_deg": float(np.degrees(np.angle(np.exp(1j * (pa - pb))))),
        "PA": float(np.vdot(a, a).real),
        "zA": za, "zB": zb,
        "PsiPhase_deg": float(np.degrees(np.angle(Z0))),
    }


def gauge(a, b, alpha_deg):
    c, s = np.cos(np.radians(alpha_deg)), np.sin(np.radians(alpha_deg))
    return c * a - s * b, s * a + c * b


def make(sp, slope, icept):
    case = _uni.two_body_base.explicit_packet_case(
        mode="gauge", packet_a=tuple(range(1, 18)), packet_b=(1, 2, 3),
        packet_a_shift=_cr0.shift_for_deg(DEG_A, slope, icept),
        packet_b_shift=_cr0.shift_for_deg(DEG_B, slope, icept))
    return (_uni.two_body_base.make_case_state(sp, case, "A", hair_enabled=True),
            _uni.two_body_base.make_case_state(sp, case, "B", hair_enabled=True))


def evolve(a, b, n_chi, n_eta, T, translate):
    omega, v = OMEGA0, 0.0
    for _ in range(T):
        if translate:
            pa, _ = _cr0.circle_position(a, n_chi, n_eta)
            pb, _ = _cr0.circle_position(b, n_chi, n_eta)
            chi = float(np.angle(np.exp(1j * (pa - pb))))
            r = float(toy.theta_from_ab(a, b, sp_g).reflection_rate)
            acc = -4.0 * np.sin(omega / 2.0) ** 2 * chi
            v += acc
            omega += (1.0 - r) * acc
            a = K.k_translate_flat(a, -v, n_chi, n_eta)
        a, b, _ = _uni.collision_step_exact(a, b, sp_g)
    return a, b


KEYS_INV = ["r", "sumS0", "sumS2", "sumAbsS2", "sumHypot", "closureAbs",
            "pair_pos_deg", "pair_z", "PsiPhaseDiff_deg", "PsiRatio"]
KEYS_DEP = ["dtheta_deg", "PA", "zA", "zB", "PsiPhase_deg"]


def table(title, T, translate, sp, n_chi, n_eta, slope, icept, out, tag):
    print(f"\n{'='*104}\n{title}\n{'='*104}")
    base = None
    rows = {}
    for ad in ALPHAS:
        a, b = make(sp, slope, icept)
        a, b = gauge(a, b, ad)
        if T:
            a, b = evolve(a, b, n_chi, n_eta, T, translate)
        o = observables(a, b, n_chi, n_eta)
        if base is None:
            base = o
        rows[ad] = {k: abs(o[k] - base[k]) for k in o}
    keys = KEYS_INV + KEYS_DEP
    print(f"  {'量':18} " + "".join(f"{('α='+str(int(x))+'°'):>12}" for x in ALPHAS[1:])
          + "   判定")
    print("  " + "-" * 100)
    verdict = {}
    for k in keys:
        d = [rows[ad][k] for ad in ALPHAS[1:]]
        inv = max(d) <= 1e-12
        verdict[k] = {"max_diff": float(max(d)), "invariant": bool(inv)}
        mark = "不変" if inv else "ゲージ依存"
        sep = "  ┄┄ 以下は対照（依存するはず）" if k == KEYS_DEP[0] else ""
        if sep:
            print("  " + "-" * 100)
        print(f"  {k:18} " + "".join(f"{x:12.3e}" for x in d) + f"   {mark}")
    out[tag] = verdict
    return verdict


def main() -> None:
    global sp_g
    t0 = time.time()
    sp_g = _uni.two_body_base.build_source_params(
        _uni.two_body_base.Params(high_n=63, recursive_collision_count=200))
    sp = sp_g
    n_chi, n_eta = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = _cr0.calibrate_shift(sp, n_chi, n_eta)
    out = {}

    v0 = table("【G1】静止状態（τ=0）——何が α に依らないか",
               0, False, sp, n_chi, n_eta, slope, icept, out, "G1_static")
    v1 = table("【G2】衝突のみ T=200（力学が G_α と可換なら不変量は保たれる）",
               200, False, sp, n_chi, n_eta, slope, icept, out, "G2_collision")
    v2 = table("【G3】CR4 本走行 T=200（並進あり——並進は a にだけ掛かる）",
               200, True, sp, n_chi, n_eta, slope, icept, out, "G3_translate")

    print(f"\n{'='*104}\n判定のまとめ\n{'='*104}")
    print(f"  {'量':18} {'静止':>10} {'衝突のみ':>10} {'並進あり':>10}")
    print("  " + "-" * 52)
    for k in KEYS_INV + KEYS_DEP:
        f = lambda v: "不変" if v[k]["invariant"] else "破れ"
        print(f"  {k:18} {f(v0):>10} {f(v1):>10} {f(v2):>10}")

    out["meta"] = {"experiment": "gauge_invariants_probe_v1",
                   "date": time.strftime("%Y-%m-%d %H:%M:%S"),
                   "alphas_deg": ALPHAS, "T": 200,
                   "threshold": 1e-12}
    p = HERE / "result_gauge_invariants_v1.json"
    p.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"\n保存 {p.name}  ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
