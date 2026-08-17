#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""波束の幅 v2——閾値ではなく「裾の重み」を掃引する

v1 の結論
---------
  * 静止パケットの時点で 10% レベル集合の連結成分は M=1 で 2、M=17 で 7。
    **主ローブは1本ではない。** 主ローブ前提の定義は初期条件から成立しない。
  * W(q) は閾値に平坦域を持たない（隣接 q 間で最大 69% 動く）。
    q=0.10 は崖の途中の一点で、恣意的な選択になる。

したがって閾値を捨てる。等振幅倍音 1..M の振幅は Dirichlet 核であり、
副ローブが 1/χ でしか落ちない。幅が閾値でなく **裾にどれだけ重みを置くか**
で決まるのはこの裾のせいで、これは系の性質であって計器の欠陥ではない。

そこで「裾の重み」を一つの指数 s に集約して掃引する。

  密度  d_i ∝ A_i^s        （A は χ 円周上の振幅）
  幅    W_s = Δχ / Σ d_i²  （逆参加数を角度に直したもの）

  s=1 は振幅の等価矩形幅（v1 の ERB）、s=2 は現行の参加率 PR（パワー版）に
  一致する。s を上げるほど裾が軽くなり、幅は中心の峰に寄る。
  閾値と違い s は連続で、どの s でも定義は well-defined。

決め方（恣意性を測定に置き換える）
-----------------------------------
  基準1 スケール則: 幅は倍音本数 M に対し 1/(2M+1) に比例すべき。
        W_s·(2M+1) が M によらず一定になる s を探す。ばらつきを CV で測る。
  基準2 衝突下の頑健性: 変動係数 std/mean が小さいこと。
  基準3 分別能: 衝突下で A と B を区別し続けられること（初期の非対称を保つ）。

使い方: python3 probe_width_definitions_v2.py
出力  : result_width_definitions_v2.json
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


_uni = _load("uni_w2", UNI / "unified_interaction_v1.py")
K = _load("kin_w2", UNI / "unified_kinetic_v1.py")
_cr0 = _load("cr0_w2", HERE / "run_cr0_control_no_theta_v2.py")
_cr1 = _load("cr1_w2", HERE / "run_cr1_kinetic_feedback_v1.py")

DEG_A, DEG_B = -30.0, +30.0
OMEGA0 = np.pi / 72.0
M_LIST = [1, 2, 3, 5, 7, 9, 13, 17, 21, 25]
S_LIST = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 8.0]


def amp_chi(psi, n_chi, n_eta):
    P = np.sum(np.abs(psi.reshape(n_chi, n_eta)) ** 2, axis=1)
    return np.sqrt(P)


def width_s(A, s, n_chi):
    """裾の重み s の幅[°]。s=1 は振幅 ERB、s=2 は現行 PR に一致する。"""
    d = A ** s
    d = d / d.sum()
    return float(1.0 / np.sum(d ** 2)) * 360.0 / n_chi


def build_state(sp, harmonics, deg, slope, icept, tag):
    case = _uni.two_body_base.explicit_packet_case(
        mode=tag, packet_a=tuple(harmonics), packet_b=tuple(harmonics),
        packet_a_shift=_cr0.shift_for_deg(deg, slope, icept),
        packet_b_shift=_cr0.shift_for_deg(deg, slope, icept))
    return _uni.two_body_base.make_case_state(sp, case, "A", hair_enabled=True)


def main() -> None:
    t0 = time.time()
    sp = _uni.two_body_base.build_source_params(
        _uni.two_body_base.Params(high_n=63, recursive_collision_count=200))
    n_chi, n_eta = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, resid = _cr0.calibrate_shift(sp, n_chi, n_eta)
    print(f"格子 n_chi={n_chi} n_eta={n_eta}  Δχ={360.0/n_chi:.4f}°")

    # --- 基準1: スケール則 --------------------------------------------
    print("\n【1】静止パケット W_s[°]（理論 720/(2M+1) と比較）")
    amps = {}
    for M in M_LIST:
        psi = build_state(sp, range(1, M + 1), DEG_A, slope, icept, f"p2_M{M}")
        amps[M] = amp_chi(psi, n_chi, n_eta)

    W = {s: {M: width_s(amps[M], s, n_chi) for M in M_LIST} for s in S_LIST}
    print("   s  " + "".join(f"{M:>8d}" for M in M_LIST))
    print("  理論 " + "".join(f"{720.0/(2*M+1):8.2f}" for M in M_LIST))
    print("  " + "-" * (5 + 8 * len(M_LIST)))
    for s in S_LIST:
        print(f"{s:5.1f} " + "".join(f"{W[s][M]:8.2f}" for M in M_LIST))

    print("\n【2】スケール則の一定性  W_s·(2M+1)（一定なら 1/(2M+1) 則に乗る）")
    print("   s  " + "".join(f"{M:>8d}" for M in M_LIST) + "     CV[%]")
    best_s, best_cv = None, 1e9
    cvs = {}
    for s in S_LIST:
        prod = np.array([W[s][M] * (2 * M + 1) for M in M_LIST])
        cv = 100.0 * prod.std() / prod.mean()
        cvs[s] = cv
        if cv < best_cv:
            best_cv, best_s = cv, s
        print(f"{s:5.1f} " + "".join(f"{p:8.1f}" for p in prod) + f"  {cv:8.2f}")
    print(f"  → スケール則に最も乗る s = {best_s}（CV {best_cv:.2f}%）")

    # 局所指数 dlnW/dln(2M+1) も見る（1/(2M+1) 則なら −1）
    print("\n【3】局所スケール指数 dlnW/dln(2M+1)（1/(2M+1) 則なら −1）")
    x = np.log(np.array([2 * M + 1 for M in M_LIST], float))
    print("   s   " + "".join(f"{M_LIST[i]}→{M_LIST[i+1]:<5d}"
                              for i in range(len(M_LIST) - 1)))
    for s in S_LIST:
        y = np.log(np.array([W[s][M] for M in M_LIST]))
        g = np.diff(y) / np.diff(x)
        print(f"{s:5.1f}  " + "".join(f"{v:8.3f}" for v in g))

    # --- 基準2・3: 衝突下 ---------------------------------------------
    print("\n【4】衝突下の時間発展（A:1〜17 / B:1〜3、T=600）")
    case = _uni.two_body_base.explicit_packet_case(
        mode="probe_width2_dyn", packet_a=tuple(range(1, 18)),
        packet_b=tuple(range(1, 4)),
        packet_a_shift=_cr0.shift_for_deg(DEG_A, slope, icept),
        packet_b_shift=_cr0.shift_for_deg(DEG_B, slope, icept))
    a = _uni.two_body_base.make_case_state(sp, case, "A", hair_enabled=True)
    b = _uni.two_body_base.make_case_state(sp, case, "B", hair_enabled=True)

    T = 600
    omega, v = OMEGA0, 0.0
    hist = {f"{nm}_{s}": [] for nm in ("A", "B") for s in S_LIST}
    for _ in range(T):
        pa, _ = _cr0.circle_position(a, n_chi, n_eta)
        pb, _ = _cr0.circle_position(b, n_chi, n_eta)
        chi = float(np.angle(np.exp(1j * (pa - pb))))
        r_now = float(_cr1.toy.theta_from_ab(a, b, sp).reflection_rate)
        acc = -4.0 * np.sin(omega / 2.0) ** 2 * chi
        v += acc
        omega += (1.0 - r_now) * acc

        for nm, psi in (("A", a), ("B", b)):
            Aa = amp_chi(psi, n_chi, n_eta)
            for s in S_LIST:
                hist[f"{nm}_{s}"].append(width_s(Aa, s, n_chi))

        a = K.k_translate_flat(a, -v, n_chi, n_eta)
        a, b, _ = _uni.collision_step_exact(a, b, sp)

    print(f"   s  {'A初期':>8} {'A平均':>8} {'A変動%':>8} "
          f"{'B初期':>8} {'B平均':>8} {'B変動%':>8}  {'分別能':>8}")
    sep = {}
    for s in S_LIST:
        A_ = np.array(hist[f"A_{s}"]); B_ = np.array(hist[f"B_{s}"])
        # 分別能: 走行中の |A−B| の平均を、初期の差で割ったもの（1 なら保持）
        d0 = abs(A_[0] - B_[0])
        keep = float(np.mean(np.abs(A_ - B_)) / d0) if d0 > 0 else float("nan")
        sep[s] = keep
        print(f"{s:5.1f} {A_[0]:8.2f} {A_.mean():8.2f} "
              f"{100*A_.std()/A_.mean():8.2f} {B_[0]:8.2f} {B_.mean():8.2f} "
              f"{100*B_.std()/B_.mean():8.2f}  {keep:8.3f}")

    out = {
        "experiment": "width_definitions_probe_v2",
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "grid": {"n_chi": n_chi, "n_eta": n_eta, "calib_resid": resid},
        "s_list": S_LIST, "m_list": M_LIST,
        "static_W": {str(s): {str(M): W[s][M] for M in M_LIST} for s in S_LIST},
        "scaling_cv": {str(s): cvs[s] for s in S_LIST},
        "best_s_by_scaling": best_s,
        "dyn": {k: [round(x, 4) for x in v] for k, v in hist.items()},
        "separability": {str(s): sep[s] for s in S_LIST},
        "note": "定義は未確定。s の決め方の材料。",
    }
    p = HERE / "result_width_definitions_v2.json"
    p.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"\n保存 {p.name}  ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
