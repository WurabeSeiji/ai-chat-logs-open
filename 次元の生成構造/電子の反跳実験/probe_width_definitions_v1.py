#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""波束の幅——定義の候補を実測して比べる（決める前の掃引）

背景
----
CR4 まで局在度は参加率 PR = 1/Σ_χ P(χ)²（パワーの逆数二乗和）で測ってきた。
これは「実効セル数」であって波束の幅ではない。次に主ローブの零点（振幅の
符号反転）で測ったが、これは **きれいな主ローブがある** ことと **中心が
決まる** ことの二つを前提に持ち込む。衝突後の波形にどちらの保証もない。

したがって幅の定義は
  (a) 中心を要求しない
  (b) 単峰（主ローブ1本）を仮定しない
の二条件を満たす必要がある。本スクリプトは候補を並べて実測で比べる。
**この段階では定義を決めない。決めるための材料を出すだけ。**

測る対象は χ 円周上の **振幅** A(χ) = sqrt(Σ_η |ψ(χ,η)|²)。
パワーではなく振幅で測る（PR がパワー由来である点が元の誤り）。

候補
----
  W(q)  レベル集合幅: A ≥ q·A_max となるセルの総角度。
        中心を要求せず、多峰でも定義できる（総支持長を測る）。
        木原案は q = 0.10。
  NC(q) 同レベル集合の連結成分数。1 なら単峰、2 以上なら主ローブ仮定が崩れている。
        —— W(q) の妥当性そのものを監視する計器。
  ERB   等価矩形幅 (Σ A)² / Σ A² × Δχ。閾値を持たない。振幅版の参加率。
  PRp   参加率（現行）1/Σ P² × Δχ。パワー由来。比較のため併記。
  CSD   円周標準偏差 sqrt(−2 ln|z|)。第1モーメントのみ。単峰を暗に仮定。

判定の見方
----------
  1. 静止パケット（倍音 1..M）では理論の主ローブ半幅 360/(2M+1) が既知。
     全幅は 2×360/(2M+1) = 720/(2M+1)。各候補がこの 1/(2M+1) 則に
     どれだけ乗るか、比例係数が M によらず一定かを見る。
  2. **閾値感度**: q を 0.02〜0.50 で掃引し、W(q) が平坦域を持つか見る。
     平坦なら閾値の選び方が結果を左右しない＝定義として使える。
     崖の途中なら q=0.10 は恣意的な選択になる。
  3. **衝突後の頑健性**: 実験条件で走らせ、各候補の時間発展を見る。
     連結成分数が 1 を離れる時刻で「主ローブ」概念が壊れる。

使い方: python3 probe_width_definitions_v1.py
出力  : result_width_definitions_v1.json
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


_uni = _load("uni_w", UNI / "unified_interaction_v1.py")
K = _load("kin_w", UNI / "unified_kinetic_v1.py")
_cr0 = _load("cr0_w", HERE / "run_cr0_control_no_theta_v2.py")
_cr1 = _load("cr1_w", HERE / "run_cr1_kinetic_feedback_v1.py")

DEG_A, DEG_B = -30.0, +30.0
OMEGA0 = np.pi / 72.0
M_LIST = [1, 2, 3, 5, 9, 17, 25]
Q_SWEEP = [0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]


# ------------------------------------------------------------- 計器（候補）

def amp_chi(psi, n_chi, n_eta):
    """χ 円周上の振幅 A(χ) = sqrt(Σ_η |ψ|²)。パワーではなく振幅。"""
    P = np.sum(np.abs(psi.reshape(n_chi, n_eta)) ** 2, axis=1)
    return np.sqrt(P)


def level_set_width(A, q):
    """レベル集合幅[°]と連結成分数。円周なので端をつなぐ。

    中心を要求しない。単峰も仮定しない。A ≥ q·A_max の総支持長を返す。
    """
    n = len(A)
    mask = A >= q * A.max()
    width = float(mask.sum()) * 360.0 / n
    # 円周上の連結成分数（全 True なら 1、全 False は起こらない）
    if mask.all():
        return width, 1
    d = np.diff(np.concatenate([mask, mask[:1]]).astype(int))
    return width, int(np.sum(d == 1))


def erb(A, n_chi):
    """等価矩形幅[°]。閾値なし。(Σ A)²/Σ A² × Δχ。"""
    return float((A.sum() ** 2) / np.sum(A ** 2)) * 360.0 / n_chi


def pr_power(psi, n_chi, n_eta):
    """参加率（現行・パワー由来）を角度[°]に直したもの。"""
    return _cr0.participation_ratio(psi, n_chi, n_eta) * 360.0 / n_chi


def circ_sd(psi, n_chi, n_eta):
    """円周標準偏差[°] = sqrt(−2 ln|z|)。"""
    _, z = _cr0.circle_position(psi, n_chi, n_eta)
    z = min(max(z, 1e-12), 1.0 - 1e-15)
    return float(np.degrees(np.sqrt(-2.0 * np.log(z))))


def measure_all(psi, n_chi, n_eta):
    A = amp_chi(psi, n_chi, n_eta)
    out = {"erb": erb(A, n_chi), "pr_power": pr_power(psi, n_chi, n_eta),
           "csd": circ_sd(psi, n_chi, n_eta)}
    for q in Q_SWEEP:
        w, nc = level_set_width(A, q)
        out[f"W{int(q*100):02d}"] = w
        out[f"NC{int(q*100):02d}"] = nc
    return out


# ------------------------------------------------------------- 走行

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
    print(f"格子 n_chi={n_chi} n_eta={n_eta}  Δχ={360.0/n_chi:.4f}°  "
          f"較正残差={resid:.2e}")

    # --- 1. 静止パケットの倍音掃引 -----------------------------------
    print("\n【1】静止パケット: 倍音 1..M の幅と 1/(2M+1) 則")
    static = {}
    for M in M_LIST:
        psi = build_state(sp, range(1, M + 1), DEG_A, slope, icept, f"probe_M{M}")
        m = measure_all(psi, n_chi, n_eta)
        m["theory_full"] = 720.0 / (2 * M + 1)
        static[M] = m

    hdr = ("  M   理論全幅   W10    NC10    W50    NC50    ERB    PRp     CSD"
           "   W10/理論")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for M in M_LIST:
        m = static[M]
        print(f"{M:3d} {m['theory_full']:9.3f} {m['W10']:7.3f} {m['NC10']:5d} "
              f"{m['W50']:7.3f} {m['NC50']:5d} {m['erb']:7.2f} "
              f"{m['pr_power']:7.2f} {m['csd']:7.2f} "
              f"{m['W10']/m['theory_full']:8.4f}")

    # --- 2. 閾値感度 --------------------------------------------------
    print("\n【2】閾値感度: W(q)[°]（上）と連結成分数（下）")
    print("   q  " + "".join(f"{M:>9d}" for M in M_LIST))
    for q in Q_SWEEP:
        k = f"W{int(q*100):02d}"
        kc = f"NC{int(q*100):02d}"
        print(f"{q:5.2f} " + "".join(f"{static[M][k]:9.3f}" for M in M_LIST)
              + "   | " + "".join(f"{static[M][kc]:2d}" for M in M_LIST))

    # 平坦性: 隣接する q の間で W が何%動くか（M ごとの最大変化率）
    print("\n   平坦性（隣接 q 間の |ΔW|/W の最大値[%]、小さいほど閾値に鈍感）")
    flat = {}
    for M in M_LIST:
        ws = [static[M][f"W{int(q*100):02d}"] for q in Q_SWEEP]
        rel = [abs(ws[i + 1] - ws[i]) / ws[i] * 100 for i in range(len(ws) - 1)]
        flat[M] = {"per_step": [round(x, 2) for x in rel], "max": max(rel)}
        print(f"   M={M:2d}: 最大 {max(rel):6.2f}%   "
              + " ".join(f"{x:5.1f}" for x in rel))

    # --- 3. 衝突下の頑健性 --------------------------------------------
    print("\n【3】衝突下の時間発展（実験条件 A:1〜17 / B:1〜3、T=600）")
    case = _uni.two_body_base.explicit_packet_case(
        mode="probe_width_dyn", packet_a=tuple(range(1, 18)),
        packet_b=tuple(range(1, 4)),
        packet_a_shift=_cr0.shift_for_deg(DEG_A, slope, icept),
        packet_b_shift=_cr0.shift_for_deg(DEG_B, slope, icept))
    a = _uni.two_body_base.make_case_state(sp, case, "A", hair_enabled=True)
    b = _uni.two_body_base.make_case_state(sp, case, "B", hair_enabled=True)

    T = 600
    omega, v = OMEGA0, 0.0
    dyn = {k: [] for k in ("t", "A_W10", "A_NC10", "A_erb", "A_pr", "A_csd",
                           "B_W10", "B_NC10", "B_erb", "B_pr", "B_csd")}
    for t in range(T):
        pa, _ = _cr0.circle_position(a, n_chi, n_eta)
        pb, _ = _cr0.circle_position(b, n_chi, n_eta)
        chi = float(np.angle(np.exp(1j * (pa - pb))))
        r_now = float(_cr1.toy.theta_from_ab(a, b, sp).reflection_rate)
        acc = -4.0 * np.sin(omega / 2.0) ** 2 * chi
        v += acc
        omega += (1.0 - r_now) * acc          # κ = 1−r（CR1 で決めた形）

        for nm, psi in (("A", a), ("B", b)):
            Aa = amp_chi(psi, n_chi, n_eta)
            w, nc = level_set_width(Aa, 0.10)
            dyn[f"{nm}_W10"].append(round(w, 4))
            dyn[f"{nm}_NC10"].append(nc)
            dyn[f"{nm}_erb"].append(round(erb(Aa, n_chi), 4))
            dyn[f"{nm}_pr"].append(round(pr_power(psi, n_chi, n_eta), 4))
            dyn[f"{nm}_csd"].append(round(circ_sd(psi, n_chi, n_eta), 4))
        dyn["t"].append(t)

        a = K.k_translate_flat(a, -v, n_chi, n_eta)
        a, b, _ = _uni.collision_step_exact(a, b, sp)

    def stat(key):
        x = np.array(dyn[key], float)
        return x[0], x.mean(), x.std(), x.min(), x.max()

    print(f"  {'計器':12} {'初期':>9} {'平均':>9} {'標準偏差':>9} "
          f"{'最小':>9} {'最大':>9}  変動[%]")
    for nm in ("A", "B"):
        for key, lab in (("W10", "W10"), ("erb", "ERB"), ("pr", "PR(power)"),
                         ("csd", "CSD")):
            k = f"{nm}_{key}"
            i0, mu, sd, lo, hi = stat(k)
            print(f"  {nm}:{lab:10} {i0:9.3f} {mu:9.3f} {sd:9.3f} "
                  f"{lo:9.3f} {hi:9.3f}  {100*sd/mu:7.3f}")
    for nm in ("A", "B"):
        nc = np.array(dyn[f"{nm}_NC10"], int)
        first = int(np.argmax(nc > 1)) if (nc > 1).any() else -1
        print(f"  {nm}: 連結成分数 NC10  最大={nc.max()}  "
              f"多峰になった步数={int((nc>1).sum())}/{T}  "
              f"初回={'なし' if first<0 else first}")

    out = {
        "experiment": "width_definitions_probe_v1",
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "grid": {"n_chi": n_chi, "n_eta": n_eta, "d_chi_deg": 360.0 / n_chi,
                 "calib_resid": resid},
        "q_sweep": Q_SWEEP,
        "static": {str(M): static[M] for M in M_LIST},
        "flatness": {str(M): flat[M] for M in M_LIST},
        "dynamic": dyn,
        "note": "定義は未確定。候補の比較材料のみ。",
    }
    p = HERE / "result_width_definitions_v1.json"
    p.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"\n保存 {p.name}  ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
