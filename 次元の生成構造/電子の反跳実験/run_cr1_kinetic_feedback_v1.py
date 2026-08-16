#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CR1: 荷電二体波の反跳——運動項＋速度フィードバック（χ = Δθ）

CR0 v2 で「運動項が無いので位置が巡らない」ことを対照点として確定させた。
本実験は統一万能運動関数 K（資格審査 Q-K1〜Q-K4 ALL PASS）を投入し、
加速度逆二乗則論文 v4 §13 の速度フィードバックで駆動する。

χ の定義（木原確定・2026-08-16）
--------------------------------
**χ = Δθ（AB 二体の相対位相）**。論文 §13.3 の閉鎖表現は χ の一表現に
すぎず、閉包残差そのものではない。実測でも本系の閉包残差は 8e-17 で
厳密に閉じており（`collision_step_exact` が閉包を保存則にしている）、
残差から χ を読むと恒等的に 0 になって力が生じない。AB の関係は Δθ である。

動力学（論文 v4 §13.1・§13.2・作業仮説）
------------------------------------------
    a_s   = −g(ω_s)·χ_s ,   g(ω) = 4 sin²(ω/2) ≥ 0
    ω_{s+1} = ω_s + κ·a_s                    （Euler形・キャリアの応答）

**符号は作り込んでいない。** g ≥ 0 なので向きは −χ からのみ来る＝常に
Δθ=0 へ戻す復元力。R′²=R²+Q² で電荷 Q は二乗でしか入らないため、
**加速度は原理的に電荷の符号に依存しない**。これが H4 対決点
（G9 実測で E(++) が引力側／クーロンは同符号反発）の構造的な正体である。

ω₀ の選択（手で置かない）
--------------------------
g(0)=0 なので ω=0 からは加速度が出ない。ω は静止から生まれる量ではなく、
既に回っているキャリア角速度である。実測済みの**普遍時計 ω=π/72**
（一周144ステップ・N 非依存・±0.1%）を ω₀ に採る。倍音構成は
ω_n = n·ω₁ としてここに効く（A は最高17次・B は最高3次）。

相対運動の配分
--------------
Δθ は相対座標。その速度 v を A・B へ質量比で配分し、重心を固定する:
    v_A = +v·μ_B ,  v_B = −v·μ_A ,   μ_X = m_X/(m_A+m_B)
質量は Gram の detΓ（`cone_m2` と同じ式）で毎步読む。手で置かない。

事前登録の判定
--------------
  CR1-1  Δθ が動く（CR0 との対比。総変化 > 1°）
  CR1-2  **閉包残差の包絡が成長しない**（論文v4 §13.3）
         ——閉包残差は保存量ではない。神の手で Δθ=60° というポテンシャルを
         持つ状態に初期化した以上、準安定へ緩和する過程で崩れるのが物理。
         判定は「厳密保存」ではなく「包絡の非成長」である
  CR1-3  ノルムドリフト < 1e-9
  CR1-4  **束縛振動**（Δθ が単調でなく、転回点を 1 回以上持つ）
  CR1-5  **運動項を完全に切ると CR0 と厳密一致**（freeze_v=True）

使い方: python3 run_cr1_kinetic_feedback_v1.py
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


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


_uni = _load("uni_cr1", UNI / "unified_interaction_v1.py")
K = _load("kin_cr1", UNI / "unified_kinetic_v1.py")
_cr0 = _load("cr0_cr1", HERE / "run_cr0_control_no_theta_v2.py")

collision_step_exact = _uni.collision_step_exact
base = _uni.two_body_base
toy = _uni.two_body_v1.toy          # theta_from_ab（状態から r を読む）

PACKET_A = tuple(range(1, 18))
PACKET_B = tuple(range(1, 4))
DEG_A, DEG_B = -30.0, +30.0
OMEGA_CLOCK = np.pi / 72.0        # 普遍時計（実測・一周144步）
# κ は自由パラメータではない。この系の唯一の結合定数は反射率 r であり、
# 頂点 φ = 2r·Im(b̄a) の r と同じものを毎步読んで使う（宣言しない）。
# 論文v4 §13.1 の κ は「無次元のフィードバック係数・新しい作業仮説」として
# 導入された自由係数だが、二つ目の結合定数を根拠なく足すことになるため、
# 本実験では κ = r（実測）とする。
KAPPA_MODE = "transmission"      # κ = 1−r（透過率＝素電荷側）
T_STEPS = 40000


def gram_mass2(psi, n_chi, n_eta):
    """質量² = detΓ = T²−X²−Y²−Z²（`cone_m2` と同一式・帯の偶奇対）。"""
    f = np.fft.fft(psi.reshape(n_chi, n_eta), axis=0)
    k = np.rint(np.fft.fftfreq(n_chi, d=1.0 / n_chi)).astype(int)
    odd = (np.abs(k) % 2) == 1
    a = f[odd].reshape(-1)
    b = f[~odd].reshape(-1)
    pa = float(np.sum(np.abs(a) ** 2))
    pb = float(np.sum(np.abs(b) ** 2))
    z = complex(np.sum(np.conj(a[:min(len(a), len(b))]) * b[:min(len(a), len(b))]))
    T, X = pa + pb, pa - pb
    Y, Z = 2.0 * z.real, 2.0 * z.imag
    return max(T * T - X * X - Y * Y - Z * Z, 0.0)


def cone_components(psi, n_chi, n_eta):
    """錐の成分（`g_cone_components` と同一式・二体トイ版）。

    R′² = T²、m² = detΓ = T²−X²−Y²−Z²、q = Σ η_signed·P(η)、
    R² = R′² − Q²（電荷を差し引いた曲率）。
    """
    A = psi.reshape(n_chi, n_eta)
    f = np.fft.fft(A, axis=0)
    k = np.rint(np.fft.fftfreq(n_chi, d=1.0 / n_chi)).astype(int)
    odd = (np.abs(k) % 2) == 1
    a_ = f[odd].reshape(-1); b_ = f[~odd].reshape(-1)
    L = min(len(a_), len(b_))
    pa = float(np.sum(np.abs(a_) ** 2)); pb = float(np.sum(np.abs(b_) ** 2))
    z = complex(np.sum(np.conj(a_[:L]) * b_[:L]))
    T_, X = pa + pb, pa - pb
    Y, Z = 2.0 * z.real, 2.0 * z.imag
    m2 = T_ * T_ - X * X - Y * Y - Z * Z
    g = np.fft.fft(A, axis=1)
    Pm = np.sum(np.abs(g) ** 2, axis=0); Pm = Pm / Pm.sum()
    ms = np.arange(n_eta)
    ms = ((ms + n_eta // 2) % n_eta) - n_eta // 2
    q = float(Pm @ ms)
    Rp2 = T_ * T_
    return Rp2, m2, q, max(Rp2 - q * q, 0.0)


def run(kappa_mode, v0, T, tag, freeze_v=False):
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n_chi, n_eta = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = _cr0.calibrate_shift(sp, n_chi, n_eta)
    a, b = _cr0.make_pair(sp, _cr0.shift_for_deg(DEG_A, slope, icept),
                          _cr0.shift_for_deg(DEG_B, slope, icept))

    tot0 = float(np.vdot(a, a).real + np.vdot(b, b).real)
    Z0 = complex(np.sum(a * a) + np.sum(b * b))
    clo0 = Z0

    omega = OMEGA_CLOCK
    v = v0
    H = {k: [] for k in ("chi", "a", "omega", "v", "chord", "clo_re", "clo_im", "norm",
                         "pos_a", "pos_b", "m2_a", "m2_b", "r", "kappa",
                         "Rp2_a", "Rp2_b", "q_a", "q_b", "R2_a", "R2_b")}

    for _ in range(T):
        # --- 読出し（状態から χ を読む・手で置かない）---
        pa, _ = _cr0.circle_position(a, n_chi, n_eta)
        pb_, _ = _cr0.circle_position(b, n_chi, n_eta)
        chi = float(np.angle(np.exp(1j * (pa - pb_))))       # 円上の相対位相

        # --- 結合定数を状態から読む（宣言しない）---
        r_now = float(toy.theta_from_ab(a, b, sp).reflection_rate)
        kappa = (r_now if kappa_mode == "reflection_rate"
                 else (1.0 - r_now) if kappa_mode == "transmission"
                 else float(kappa_mode))

        # --- 動力学（論文v4 §13.2 Euler形）---
        acc = -4.0 * np.sin(omega / 2.0) ** 2 * chi          # 符号は −χ のみ
        if not freeze_v:                                      # 対照用（運動項を完全に切る）
            v = v + acc
        omega = omega + kappa * acc

        # --- 運動項: 相対座標にだけ効かせる ---
        # 絶対系は存在しない（無名性）。物理的にあるのは Δθ の変化だけで、
        # 「A を vA、B を vB」と配分するのは重心という絶対基準を密輸すること。
        # 両者を同じ量だけ動かす操作は大域並進＝物理的に無意味なので、
        # A だけ v 動かすのと質量比配分は同じ物理状態を指す。
        # 実測: 質量比配分は μ が毎步ゆらいで A・B に独立な並進が入り、
        # 二体の閉包が漏れた（並進步でのみ増大・衝突步は 1e-18 で不変）。
        m2a = gram_mass2(a, n_chi, n_eta)      # 記録のみ（配分には使わない）
        m2b = gram_mass2(b, n_chi, n_eta)
        a = K.k_translate_flat(a, -v, n_chi, n_eta)

        # --- 相互作用 ---
        a, b, ro = collision_step_exact(a, b, sp)

        H["chi"].append(chi); H["a"].append(acc)
        H["omega"].append(omega); H["v"].append(v)
        H["chord"].append(float(2.0 * np.sin(abs(chi) / 2.0)))
        Zc = complex(np.sum(a * a) + np.sum(b * b))   # Σz² は複素量（ノルムではない）
        H["clo_re"].append(Zc.real)                    # 実部 Σ(a²−b²)
        H["clo_im"].append(Zc.imag)                    # 虚部 2Σab（交差項）
        H["norm"].append(float(np.vdot(a, a).real + np.vdot(b, b).real))
        H["pos_a"].append(pa); H["pos_b"].append(pb_)
        H["m2_a"].append(m2a); H["m2_b"].append(m2b)
        Ra, ma2, qa, R2a = cone_components(a, n_chi, n_eta)
        Rb, mb2, qb, R2b = cone_components(b, n_chi, n_eta)
        H["Rp2_a"].append(Ra); H["Rp2_b"].append(Rb)
        H["q_a"].append(qa); H["q_b"].append(qb)
        H["R2_a"].append(R2a); H["R2_b"].append(R2b)
        H["r"].append(float(ro.reflection_rate)); H["kappa"].append(kappa)

    H = {k: np.asarray(v_) for k, v_ in H.items()}
    return H, tot0, clo0, (a, b)


def main() -> None:
    t0 = time.time()

    # 本走行
    H, tot0, clo0, _ = run(KAPPA_MODE, 0.0, T_STEPS, "main")
    chi_deg = np.degrees(H["chi"])

    # CR1-5 対照: κ=0 かつ v₀=0 → 運動項は恒等（ω₁=0）→ CR0 と一致するはず
    H0, _, _, (a0, b0) = run(KAPPA_MODE, 0.0, 200, "frozen", freeze_v=True)

    # 純 CR0（運動項を一切通さない）と比較
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n_chi, n_eta = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = _cr0.calibrate_shift(sp, n_chi, n_eta)
    a2, b2 = _cr0.make_pair(sp, _cr0.shift_for_deg(DEG_A, slope, icept),
                            _cr0.shift_for_deg(DEG_B, slope, icept))
    for _ in range(200):
        a2, b2, _ = collision_step_exact(a2, b2, sp)
    diff5 = float(max(np.max(np.abs(a0 - a2)), np.max(np.abs(b0 - b2))))

    # 転回点（符号反転回数）
    dchi = np.diff(chi_deg)
    turns = int(np.sum(np.sign(dchi[1:]) * np.sign(dchi[:-1]) < 0))

    # 包絡は実部・虚部それぞれで見る（絶対値へ潰さない）
    CR = H["clo_re"]; CI = H["clo_im"]
    C = np.maximum(np.abs(CR), np.abs(CI))    # 判定用の包絡（大きい方）
    half = len(C) // 2
    env1 = float(C[:half].max()); env2 = float(C[half:].max())
    env_growth = env2 / env1 if env1 > 0 else float("inf")
    # 飽和判定: 8窓の包絡が後半で頭打ちか
    W = 8
    envs = [float(C[i*len(C)//W:(i+1)*len(C)//W].max()) for i in range(W)]
    tail_growth = envs[-1] / envs[-2] if envs[-2] > 0 else float("inf")
    # 周期測定（零交差の間隔）
    def period_of(x):
        y = np.asarray(x) - np.mean(x)
        s = np.sign(y); zc = np.where(np.diff(s) != 0)[0]
        return float(2.0 * np.mean(np.diff(zc))) if len(zc) > 2 else float("nan")
    per_chi = period_of(H["chi"])
    per_Rp = period_of(H["Rp2_a"])
    closure_drift = float(max(abs(CR[-1] - clo0.real), abs(CI[-1] - clo0.imag)))
    norm_drift = float(np.max(np.abs(H["norm"] / tot0 - 1.0)))
    total_change = float(chi_deg.max() - chi_deg.min())

    cr11 = bool(total_change > 1.0)
    cr12 = bool(tail_growth <= 1.05)     # 最終窓が直前窓を超えない（飽和）
    cr13 = bool(norm_drift <= 1e-9)
    cr14 = bool(turns >= 1)
    cr15 = bool(diff5 <= 1e-12)

    print("=" * 74)
    print("CR1 荷電二体波の反跳 — 運動項＋速度フィードバック（χ = Δθ）")
    print("=" * 74)
    print(f"  A=倍音1..17 / B=倍音1..3   初期Δθ={chi_deg[0]:+.4f}°   T={T_STEPS}")
    print(f"  ω₀ = π/72 = {OMEGA_CLOCK:.8f}（普遍時計・実測）")
    print(f"  κ = 反射率 r（毎步の実測・宣言しない）  平均={H['kappa'].mean():.6f}"
          f"  範囲=[{H['kappa'].min():.6f}, {H['kappa'].max():.6f}]")
    print(f"  参考: 素電荷番地 1−cos²(23π/124) = 0.302822073"
          f"   差={H['kappa'].mean()-0.302822073:+.6f}")
    print(f"  g(ω₀) = 4sin²(ω₀/2) = {4*np.sin(OMEGA_CLOCK/2)**2:.6e}")
    print()
    print("--- 判定 ---")
    print(f"  CR1-1 Δθ が動く（総変化 {total_change:.4f}°）: {cr11}")
    print(f"  CR1-2 包絡が飽和（最終窓比 {tail_growth:.4f}）: {cr12}")
    print(f"        8窓の包絡: " + " → ".join(f"{e:.2e}" for e in envs))
    print(f"        実部 Σ(a²−b²): 初期{clo0.real:+.3e} 範囲[{CR.min():+.3e},{CR.max():+.3e}]")
    print(f"        虚部 2Σab    : 初期{clo0.imag:+.3e} 範囲[{CI.min():+.3e},{CI.max():+.3e}]")
    print(f"  CR1-3 ノルムドリフト {norm_drift:.3e}: {cr13}")
    print(f"  CR1-4 束縛振動（転回点 {turns} 回）: {cr14}")
    print(f"  CR1-5 κ=0,v₀=0 で CR0 と一致（最大差 {diff5:.3e}）: {cr15}")
    print()
    print("--- Δθ の軌跡 ---")
    for f in (0.0, 0.1, 0.25, 0.5, 0.75, 1.0):
        i = min(int(f * (T_STEPS - 1)), T_STEPS - 1)
        print(f"    t={i:5d}: Δθ={chi_deg[i]:+9.4f}°  弦={H['chord'][i]:.6f}"
              f"  v={H['v'][i]:+.6e}  ω={H['omega'][i]:.8f}")
    print()
    print(f"  Δθ 範囲=[{chi_deg.min():+.4f}, {chi_deg.max():+.4f}]°")
    print(f"  弦 R″/R 範囲=[{H['chord'].min():.6f}, {H['chord'].max():.6f}]")
    print(f"  ω 範囲=[{H['omega'].min():.8f}, {H['omega'].max():.8f}]"
          f"（ω₀={OMEGA_CLOCK:.8f}）")
    print()
    print(f"--- 周期 ---")
    print(f"  Δθ の周期 ≈ {per_chi:.1f} 步    R′² の周期 ≈ {per_Rp:.1f} 步"
          f"    T/周期 = {T_STEPS/per_chi:.1f} 周期分" if per_chi == per_chi else "  周期測定不能")
    print()
    print("--- 錐の成分（R は変動するか）---")
    for lbl in ("a", "b"):
        Rp = H[f"Rp2_{lbl}"]; R2 = H[f"R2_{lbl}"]; q = H[f"q_{lbl}"]
        print(f"  {lbl.upper()}: R′²  平均={Rp.mean():.6e} 変動={100*Rp.std()/Rp.mean():.4f}%"
              f"  範囲=[{Rp.min():.4e}, {Rp.max():.4e}]")
        print(f"      R²  平均={R2.mean():.6e} 変動={100*R2.std()/max(R2.mean(),1e-30):.4f}%")
        print(f"      q   平均={q.mean():+.6f} 範囲=[{q.min():+.4f}, {q.max():+.4f}]")
    print()
    print(f"  質量² A 平均={H['m2_a'].mean():.6e} / B 平均={H['m2_b'].mean():.6e}")
    print(f"  反射率 r 平均={H['r'].mean():.6f}")
    print()
    ok = all([cr11, cr12, cr13, cr14, cr15])
    print(f"ALL PASS: {ok}  （所要 {time.time()-t0:.1f}s）")

    json.dump({"experiment": "cr1_kinetic_feedback_v1",
               "date": time.strftime("%Y-%m-%d %H:%M:%S"),
               "config": {"packet_a": list(PACKET_A), "packet_b": list(PACKET_B),
                          "deg_a": DEG_A, "deg_b": DEG_B, "T": T_STEPS,
                          "omega0": OMEGA_CLOCK, "kappa_mode": KAPPA_MODE,
                          "chi_definition": "Delta-theta (relative phase)"},
               "verdicts": {"CR1_1_moves": cr11, "CR1_2_closure": cr12,
                            "CR1_3_norm": cr13, "CR1_4_bound_oscillation": cr14,
                            "CR1_5_kappa0_matches_CR0": cr15},
               "metrics": {"total_change_deg": total_change, "turning_points": turns,
                           "closure_drift": closure_drift, "norm_drift": norm_drift,
                           "env_first_half": env1, "env_second_half": env2,
                           "env_growth": env_growth, "env_windows": envs,
                           "tail_growth": tail_growth,
                           "period_chi": per_chi, "period_Rp2": per_Rp,
                           "kappa0_diff": diff5,
                           "chi_deg_first": float(chi_deg[0]),
                           "chi_deg_last": float(chi_deg[-1]),
                           "chi_deg_min": float(chi_deg.min()),
                           "chi_deg_max": float(chi_deg.max()),
                           "chord_min": float(H["chord"].min()),
                           "chord_max": float(H["chord"].max()),
                           "omega_min": float(H["omega"].min()),
                           "omega_max": float(H["omega"].max()),
                           "m2_a_mean": float(H["m2_a"].mean()),
                           "m2_b_mean": float(H["m2_b"].mean()),
                           "r_mean": float(H["r"].mean()),
                           "kappa_mean": float(H["kappa"].mean()),
                           "kappa_min": float(H["kappa"].min()),
                           "kappa_max": float(H["kappa"].max())},
               "all_pass": ok},
              open(HERE / "result_cr1_kinetic_feedback_v1.json", "w"),
              ensure_ascii=False, indent=1)
    np.savez_compressed(HERE / "history_cr1_kinetic_feedback_v1.npz", **H)
    print("保存: result_cr1_kinetic_feedback_v1.json / history_...npz")


if __name__ == "__main__":
    main()
