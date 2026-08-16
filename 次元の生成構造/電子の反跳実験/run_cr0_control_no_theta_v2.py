#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CR0 v2: 荷電二体波の反跳——θ書き込みなしの対照点（較正版）

v1 からの変更（v1 の未解決3件を処理する）
------------------------------------------
(1) **shift の較正**: `packet_shift` は円周の度ではない。v1 では −30 を
    指定してピークが −104.766° に着地し、実測 Δθ が 150.76°（指定 60°）
    になっていた。本版は **走行前に2点で実測して線形係数を求め**、
    指定角に着地する shift を逆算する。係数は直書きしない。
(2) **局在度は参加率 PR で測る**: PR = 1/Σ_χ P(χ)²（実効セル数）。
    v1 の第1モーメント |z| は混合で A・B とも 0.551 に潰れ、
    設定した局在の非対称を追えなかった（T=0 では 0.941/0.667 と差が出る）。
    |z| は副次指標として併記する。
(3) **巻きスペクトルの時間発展を保存**: v1 で合計巻きが 2.28e-02 動いた
    原因を切り分ける。二つの候補を同時に検定する。
      (a) 和則ウォーク: m* = 2·m_pump − m_seed。m_A=1, m_B=2 なら
          2·1−2 = 0 と 2·2−1 = 3 が生成される。P(0)・P(3) の成長を見る。
      (b) η の Nyquist 折返し: ne=16 の端ビン（m=−8, +7）にパワーが
          溜まるか。v8 で「巻き数の真の保存則は mod ne の巡回保存」と
          同定された前例がある。

構成（木原仕様・v1 と同一）
---------------------------
  A: 等振幅倍音 1..17（局在性 高）  B: 等振幅倍音 1..3（局在性 低）
  円周上の初期位置: A = −30°, B = +30°（相対位相差 60°）
  観測用の波は置かない。力学は collision_step_exact のみ。θ は書き込まない。

事前登録の判定
--------------
  CR0-1  位置が巡らない（総回転 |turns| < 1）
  CR0-2  閉塞ドリフト < 1e-12
  CR0-3  ノルムドリフト < 1e-12
  CR0-5  φ が両符号を持つ
  CR0-6  **較正後の初期 Δθ が 60° ± 1°**（v2 で追加）
  CR0-7  **合計巻きの変化が、和則生成で説明できる**（v2 で追加）
         判定: ΔP(m=0)+ΔP(m=3) が |Δ| 総量の 50% 以上を占める
         → 和則ウォーク支持。端ビン優位なら折返し支持。

使い方: python3 run_cr0_control_no_theta_v2.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
UNI = HERE.parent / "統一万能関数_v1" / "unified_interaction_v1.py"

_spec = importlib.util.spec_from_file_location("uni_cr0v2", UNI)
_uni = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = _uni
_spec.loader.exec_module(_uni)

collision_step_exact = _uni.collision_step_exact
base = _uni.two_body_base

PACKET_A = tuple(range(1, 18))
PACKET_B = tuple(range(1, 4))
DEG_A = -30.0
DEG_B = +30.0
T_STEPS = 2000
SPEC_EVERY = 10          # 巻きスペクトルの保存間隔


# ------------------------------------------------------------------ 読出し

def power_chi(psi, n_chi, n_eta):
    P = np.sum(np.abs(psi.reshape(n_chi, n_eta)) ** 2, axis=1)
    return P / P.sum()


def peak_deg(psi, n_chi, n_eta):
    """パワー分布のピーク位置（度・−180..180）。較正に使う。"""
    P = power_chi(psi, n_chi, n_eta)
    d = float(np.argmax(P)) * 360.0 / n_chi
    return d - 360.0 if d > 180.0 else d


def circle_position(psi, n_chi, n_eta):
    """円周第1モーメント: (重心角[rad], |z|)。"""
    P = power_chi(psi, n_chi, n_eta)
    w = np.exp(2j * np.pi * np.arange(n_chi) / n_chi)
    z = complex(P @ w)
    return float(np.angle(z)), float(abs(z))


def participation_ratio(psi, n_chi, n_eta):
    """参加率 PR = 1/Σ P² （実効セル数）。局在度の正規の計器。"""
    P = power_chi(psi, n_chi, n_eta)
    return float(1.0 / np.sum(P ** 2))


def winding_spectrum(psi, n_chi, n_eta):
    """巻き m のスペクトル P(m)（規格化）と符号付き巻き数。"""
    f = np.fft.fft(psi.reshape(n_chi, n_eta), axis=1)
    P = np.sum(np.abs(f) ** 2, axis=0)
    P = P / P.sum()
    m = np.arange(n_eta)
    m_signed = ((m + n_eta // 2) % n_eta) - n_eta // 2
    return P, m_signed, float(P @ m_signed)


def band_split(psi, n_chi, n_eta):
    f = np.fft.fft(psi.reshape(n_chi, n_eta), axis=0)
    P = np.sum(np.abs(f) ** 2, axis=1)
    k = np.rint(np.fft.fftfreq(n_chi, d=1.0 / n_chi)).astype(int)
    odd = (np.abs(k) % 2) == 1
    return float(P[odd].sum()), float(P[~odd].sum())


def unwrap_turns(series):
    u = np.unwrap(np.asarray(series))
    return float((u[-1] - u[0]) / (2.0 * np.pi))


# ------------------------------------------------------------------ 較正

def make_pair(sp, shift_a, shift_b):
    case = base.explicit_packet_case(
        mode="cr0v2_charged_two_body",
        packet_a=PACKET_A, packet_b=PACKET_B,
        packet_a_shift=shift_a, packet_b_shift=shift_b,
    )
    a = base.make_case_state(sp, case, "A", hair_enabled=True)
    b = base.make_case_state(sp, case, "B", hair_enabled=True)
    return a, b


def calibrate_shift(sp, n_chi, n_eta, n_probe=17, span=float(np.pi) / 2.0):
    """shift → 円周角 の写像を実測で決める（係数を直書きしない）。

    実測の結論: **shift はラジアンである**。
        重心角 = π + shift   （slope = 180/π = 57.295780 度/単位、1周 = 2π）
    ただしこれを直書きせず毎回測る。掃引幅は 1/4 周（π/2）に取る——
    初走行で span=128（約20周）を17点で刻み、1点あたり 458° 進んで
    エイリアシングを起こし、掃引点上だけ残差 2e-13 で「完璧に」合う
    偽の傾き 12.2958 を得た（真の傾きとの差はちょうど 360/8）。
    較正は Δθ を測るのと同じ計器（重心）で行う。argmax は等振幅パケットの
    副ローブ間で跳ぶため使えない。
    """
    xs = np.linspace(0.0, span, n_probe)
    ys = []
    for s in xs:
        a, _ = make_pair(sp, s, s)
        # 較正は **重心**（Δθ を測るのと同じ計器）で行う。argmax は
        # 等振幅パケットの副ローブ間で跳ぶため較正に使えない（v2 初走行の失敗）。
        ys.append(circle_position(a, n_chi, n_eta)[0])
    yd = np.degrees(np.unwrap(np.asarray(ys)))
    slope, intercept = np.polyfit(xs, yd, 1)
    resid = float(np.max(np.abs(yd - (slope * xs + intercept))))
    return float(slope), float(intercept), resid


def shift_for_deg(deg, slope, intercept):
    """指定角に着地する shift（|s| が最小の枝を返す）。"""
    s0 = (deg - intercept) / slope
    period = 360.0 / slope
    return float(s0 - np.round(s0 / period) * period)


# ------------------------------------------------------------------ 本体

def main() -> None:
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n_chi, n_eta = int(sp.chi_grid_n), int(sp.eta_grid_n)

    slope, intercept, cal_resid = calibrate_shift(sp, n_chi, n_eta)
    shift_a = shift_for_deg(DEG_A, slope, intercept)
    shift_b = shift_for_deg(DEG_B, slope, intercept)

    a, b = make_pair(sp, shift_a, shift_b)

    pk_a0, pk_b0 = peak_deg(a, n_chi, n_eta), peak_deg(b, n_chi, n_eta)
    pos_a0, _ = circle_position(a, n_chi, n_eta)
    pos_b0, _ = circle_position(b, n_chi, n_eta)
    dtheta0 = float(np.degrees(np.angle(np.exp(1j * (pos_a0 - pos_b0)))))
    pr_a0 = participation_ratio(a, n_chi, n_eta)
    pr_b0 = participation_ratio(b, n_chi, n_eta)

    Pm_a0, m_signed, q_a0 = winding_spectrum(a, n_chi, n_eta)
    Pm_b0, _, q_b0 = winding_spectrum(b, n_chi, n_eta)

    tot0 = float(np.vdot(a, a).real + np.vdot(b, b).real)
    clo0 = abs(complex(np.sum(a * a) + np.sum(b * b)))

    keys = ("theta", "r", "phi_min", "phi_max", "phi_pos_frac",
            "pos_a", "pos_b", "dtheta", "chord", "closure", "norm",
            "q_a", "q_b", "q_tot", "pr_a", "pr_b", "z_a", "z_b",
            "pf_a", "pb_a", "pf_b", "pb_b")
    hist = {k: [] for k in keys}
    spec_a, spec_b, spec_t = [], [], []

    for j in range(T_STEPS):
        a, b, ro = collision_step_exact(a, b, sp)

        r = float(ro.reflection_rate)
        phi = 2.0 * r * np.imag(np.conj(b) * a)
        pa, za = circle_position(a, n_chi, n_eta)
        pb_, zb = circle_position(b, n_chi, n_eta)
        d = float(np.angle(np.exp(1j * (pa - pb_))))
        Pa, _, qa = winding_spectrum(a, n_chi, n_eta)
        Pb, _, qb = winding_spectrum(b, n_chi, n_eta)
        poa, pea = band_split(a, n_chi, n_eta)
        pob, peb = band_split(b, n_chi, n_eta)

        hist["theta"].append(float(ro.theta)); hist["r"].append(r)
        hist["phi_min"].append(float(phi.min()))
        hist["phi_max"].append(float(phi.max()))
        hist["phi_pos_frac"].append(float((phi > 0).mean()))
        hist["pos_a"].append(pa); hist["pos_b"].append(pb_)
        hist["dtheta"].append(d)
        hist["chord"].append(float(2.0 * np.sin(abs(d) / 2.0)))
        hist["closure"].append(abs(complex(np.sum(a * a) + np.sum(b * b))))
        hist["norm"].append(float(np.vdot(a, a).real + np.vdot(b, b).real))
        hist["q_a"].append(qa); hist["q_b"].append(qb)
        hist["q_tot"].append(qa + qb)
        hist["pr_a"].append(participation_ratio(a, n_chi, n_eta))
        hist["pr_b"].append(participation_ratio(b, n_chi, n_eta))
        hist["z_a"].append(za); hist["z_b"].append(zb)
        hist["pf_a"].append(poa); hist["pb_a"].append(pea)
        hist["pf_b"].append(pob); hist["pb_b"].append(peb)

        if j % SPEC_EVERY == 0:
            spec_a.append(Pa); spec_b.append(Pb); spec_t.append(j)

    H = {k: np.asarray(v) for k, v in hist.items()}
    SA = np.asarray(spec_a); SB = np.asarray(spec_b)

    # --- 判定 ---
    turns_a = unwrap_turns(H["pos_a"]); turns_b = unwrap_turns(H["pos_b"])
    closure_drift = float(abs(H["closure"][-1] - clo0))
    norm_drift = float(np.max(np.abs(H["norm"] / tot0 - 1.0)))
    dq_tot = float(np.max(np.abs(H["q_tot"] - (q_a0 + q_b0))))
    both_signs = bool(np.any(H["phi_pos_frac"] > 0.0)
                      and np.any(H["phi_pos_frac"] < 1.0))

    # CR0-7: 巻きの移動先を和則と折返しで分解する
    dPa = SA[-1] - SA[0]; dPb = SB[-1] - SB[0]
    dP = np.abs(dPa) + np.abs(dPb)
    tot_move = float(dP.sum())
    idx = {int(ms): i for i, ms in enumerate(m_signed)}
    init_m = [int(sp.m_A), int(sp.m_B)]      # 初期占有（A↔B のチャネル間移動）
    sumrule_m = [2 * int(sp.m_A) - int(sp.m_B), 2 * int(sp.m_B) - int(sp.m_A)]
    nyquist_m = [int(m_signed.min()), int(m_signed.max())]
    # 「新規生成」= 初期に占有されていなかった巻き
    new_mask = np.array([int(ms) not in init_m for ms in m_signed])
    new_move = float(dP[new_mask].sum())
    sr = float(sum(dP[idx[m]] for m in sumrule_m if m in idx))
    ny = float(sum(dP[idx[m]] for m in nyquist_m if m in idx))
    sr_frac = sr / new_move if new_move > 0 else float("nan")
    ny_frac = ny / new_move if new_move > 0 else float("nan")
    exch = tot_move - new_move

    cr01 = bool(abs(turns_a) < 1.0 and abs(turns_b) < 1.0)
    cr02 = bool(closure_drift <= 1e-12 * max(tot0, 1.0))
    cr03 = bool(norm_drift <= 1e-12)
    cr05 = both_signs
    cr06 = bool(abs(abs(dtheta0) - abs(DEG_B - DEG_A)) <= 1.0)
    cr07 = bool(sr_frac >= 0.5)

    print("=" * 74)
    print("CR0 v2 荷電二体波の反跳 — θ書き込みなしの対照点（較正版）")
    print("=" * 74)
    print(f"構成: A=倍音1..17（{len(PACKET_A)}本） / B=倍音1..3（{len(PACKET_B)}本）"
          f"  χ={n_chi} η={n_eta} T={T_STEPS}")
    print(f"      搬送波 q_A={sp.q_A} q_B={sp.q_B} / 巻き m_A={sp.m_A} m_B={sp.m_B}")
    print()
    print("--- shift の較正（掃引17点・実測・直書きなし） ---")
    print(f"  角度 = {slope:.6f}·shift + {intercept:+.4f}  (mod 360)")
    print(f"  最大残差 = {cal_resid:.4e}°   1周に必要な shift = {360/slope:.4f}")
    print(f"  → 指定 A={DEG_A}° B={DEG_B}° に対し "
          f"shift_A={shift_a:.4f} shift_B={shift_b:.4f}")
    print(f"  着地: A ピーク={pk_a0:+.4f}° 重心={np.degrees(pos_a0):+.4f}°"
          f" / B ピーク={pk_b0:+.4f}° 重心={np.degrees(pos_b0):+.4f}°")
    print(f"  初期 Δθ = {dtheta0:+.4f}°（指定 {DEG_A-DEG_B:+.1f}°）")
    print()
    print("--- 局在度（T=0） ---")
    print(f"  A: PR={pr_a0:8.2f}セル ({100*pr_a0/n_chi:5.2f}%)   B: PR={pr_b0:8.2f}セル"
          f" ({100*pr_b0/n_chi:5.2f}%)   比={pr_b0/pr_a0:.2f}倍")
    print()
    print("--- 判定 ---")
    print(f"  CR0-1 巡らない（A={turns_a:+.4f}周 B={turns_b:+.4f}周）: {cr01}")
    print(f"  CR0-2 閉塞ドリフト {closure_drift:.3e}: {cr02}")
    print(f"  CR0-3 ノルムドリフト {norm_drift:.3e}: {cr03}")
    print(f"  CR0-5 φ が両符号: {cr05}")
    print(f"  CR0-6 初期Δθ が {abs(DEG_B-DEG_A):.0f}°±1°: {cr06}")
    print(f"  CR0-7 新規巻きを和則が説明（和則 {sr_frac:.4f} / 端ビン {ny_frac:.4f}）: {cr07}")
    print(f"        新規生成 {new_move:.6f} / チャネル間交換 {exch:.6f}"
          f" / 和則予測 m*={sumrule_m}")
    print()
    print("--- 巻きの移動（P(m) の T=0 → T=末 の変化・A+B 合算） ---")
    order = np.argsort(m_signed)
    for i in order:
        m = int(m_signed[i]); dv = float(dP[i])
        tag = ""
        if m in init_m: tag = " ← 初期占有（交換）"
        if m in sumrule_m: tag = " ← 和則 m*=2m_p−m_s"
        if m in nyquist_m: tag = " ← 端ビン(Nyquist)"
        if abs(dv) > 1e-6:
            print(f"    m={m:+3d}  |ΔP|={dv:.6f}{tag}")
    print(f"  合計 |ΔP| = {tot_move:.6f}   合計巻き変化 Δq_tot={dq_tot:.3e}")
    print()
    print("--- 記録 ---")
    print(f"  Δθ: 初期={np.degrees(H['dtheta'][0]):+.3f}° "
          f"最終={np.degrees(H['dtheta'][-1]):+.3f}° "
          f"範囲=[{np.degrees(H['dtheta'].min()):+.3f}, "
          f"{np.degrees(H['dtheta'].max()):+.3f}]°")
    print(f"  弦 R″/R: 平均={H['chord'].mean():.6f} "
          f"範囲=[{H['chord'].min():.6f}, {H['chord'].max():.6f}]")
    print(f"  PR: A 平均={H['pr_a'].mean():.2f} (T=0 {pr_a0:.2f}) / "
          f"B 平均={H['pr_b'].mean():.2f} (T=0 {pr_b0:.2f})")
    print(f"  |z|: A 平均={H['z_a'].mean():.6f} / B 平均={H['z_b'].mean():.6f}")
    print(f"  反射率 r: 平均={H['r'].mean():.6f} "
          f"範囲=[{H['r'].min():.6f}, {H['r'].max():.6f}]")
    print(f"  帯パリティ Pf/(Pf+Pb): A={(H['pf_a']/(H['pf_a']+H['pb_a'])).mean():.6f}"
          f" B={(H['pf_b']/(H['pf_b']+H['pb_b'])).mean():.6f}")
    print()
    print(f"ALL PASS: {all([cr01, cr02, cr03, cr05, cr06, cr07])}"
          f"  （所要 {time.time()-t0:.1f}s）")

    out = {
        "experiment": "cr0_control_no_theta_v2",
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "config": {"packet_a": list(PACKET_A), "packet_b": list(PACKET_B),
                   "deg_a": DEG_A, "deg_b": DEG_B, "T": T_STEPS,
                   "chi_grid_n": n_chi, "eta_grid_n": n_eta,
                   "q_A": float(sp.q_A), "q_B": float(sp.q_B),
                   "m_A": int(sp.m_A), "m_B": int(sp.m_B)},
        "calibration": {"slope_deg_per_shift": slope, "intercept_deg": intercept,
                        "max_residual_deg": cal_resid,
                        "shift_a": shift_a, "shift_b": shift_b,
                        "peak_a0": pk_a0, "peak_b0": pk_b0,
                        "dtheta0_deg": dtheta0},
        "verdicts": {"CR0_1_no_winding": cr01, "CR0_2_closure": cr02,
                     "CR0_3_norm": cr03, "CR0_5_phi_both_signs": cr05,
                     "CR0_6_dtheta_60": cr06, "CR0_7_sumrule": cr07},
        "metrics": {
            "turns_a": turns_a, "turns_b": turns_b,
            "closure_drift": closure_drift, "norm_drift": norm_drift,
            "dq_tot": dq_tot, "q_a0": q_a0, "q_b0": q_b0,
            "pr_a0": pr_a0, "pr_b0": pr_b0,
            "pr_a_mean": float(H["pr_a"].mean()),
            "pr_b_mean": float(H["pr_b"].mean()),
            "z_a_mean": float(H["z_a"].mean()),
            "z_b_mean": float(H["z_b"].mean()),
            "sumrule_frac": sr_frac, "nyquist_frac": ny_frac,
            "total_dP": tot_move, "new_dP": new_move, "exchange_dP": exch,
            "sumrule_m": sumrule_m,
            "dtheta_deg_first": float(np.degrees(H["dtheta"][0])),
            "dtheta_deg_last": float(np.degrees(H["dtheta"][-1])),
            "dtheta_deg_min": float(np.degrees(H["dtheta"].min())),
            "dtheta_deg_max": float(np.degrees(H["dtheta"].max())),
            "chord_mean": float(H["chord"].mean()),
            "chord_min": float(H["chord"].min()),
            "chord_max": float(H["chord"].max()),
            "r_mean": float(H["r"].mean()),
            "pf_frac_a": float((H["pf_a"]/(H["pf_a"]+H["pb_a"])).mean()),
            "pf_frac_b": float((H["pf_b"]/(H["pf_b"]+H["pb_b"])).mean()),
        },
        "winding": {"m_signed": m_signed.tolist(),
                    "P_a_initial": Pm_a0.tolist(), "P_b_initial": Pm_b0.tolist(),
                    "P_a_final": SA[-1].tolist(), "P_b_final": SB[-1].tolist(),
                    "dP_abs_sum": dP.tolist()},
    }
    (HERE / "result_cr0_control_no_theta_v2.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    np.savez_compressed(HERE / "history_cr0_control_no_theta_v2.npz",
                        spec_t=np.asarray(spec_t), spec_a=SA, spec_b=SB,
                        m_signed=m_signed, **H)
    print("保存: result_cr0_control_no_theta_v2.json / history_...v2.npz")


if __name__ == "__main__":
    main()
