#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CR0 予備: 荷電二体波の反跳と合成加速度の分離——θ書き込みなしの対照点

位置づけ
--------
「AB二体閉鎖位相系における未来位相位置加速度写像と調和閉鎖による逆二乗則」
（Concept DOI 10.5281/zenodo.21441081）の実験系を、統一万能関数の正本の上で
再構成する準備段階。本稿の実験は **θ を一度も書き込まない**。

なぜ書き込まないか: 現行の力学 F には位相前進（運動項）が存在しない。
弾性部 rotate_ab も媒介頂点も AB チャネル空間の **実回転** であり、
状態に e^{iφ} を掛けない。したがって合成位相は混合の副作用として振れるが
巻数は上がらない。この「巡らないこと」を先に対照点として記録しておかないと、
のちに運動項（調和閉鎖分散 ω_n = n·ω₁）を昇格させたとき、その寄与を
分離できなくなる。ゼロ点対照と同じ論理である。

構成（木原仕様）
----------------
  A: 等振幅倍音 1..17（局在性 高・最小位相セル 2π/17 = 21.2°）
  B: 等振幅倍音 1..3 （局在性 低・最小位相セル 2π/3  = 120°）
  円周上の初期位置: A = −30°, B = +30°（相対位相差 60°）
  観測用の波は置かない（直接読み出す）

力学は `collision_step_exact`（二体正本・再輸出）のみ。追加も改変もしない。

事前登録の判定
--------------
  CR0-1  合成位相 arg(Σa)・arg(Σb) の巻数が 0（＝巡らない）
         → 運動項の欠落を実証する。**これは失敗ではなく本実験の主目的**
  CR0-2  閉塞ドリフト < 1e-12（Σa²+Σb² の保存）
  CR0-3  ノルムドリフト < 1e-12
  CR0-4  巻き m スペクトル（η方向FFT）の主巻きが保存
         → 電荷 Q = m/3 が動かないこと
  CR0-5  φ = 2r·Im(conj(b)·a) は要素ごとに符号を持つ（単一符号に潰れない）

記録のみ（判定しない）
  Δθ(t)・弦 R″=2R·sin(Δθ/2)・r(t)・Pf/Pb・帯パリティ分率・巻きスペクトル

使い方: python3 run_cr0_control_no_theta_v1.py
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

# 力学は統一万能関数の再輸出のみを使う（恒久ルール）
_spec = importlib.util.spec_from_file_location("uni_cr0", UNI)
_uni = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = _uni
_spec.loader.exec_module(_uni)

collision_step_exact = _uni.collision_step_exact
base = _uni.two_body_base          # System A（状態生成の正本）
toy = _uni.two_body_v1.toy if hasattr(_uni.two_body_v1, "toy") else None

PACKET_A = tuple(range(1, 18))     # 1..17
PACKET_B = tuple(range(1, 4))      # 1..3
DEG_A = -30.0
DEG_B = +30.0
T_STEPS = 2000


# ------------------------------------------------------------------ 構成

def build_pair():
    """A・B を指定倍音・指定初期位置で構成する。

    等振幅は default_weights（全て 1.0）で既定。振幅を手で置かない。
    円周上の位置は packet shift（χ格子上の平行移動）で与える。
    """
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n_chi = int(sp.chi_grid_n)

    shift_a = DEG_A / 360.0 * n_chi
    shift_b = DEG_B / 360.0 * n_chi

    case = base.explicit_packet_case(
        mode="cr0_charged_two_body",
        packet_a=PACKET_A, packet_b=PACKET_B,
        packet_a_shift=shift_a, packet_b_shift=shift_b,
    )
    a = base.make_case_state(sp, case, "A", hair_enabled=True)
    b = base.make_case_state(sp, case, "B", hair_enabled=True)
    return a, b, sp, case


# ------------------------------------------------------------------ 読出し

def winding_spectrum(psi, n_chi, n_eta):
    """巻き m のスペクトルと符号付き巻き数（万能巻数と同一の符号規約）。

    二体トイでは η は毛円の **位置** であり、巻きは exp(i·m·η) として
    位相に乗る。したがって m を読むには η 方向の FFT が要る
    （多体エンジンの C2 は既に η がフーリエ添字なので、そこでは不要）。

    符号は剰余の折返しだけで作る（条件分岐を置かない）:
        m_signed = ((m + n_eta//2) % n_eta) - n_eta//2
    """
    f = np.fft.fft(psi.reshape(n_chi, n_eta), axis=1)
    P = np.sum(np.abs(f) ** 2, axis=0)
    P = P / P.sum()
    m = np.arange(n_eta)
    m_signed = ((m + n_eta // 2) % n_eta) - n_eta // 2
    return P, m_signed, float(P @ m_signed)


def band_split(psi, n_chi, n_eta):
    """帯 k の奇偶パワー（χ方向FFT）。搬送波込みの生スペクトルで測る。"""
    f = np.fft.fft(psi.reshape(n_chi, n_eta), axis=0)
    P = np.sum(np.abs(f) ** 2, axis=1)
    k = np.rint(np.fft.fftfreq(n_chi, d=1.0 / n_chi)).astype(int)
    odd = (np.abs(k) % 2) == 1
    po = float(P[odd].sum())
    pe = float(P[~odd].sum())
    return po, pe


def circle_position(psi, n_chi, n_eta):
    """円周上の位置 = パワー分布の第1円周モーメント（重心）。

    **arg(Σψ) は位置ではない**（v1 初走行の計器誤り）。合成振幅の偏角は
    搬送波 exp(i·q·p0·χ) の符号差をそのまま拾うため、q_A=+1・q_B=−1 の
    構成では初期差が 180° と出てしまい、指定した 60° を測れない。

    位置はレジスタ軸（χ）に住む。パワー P(χ)=Σ_η|ψ|² の第1モーメント
        z = Σ_χ P(χ)·e^{2πiχ/n}
    の偏角が円周上の重心位置であり、|z| が局在度（1 に近いほど局在）。
    これは統一読出しの位置規約（円上の値で返し巻数展開は呼出側）と同じ。
    """
    P = np.sum(np.abs(psi.reshape(n_chi, n_eta)) ** 2, axis=1)
    w = np.exp(2j * np.pi * np.arange(n_chi) / n_chi)
    z = complex(P @ w) / float(P.sum())
    return float(np.angle(z)), float(abs(z))


def resultant_phase(psi):
    """合成振幅の偏角 arg(Σψ)（位置ではない・参考記録用）。"""
    return float(np.angle(np.sum(psi)))


def unwrap_turns(series):
    """円上の位相列から総回転量（回転数）を求める。判定用。"""
    u = np.unwrap(np.asarray(series))
    return float((u[-1] - u[0]) / (2.0 * np.pi))


# ------------------------------------------------------------------ 本体

def main() -> None:
    t0 = time.time()
    a, b, sp, case = build_pair()
    n_chi, n_eta = int(sp.chi_grid_n), int(sp.eta_grid_n)

    tot0 = float(np.vdot(a, a).real + np.vdot(b, b).real)
    clo0 = abs(complex(np.sum(a * a) + np.sum(b * b)))

    Pm_a0, m_signed, q_a0 = winding_spectrum(a, n_chi, n_eta)
    Pm_b0, _, q_b0 = winding_spectrum(b, n_chi, n_eta)

    hist = {k: [] for k in
            ("theta", "r", "phi_min", "phi_max", "phi_mean", "phi_pos_frac",
             "arg_a", "arg_b", "dtheta", "chord", "closure", "norm",
             "q_a", "q_b", "q_tot", "pf_a", "pb_a", "pf_b", "pb_b",
             "pos_a", "pos_b", "loc_a", "loc_b", "res_a", "res_b")}

    for _ in range(T_STEPS):
        a, b, ro = collision_step_exact(a, b, sp)

        r = float(ro.reflection_rate)
        phi = 2.0 * r * np.imag(np.conj(b) * a)
        pa, la = circle_position(a, n_chi, n_eta)      # ★ 位置＝円周第1モーメント
        pb_, lb = circle_position(b, n_chi, n_eta)
        aa, bb = pa, pb_
        d = float(np.angle(np.exp(1j * (aa - bb))))     # 円上の差（巻かない）
        ra_, rb_ = resultant_phase(a), resultant_phase(b)   # 参考記録

        tot = float(np.vdot(a, a).real + np.vdot(b, b).real)
        clo = abs(complex(np.sum(a * a) + np.sum(b * b)))
        _, _, qa = winding_spectrum(a, n_chi, n_eta)
        _, _, qb = winding_spectrum(b, n_chi, n_eta)
        poa, pea = band_split(a, n_chi, n_eta)
        pob, peb = band_split(b, n_chi, n_eta)

        hist["theta"].append(float(ro.theta))
        hist["r"].append(r)
        hist["phi_min"].append(float(phi.min()))
        hist["phi_max"].append(float(phi.max()))
        hist["phi_mean"].append(float(phi.mean()))
        hist["phi_pos_frac"].append(float((phi > 0).mean()))
        hist["arg_a"].append(aa)
        hist["arg_b"].append(bb)
        hist["dtheta"].append(d)
        hist["chord"].append(float(2.0 * np.sin(abs(d) / 2.0)))   # R=1 規格
        hist["closure"].append(clo)
        hist["norm"].append(tot)
        hist["q_a"].append(qa)
        hist["q_b"].append(qb)
        hist["q_tot"].append(qa + qb)
        hist["pos_a"].append(pa); hist["pos_b"].append(pb_)
        hist["loc_a"].append(la); hist["loc_b"].append(lb)
        hist["res_a"].append(ra_); hist["res_b"].append(rb_)
        hist["pf_a"].append(poa); hist["pb_a"].append(pea)
        hist["pf_b"].append(pob); hist["pb_b"].append(peb)

    H = {k: np.asarray(v) for k, v in hist.items()}

    turns_a = unwrap_turns(H["arg_a"])
    turns_b = unwrap_turns(H["arg_b"])
    closure_drift = float(abs(H["closure"][-1] - clo0))
    norm_drift = float(np.max(np.abs(H["norm"] / tot0 - 1.0)))
    dq_a = float(np.max(np.abs(H["q_a"] - q_a0)))
    dq_b = float(np.max(np.abs(H["q_b"] - q_b0)))
    dq_tot = float(np.max(np.abs(H["q_tot"] - (q_a0 + q_b0))))
    both_signs = bool(np.any(H["phi_pos_frac"] > 0.0) and
                      np.any(H["phi_pos_frac"] < 1.0))

    cr01 = bool(abs(turns_a) < 1.0 and abs(turns_b) < 1.0)
    cr02 = bool(closure_drift <= 1e-12 * max(tot0, 1.0))
    cr03 = bool(norm_drift <= 1e-12)
    cr04 = bool(dq_tot <= 1e-9)     # 合計巻き（個別はチャネル間で移動する）
    cr05 = both_signs

    print("=" * 72)
    print("CR0 荷電二体波の反跳と合成加速度の分離 — θ書き込みなしの対照点")
    print("=" * 72)
    print(f"構成: A=倍音{PACKET_A[0]}..{PACKET_A[-1]}（{len(PACKET_A)}本）"
          f" / B=倍音{PACKET_B[0]}..{PACKET_B[-1]}（{len(PACKET_B)}本）")
    print(f"      初期位置 A={DEG_A}°  B={DEG_B}°  相対位相差={DEG_B-DEG_A}°")
    print(f"      χ格子={n_chi}  η格子={n_eta}  T={T_STEPS}")
    print(f"      搬送波 q_A={sp.q_A} q_B={sp.q_B} / 巻き m_A={sp.m_A} m_B={sp.m_B}")
    print()
    print("--- 初期の巻きスペクトル（電荷 Q=m/3 の読み） ---")
    for lbl, P in (("A", Pm_a0), ("B", Pm_b0)):
        top = int(np.argmax(P))
        print(f"  {lbl}: 主巻き m={m_signed[top]:+d}  純度={P[top]:.6f}"
              f"  符号付き巻数={(q_a0 if lbl=='A' else q_b0):+.6f}"
              f"  → Q={m_signed[top]/3.0:+.4f}"
              f"（{'可読' if m_signed[top] % 3 == 0 else '分数・単独不可読'}）")
    print()
    print("--- 判定 ---")
    print(f"  CR0-1 巡らない（総回転 A={turns_a:+.4f}周 B={turns_b:+.4f}周）: {cr01}")
    print(f"  CR0-2 閉塞ドリフト {closure_drift:.3e}: {cr02}")
    print(f"  CR0-3 ノルムドリフト {norm_drift:.3e}: {cr03}")
    print(f"  CR0-4 合計巻き保存 Δq_tot={dq_tot:.3e}: {cr04}"
          f"   （個別 Δq_A={dq_a:.3e} Δq_B={dq_b:.3e}＝チャネル間移動）")
    print(f"  CR0-5 φ が両符号を持つ: {cr05}")
    print()
    print("--- 記録（対照点） ---")
    print(f"  局在度|z|: A 平均={H['loc_a'].mean():.6f} B 平均={H['loc_b'].mean():.6f}")
    print(f"  合成振幅偏角の総回転（位置ではない・参考）: "
          f"A={unwrap_turns(H['res_a']):+.2f}周 B={unwrap_turns(H['res_b']):+.2f}周")
    print(f"  Δθ:  初期={np.degrees(H['dtheta'][0]):+.3f}°  "
          f"最終={np.degrees(H['dtheta'][-1]):+.3f}°  "
          f"範囲=[{np.degrees(H['dtheta'].min()):+.3f}, "
          f"{np.degrees(H['dtheta'].max()):+.3f}]°")
    print(f"  弦 R″/R: 平均={H['chord'].mean():.6f}  "
          f"範囲=[{H['chord'].min():.6f}, {H['chord'].max():.6f}]")
    print(f"  反射率 r: 平均={H['r'].mean():.6f}  "
          f"範囲=[{H['r'].min():.6f}, {H['r'].max():.6f}]")
    print(f"  θ: 平均={H['theta'].mean():.6f} rad")
    print(f"  φ: 最小={H['phi_min'].min():.3e} 最大={H['phi_max'].max():.3e}"
          f"  正の割合 平均={H['phi_pos_frac'].mean():.4f}")
    print(f"  帯パリティ分率 A: Pf/(Pf+Pb)="
          f"{(H['pf_a']/(H['pf_a']+H['pb_a'])).mean():.6f}")
    print(f"  帯パリティ分率 B: Pf/(Pf+Pb)="
          f"{(H['pf_b']/(H['pf_b']+H['pb_b'])).mean():.6f}")
    print()
    print(f"ALL PASS: {all([cr01, cr02, cr03, cr04, cr05])}"
          f"  （所要 {time.time()-t0:.1f}s）")

    out = {
        "experiment": "cr0_control_no_theta_v1",
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "config": {
            "packet_a": list(PACKET_A), "packet_b": list(PACKET_B),
            "deg_a": DEG_A, "deg_b": DEG_B, "T": T_STEPS,
            "chi_grid_n": n_chi, "eta_grid_n": n_eta,
            "q_A": float(sp.q_A), "q_B": float(sp.q_B),
            "m_A": int(sp.m_A), "m_B": int(sp.m_B),
        },
        "verdicts": {"CR0_1_no_winding": cr01, "CR0_2_closure": cr02,
                     "CR0_3_norm": cr03, "CR0_4_winding_conserved": cr04,
                     "CR0_5_phi_both_signs": cr05},
        "metrics": {
            "turns_a": turns_a, "turns_b": turns_b,
            "closure_drift": closure_drift, "norm_drift": norm_drift,
            "dq_a": dq_a, "dq_b": dq_b, "dq_tot": dq_tot,
            "q_a0": q_a0, "q_b0": q_b0,
            "dtheta_deg_first": float(np.degrees(H["dtheta"][0])),
            "dtheta_deg_last": float(np.degrees(H["dtheta"][-1])),
            "dtheta_deg_min": float(np.degrees(H["dtheta"].min())),
            "dtheta_deg_max": float(np.degrees(H["dtheta"].max())),
            "chord_mean": float(H["chord"].mean()),
            "chord_min": float(H["chord"].min()),
            "chord_max": float(H["chord"].max()),
            "r_mean": float(H["r"].mean()), "r_min": float(H["r"].min()),
            "r_max": float(H["r"].max()), "theta_mean": float(H["theta"].mean()),
            "phi_min": float(H["phi_min"].min()),
            "phi_max": float(H["phi_max"].max()),
            "phi_pos_frac_mean": float(H["phi_pos_frac"].mean()),
            "loc_a_mean": float(H["loc_a"].mean()),
            "loc_b_mean": float(H["loc_b"].mean()),
            "turns_resultant_a": unwrap_turns(H["res_a"]),
            "turns_resultant_b": unwrap_turns(H["res_b"]),
            "pf_frac_a": float((H["pf_a"] / (H["pf_a"] + H["pb_a"])).mean()),
            "pf_frac_b": float((H["pf_b"] / (H["pf_b"] + H["pb_b"])).mean()),
        },
        "winding_spectrum_initial": {
            "m_signed": m_signed.tolist(),
            "P_a": Pm_a0.tolist(), "P_b": Pm_b0.tolist(),
        },
    }
    (HERE / "result_cr0_control_no_theta_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    np.savez_compressed(HERE / "history_cr0_control_no_theta_v1.npz", **H)
    print("保存: result_cr0_control_no_theta_v1.json / history_...npz")


if __name__ == "__main__":
    main()
