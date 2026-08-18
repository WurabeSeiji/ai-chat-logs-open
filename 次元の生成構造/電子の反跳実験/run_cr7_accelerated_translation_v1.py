#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CR7: 加速する並進——エネルギーを増やすと R′ は動くか

狙い
----
CR6（並進のみ・等速）では、保存量がすべてビット単位で不変だった
（Σ|a|²+Σ|b|²・R²・R′a²・R′b²・|z_A| のいずれも初期からの差 0.000e+00）。

そこで速度に一定増分の加速度を与える。1 步あたりの並進量を

    step(t) = STEP0 + ACCEL·t   ［度］（A に +、B に −）

とする。速度が増えるので、系のエネルギーが増えれば R′ が動くはずである。

予想（CR6 の実測にもとづく）: **R′ は動かない。**
並進は k 空間の位相ランプ A_k → A_k·e^{ikω} であり、これは複素関数の剛体的な
平行移動にすぎない。Σ(Re a)² は添字の付け替えで不変なので、ω をどれだけ
大きくしても値は変わらない。

  → 予想どおりなら、現在の運動項は**運動エネルギーを状態に入れていない**。
    速度は波に担われず、毎步外から与える外部パラメータになっている。
  → 予想が外れて R′ が動くなら、並進の実装が単なる平行移動ではない。

どちらでも所見になる。これが本実験の意味。

構成（CR6 からの差分は並進量だけ）
----------------------------------
  (1) 初期化: A・B とも倍音 1〜17 の等振幅（完全な局在性）。
      初期位置は A = −30°, B = +30°。CR6 と同一。
  (2) 相互作用: collision_step_exact を**呼ばない**（本文中でコメントアウト）。
  (3) 並進: 1 步あたり A に +(STEP0+ACCEL·t)°, B に −(STEP0+ACCEL·t)°。
      初期 ±1°/步、1 步ごとに ±0.1°/步 ずつ増える。

記録するもの
------------
  cloRe / cloIm : Σa²+Σb²（非共役）の実部・虚部
  clo_a / clo_b : 各体の Σa²・Σb²（非共役）の実部・虚部
  Rp2A / Rp2B   : R′² = Σx²（複素にまとめる前の実成分の二乗和）
  R2            : R′a² + R′b²
  pqA / pqB     : Σp·q（閉包の虚部条件の残差）
  norm          : Σ|a|²+Σ|b|²（エネルギー相当・一定のはず）
  chi / turns   : 相対位相と累積巻数
  波形・包絡     : CR5 と同じ形式（図はそのまま流用できる）

出力: cr7_accelerated_data_v1.json / cr7_accelerated_v1.html

使い方: python3 run_cr7_accelerated_translation_v1.py
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
SHELL = HERE / "cr_waveform_shell_v1.html"


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


_uni = _load("uni_cr7", UNI / "unified_interaction_v1.py")
K = _load("kin_cr7", UNI / "unified_kinetic_v1.py")
_cr0 = _load("cr0_cr7", HERE / "run_cr0_control_no_theta_v2.py")
_cr1 = _load("cr1_cr7", HERE / "run_cr1_kinetic_feedback_v1.py")

PACKET = tuple(range(1, 18))      # A・B とも倍音 1〜17（等振幅・完全な局在性）
DEG_A, DEG_B = -30.0, +30.0       # 初期位置は現行のまま
STEP0 = 1.0                       # 初期の並進量［度/步］。A に +、B に −
ACCEL = 0.1                       # 1 步ごとの増分［度/步²］（一定加速度）
T_STEPS = 400
CHI_KEEP = 256


def centered_profiles(psi, n_chi, n_eta, keep):
    """CR5 と同一（η=0 スライスを中心位相 0° に揃える）。"""
    A2 = psi.reshape(n_chi, n_eta)
    env = np.sqrt(np.sum(np.abs(A2) ** 2, axis=1))
    wave_c = A2[:, 0].copy()
    ang, _ = _cr0.circle_position(psi, n_chi, n_eta)
    shift = int(np.rint(ang / (2.0 * np.pi) * n_chi))
    env = np.roll(env, -shift)
    wave_c = np.roll(wave_c, -shift)
    carrier = float(np.angle(wave_c[0])) if abs(wave_c[0]) > 0 else 0.0
    wave_c = wave_c * np.exp(-1j * carrier)
    env = np.roll(env, n_chi // 2)
    wave_c = np.roll(wave_c, n_chi // 2)
    step = n_chi // keep
    return wave_c.real[::step], env[::step], carrier


def sq(psi):
    """(非共役二乗和 Σz², R′²=Σx², Σp·q)。"""
    z = complex(np.sum(psi * psi))
    p, q = psi.real, psi.imag
    return z, float(np.sum(p * p)), float(np.sum(p * q))


def calibrate_omega1(sp, n_chi, n_eta, slope, icept, deg):
    """指定角［度］/步 の並進を与える omega1 を実測で決める（符号も含めて）。

    k_translate_flat の符号規約を直書きせず、1 步進めて重心角の変化を測る。
    """
    case = _uni.two_body_base.explicit_packet_case(
        mode="cr7_cal", packet_a=PACKET, packet_b=PACKET,
        packet_a_shift=_cr0.shift_for_deg(0.0, slope, icept),
        packet_b_shift=_cr0.shift_for_deg(0.0, slope, icept))
    s = _uni.two_body_base.make_case_state(sp, case, "A", hair_enabled=True)
    probe = np.radians(deg)
    p0, _ = _cr0.circle_position(s, n_chi, n_eta)
    p1, _ = _cr0.circle_position(
        K.k_translate_flat(s, probe, n_chi, n_eta), n_chi, n_eta)
    moved = float(np.degrees(np.angle(np.exp(1j * (p1 - p0)))))
    sign = 1.0 if moved > 0 else -1.0
    return probe * sign, abs(moved)


def main() -> None:
    t0 = time.time()
    sp = _uni.two_body_base.build_source_params(
        _uni.two_body_base.Params(high_n=63, recursive_collision_count=200))
    n_chi, n_eta = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = _cr0.calibrate_shift(sp, n_chi, n_eta)

    om1_per_deg, moved = calibrate_omega1(sp, n_chi, n_eta, slope, icept, 1.0)
    print(f"格子 n_chi={n_chi} n_eta={n_eta}")
    print(f"並進の較正: omega1={om1_per_deg:+.8f} rad で 1 步 {moved:.6f}° 動く "
          f"（指定 1.0°、差 {abs(moved-1.0):.2e}）")
    print(f"並進量: step(t) = {STEP0} + {ACCEL}·t ［度/步］"
          f"  → t=0 で {STEP0}°、t={T_STEPS-1} で {STEP0+ACCEL*(T_STEPS-1):.1f}°")

    case = _uni.two_body_base.explicit_packet_case(
        mode="cr7_accelerated", packet_a=PACKET, packet_b=PACKET,
        packet_a_shift=_cr0.shift_for_deg(DEG_A, slope, icept),
        packet_b_shift=_cr0.shift_for_deg(DEG_B, slope, icept))
    a = _uni.two_body_base.make_case_state(sp, case, "A", hair_enabled=True)
    b = _uni.two_body_base.make_case_state(sp, case, "B", hair_enabled=True)

    theta0 = float(_cr1.toy.theta_from_ab(a, b, sp).theta)   # 記録のみ・適用しない
    H = {k: [] for k in ("chi", "posA", "posB", "r_readonly",
                         "cloRe", "cloIm", "cloReA", "cloImA", "cloReB", "cloImB",
                         "Rp2A", "Rp2B", "R2", "pqA", "pqB", "norm", "za", "zb",
                         "step_deg")}
    WA, WB, EA, EB = [], [], [], []

    for t in range(T_STEPS):
        pa, za = _cr0.circle_position(a, n_chi, n_eta)
        pb, zb = _cr0.circle_position(b, n_chi, n_eta)
        H["posA"].append(float(pa)); H["posB"].append(float(pb))
        H["chi"].append(float(np.angle(np.exp(1j * (pa - pb)))))
        # θ は読むだけ。回転は適用しない（相互作用を切っているため）
        H["r_readonly"].append(
            round(float(_cr1.toy.theta_from_ab(a, b, sp).reflection_rate), 8))

        za_, ra2, pqa = sq(a)
        zb_, rb2, pqb = sq(b)
        Z = za_ + zb_
        H["cloRe"].append(float(f"{Z.real:.6g}")); H["cloIm"].append(float(f"{Z.imag:.6g}"))
        H["cloReA"].append(float(f"{za_.real:.6g}")); H["cloImA"].append(float(f"{za_.imag:.6g}"))
        H["cloReB"].append(float(f"{zb_.real:.6g}")); H["cloImB"].append(float(f"{zb_.imag:.6g}"))
        H["Rp2A"].append(round(ra2, 10)); H["Rp2B"].append(round(rb2, 10))
        H["R2"].append(round(ra2 + rb2, 10))
        H["pqA"].append(float(f"{pqa:.6g}")); H["pqB"].append(float(f"{pqb:.6g}"))
        H["norm"].append(round(float(np.vdot(a, a).real + np.vdot(b, b).real), 12))
        H["za"].append(round(za, 6)); H["zb"].append(round(zb, 6))

        wa, ea, _ = centered_profiles(a, n_chi, n_eta, CHI_KEEP)
        wb, eb, _ = centered_profiles(b, n_chi, n_eta, CHI_KEEP)
        WA.append(wa); EA.append(ea); WB.append(wb); EB.append(eb)

        # --- 並進のみ。A に +step(t)、B に −step(t)（一定加速度）---
        step_t = STEP0 + ACCEL * t
        H["step_deg"].append(round(step_t, 6))
        om_t = om1_per_deg * step_t
        a = K.k_translate_flat(a, +om_t, n_chi, n_eta)
        b = K.k_translate_flat(b, -om_t, n_chi, n_eta)
        # --- 相互作用は無効化（対照点） ---
        # a, b, _ = _uni.collision_step_exact(a, b, sp)

    WA = np.array(WA); WB = np.array(WB); EA = np.array(EA); EB = np.array(EB)
    scale = float(max(EA.max(), EB.max()))
    wscale = float(max(np.abs(WA).max(), np.abs(WB).max()))
    rnd = lambda X, s_: [[round(float(x), 4) for x in row] for row in (X / s_)]

    chi = np.array(H["chi"])
    turns = float((np.unwrap(chi)[-1] - np.unwrap(chi)[0]) / (2 * np.pi))

    def rng(k):
        v = np.array(H[k], float)
        return float(v.min()), float(v.max()), float(np.max(np.abs(v - v[0])))

    print(f"\n【判定】相互作用なし・加速する並進 T={T_STEPS}")
    print(f"  Δθ: {np.degrees(chi[0]):+.3f}° → 累積 {turns:+.4f} 周 "
          f"（予測 {sum(2*(STEP0+ACCEL*t) for t in range(T_STEPS))/360.0:+.4f} 周）")
    print(f"  {'量':14} {'初期':>14} {'最小':>14} {'最大':>14} {'初期からの最大差':>16}")
    for k, lab in (("norm", "Σ|a|²+Σ|b|²"), ("R2", "R²=R′a²+R′b²"),
                   ("Rp2A", "R′a²"), ("Rp2B", "R′b²"),
                   ("cloRe", "Σz² 実部"), ("cloIm", "Σz² 虚部"),
                   ("cloReA", "Σa² 実部"), ("cloImA", "Σa² 虚部"),
                   ("pqA", "Σp·q (A)"), ("za", "|z_A|")):
        lo, hi, dmax = rng(k)
        print(f"  {lab:14} {H[k][0]:14.6e} {lo:14.6e} {hi:14.6e} {dmax:16.3e}")

    out = {
        "experiment": "cr7_accelerated_translation_v1", "case": "accelerated",
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "config": {"packet_a": list(PACKET), "packet_b": list(PACKET),
                   "label": "A・B ともに倍音1〜17（等振幅）",
                   "note": f"相互作用なし・加速する並進（step={STEP0}+{ACCEL}·t 度/步）",
                   "deg_a": DEG_A, "deg_b": DEG_B, "T": T_STEPS, "every": 1,
                   "step0": STEP0, "accel": ACCEL, "omega1_per_deg": om1_per_deg,
                   "collision": False, "kappa": None,
                   "chi_keep": CHI_KEEP, "n_chi": n_chi, "n_eta": n_eta,
                   "theta0": theta0, "r0": H["r_readonly"][0],
                   "period_pred": None, "turns": turns,
                   "norm_scale": scale, "wave_scale": wscale},
        "chi": [round(x, 8) for x in H["chi"]],
        "dchi": [float(f"{x:.6g}") for x in np.gradient(chi)],
        "mid": [0.0] * T_STEPS,
        "r": H["r_readonly"], "v": [0.0] * T_STEPS, "omega": [0.0] * T_STEPS,
        "cloRe": H["cloRe"], "cloIm": H["cloIm"],
        "cloReA": H["cloReA"], "cloImA": H["cloImA"],
        "cloReB": H["cloReB"], "cloImB": H["cloImB"],
        "Rp2A": H["Rp2A"], "Rp2B": H["Rp2B"], "R2": H["R2"],
        "pqA": H["pqA"], "pqB": H["pqB"], "norm": H["norm"],
        "step_deg": H["step_deg"],
        "posA": [round(x, 8) for x in H["posA"]],
        "posB": [round(x, 8) for x in H["posB"]],
        "PA": [round(x, 8) for x in H["Rp2A"]],
        "PB": [round(x, 8) for x in H["Rp2B"]],
        "za": H["za"], "zb": H["zb"],
        "waveA": rnd(WA, wscale), "waveB": rnd(WB, wscale),
        "envA": rnd(EA, scale), "envB": rnd(EB, scale),
    }
    p = HERE / "cr7_accelerated_data_v1.json"
    data = json.dumps(out, ensure_ascii=False, separators=(",", ":"))
    p.write_text(data, encoding="utf-8")
    print(f"\n保存 {p.name}  {p.stat().st_size//1024}KB")

    # ---- 図（CR5 のシェルをそのまま流用）------------------------------
    sub = ("<b>加速する並進。相互作用は切ってある。</b>"
           "collision_step_exact を呼ばず、並進だけを掛けている。"
           f"1 步あたりの並進量は step(t) = {STEP0} + {ACCEL}·t ［度］で、"
           f"t=0 で {STEP0}°、t={T_STEPS-1} で {STEP0+ACCEL*(T_STEPS-1):.1f}° まで増える"
           "（A に +、B に −）。速度が増えるので系のエネルギーが増えれば "
           "<b>R′ が動くはず</b>——動かなければ、運動項が運動エネルギーを"
           "状態に入れていないことになる。"
           f"<br><br>構成: <b>A・B ともに倍音1〜17（等振幅）</b>／"
           f"初期位置 A={DEG_A}° B={DEG_B}°／T={T_STEPS} 步・1步刻み／"
           f"累積 {turns:+.4f} 周")
    note = ("読出しの r は記録しているだけで、回転には適用していない"
            "（相互作用を無効化しているため）。"
            "<br><br>"
            "○ は固定サイズで大きさは何も表さない。波形は η=0 のスライスを"
            "中心位相 0° に揃えたもの、細線は包絡。正規化は A・B・全時刻で共通。")
    html = (SHELL.read_text(encoding="utf-8")
            .replace("__TITLE__", "加速する並進 — 相互作用なし・一定加速度")
            .replace("__EYEBROW__", "CR7 · 加速")
            .replace("__SUB__", sub)
            .replace("__NOTE__", note)
            .replace("__DATA__", data))
    ph = HERE / "cr7_accelerated_v1.html"
    ph.write_text(html, encoding="utf-8")
    print(f"保存 {ph.name}  {ph.stat().st_size//1024}KB  ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
