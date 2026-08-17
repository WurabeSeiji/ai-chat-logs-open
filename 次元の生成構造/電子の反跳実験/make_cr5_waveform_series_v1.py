#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CR5: 合成波形を記録する（円の大きさを廃止し、波形そのものを出す）

CR4 までの図の問題（実測で確認）
--------------------------------
図の円の半径は A の参加率 PR ひとつから作られていた:

    rad = 8 + 26*(PR[i]-PR_MIN)/(PR_MAX-PR_MIN)

  (1) この rad が A の円と B の円の**両方**に使われていた。B 用の半径は
      存在せず、記録データにも B の PR が無い。B の円は A の値を B の位置に
      描いたもので、何も表していなかった。
  (2) PR_MIN/PR_MAX はその走行の実測範囲なので、円の大小はケース内の相対値。
      ケース間で比較できない。
  (3) PR はパワー由来。さらに像が 180° 間隔で割れる構成では二つの像を
      合わせて数えるので、広がりとして読めない
      （probe_harmonic_composition_v1 T4 で確認）。

本版の方針
----------
  * 円は**固定サイズ**にする。大きさで量を表さない。
  * 代わりに **A・B の合成波形そのもの**を記録して描く。スカラーに潰さない
    ので「幅をどう定義するか」の選択が要らず、像が割れていればそれが見える。

合成波形の作り方（両方とも記録する）
------------------------------------
  wave : Σ_η ψ(χ,η) の実部。倍音の合成そのもの（符号つき）。
  env  : sqrt(Σ_η |ψ(χ,η)|²)。振幅の包絡（非負）。

  それぞれについて **中心位相を 0° に揃える**:
    (a) 位置の中心 = 円周第1モーメントの偏角（Δθ を測るのと同じ計器）を
        χ=0 に回す。A と B を重ねて形を比べられるようにするため。
    (b) 中心での搬送波位相を割り戻す。これをしないと波形が毎步回って
        形が読めない。割り戻す量は記録する（carrier_a / carrier_b）。

  正規化は **A・B・全時刻で共通の一つの係数**（包絡の最大値）で行う。
  片方だけで正規化すると A↔B のパワー移動が見えなくなるため。

出力: cr5_waveform_<tag>_data_v1.json
図  : build_cr5_pages_v1.py が cr_waveform_shell_v1.html と組み合わせる

使い方:
  python3 make_cr5_waveform_series_v1.py            # 全ケース
  python3 make_cr5_waveform_series_v1.py case3_3    # 個別
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


_uni = _load("uni_cr5", UNI / "unified_interaction_v1.py")
K = _load("kin_cr5", UNI / "unified_kinetic_v1.py")
_cr0 = _load("cr0_cr5", HERE / "run_cr0_control_no_theta_v2.py")
_cr1 = _load("cr1_cr5", HERE / "run_cr1_kinetic_feedback_v1.py")

CASES = {
    "case17_3":  {"a": tuple(range(1, 18)), "b": tuple(range(1, 4)),
                  "label": "A:倍音1〜17 / B:倍音1〜3", "note": "元の実験条件（非対称）"},
    "case3_3":   {"a": tuple(range(1, 4)),  "b": tuple(range(1, 4)),
                  "label": "A・B ともに倍音1〜3", "note": "対称・少倍音"},
    "case17_17": {"a": tuple(range(1, 18)), "b": tuple(range(1, 18)),
                  "label": "A・B ともに倍音1〜17", "note": "対称・多倍音"},
    "even":      {"a": (1, 2, 3, 4, 6, 8, 10, 12, 14, 16),
                  "b": (1, 2, 3, 4, 6, 8, 10, 12, 14, 16),
                  "label": "A・B ともに倍音1,2,3,4,6,8,10,12,14,16",
                  "note": "偶数優位（r=1/20・13.89 步/周期・像が180°で2つ）"},
    "even48":    {"a": tuple([1, 2, 3, 4] + list(range(6, 49, 2))),
                  "b": tuple([1, 2, 3, 4] + list(range(6, 49, 2))),
                  "label": "A・B ともに倍音1,2,3,4,6,8,…,48（26本）",
                  "note": "偶数優位を帯域48まで拡大（r=1/52・22.58 步/周期）"},
    "even96":    {"a": tuple([1, 2, 3, 4] + list(range(6, 97, 2))),
                  "b": tuple([1, 2, 3, 4] + list(range(6, 97, 2))),
                  "label": "A・B ともに倍音1,2,3,4,6,8,…,96（50本）",
                  "note": "偶数優位を帯域96まで拡大（r=1/100・31.36 步/周期）"},
}
DEG_A, DEG_B = -30.0, +30.0
T_STEPS = 400          # 1步刻み。支配周期は 4〜14 步なので 30〜90 周期ぶん。
OMEGA0 = np.pi / 72.0
CHI_KEEP = 256         # 記録する χ 点数（512 の 1/2 = 1.41°刻み）


def dominant_winding(psi, n_chi, n_eta):
    """η 巻き m のうちパワー最大のビン番号（FFT の生添字）を返す。"""
    f = np.fft.fft(psi.reshape(n_chi, n_eta), axis=1)
    return int(np.argmax(np.sum(np.abs(f) ** 2, axis=0)))


def centered_profiles(psi, n_chi, n_eta, keep, mbin=None):
    """中心位相を 0° に揃えた (合成波形, 包絡, 中心角[rad], 搬送波位相, m1占有率)。

    合成波形の取り方（三つ試して二つを棄却した）
      × Σ_η ψ の実部: η 方向の単純和は m=0 成分だけを拾う。この系の状態は
        m=+1 と m=+2 にしか無いので **恒等的に 0**（実測 |Σ_η ψ|max = 0.000000）。
      × 支配巻き m ひとつへの射影: m=+1 と m=+2 の間を毎步往復するため、
        固定 m では占有率が 1.0000〜0.0033 まで落ちて波形が消える。
        （実測: 最大占有の m は +1 と +2 を行き来し、両者の和は常に 1.00000。
          拡散ではなく二モード交換。）
      ○ **η=0 のスライス** ψ(χ,0): 射影も選択も入らない実際の場の値。
        全ての巻き成分をそのまま含む。

    包絡 env = sqrt(Σ_η |ψ|²) は η に依存しない量なので、そのまま併記する。
    """
    A2 = psi.reshape(n_chi, n_eta)
    env = np.sqrt(np.sum(np.abs(A2) ** 2, axis=1))     # 包絡（η 非依存）
    f = np.fft.fft(A2, axis=1)
    Pm = np.sum(np.abs(f) ** 2, axis=0)
    frac = float(Pm[1] / Pm.sum())                     # m=+1 の占有率（記録用）
    wave_c = A2[:, 0].copy()                           # η=0 スライス（実際の場）

    # (a) 位置の中心（Δθ と同じ計器＝円周第1モーメント）を χ=0 へ
    ang, _ = _cr0.circle_position(psi, n_chi, n_eta)
    shift = int(np.rint(ang / (2.0 * np.pi) * n_chi))
    env = np.roll(env, -shift)
    wave_c = np.roll(wave_c, -shift)

    # (b) 中心での搬送波位相を割り戻す（割り戻さないと毎步回って形が読めない）
    carrier = float(np.angle(wave_c[0])) if abs(wave_c[0]) > 0 else 0.0
    wave_c = wave_c * np.exp(-1j * carrier)


    # 中心を配列の真ん中へ置き、χ を −180..+180 で見せる
    env = np.roll(env, n_chi // 2)
    wave_c = np.roll(wave_c, n_chi // 2)

    step = n_chi // keep
    return (wave_c.real[::step], env[::step], float(ang), carrier, frac)


def run_case(tag):
    cfg = CASES[tag]
    sp = _uni.two_body_base.build_source_params(
        _uni.two_body_base.Params(high_n=63, recursive_collision_count=200))
    n_chi, n_eta = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = _cr0.calibrate_shift(sp, n_chi, n_eta)
    case = _uni.two_body_base.explicit_packet_case(
        mode="cr5_" + tag, packet_a=cfg["a"], packet_b=cfg["b"],
        packet_a_shift=_cr0.shift_for_deg(DEG_A, slope, icept),
        packet_b_shift=_cr0.shift_for_deg(DEG_B, slope, icept))
    a = _uni.two_body_base.make_case_state(sp, case, "A", hair_enabled=True)
    b = _uni.two_body_base.make_case_state(sp, case, "B", hair_enabled=True)

    theta0 = float(_cr1.toy.theta_from_ab(a, b, sp).theta)
    omega, v = OMEGA0, 0.0
    H = {k: [] for k in ("chi", "mid", "r", "v", "omega", "cloRe", "cloIm",
                         "PA", "PB", "za", "zb", "ca", "cb", "fa", "fb")}
    WA, WB, EA, EB = [], [], [], []

    for _ in range(T_STEPS):
        pa, za = _cr0.circle_position(a, n_chi, n_eta)
        pb, zb = _cr0.circle_position(b, n_chi, n_eta)
        chi = float(np.angle(np.exp(1j * (pa - pb))))
        mid = float(np.angle(np.exp(1j * pa) + np.exp(1j * pb)))
        r_now = float(_cr1.toy.theta_from_ab(a, b, sp).reflection_rate)
        acc = -4.0 * np.sin(omega / 2.0) ** 2 * chi
        v += acc
        omega += (1.0 - r_now) * acc          # κ = 1−r（CR1 で決めた形）

        wa, ea, _, ca, fa = centered_profiles(a, n_chi, n_eta, CHI_KEEP)
        wb, eb, _, cb, fb = centered_profiles(b, n_chi, n_eta, CHI_KEEP)
        WA.append(wa); EA.append(ea); WB.append(wb); EB.append(eb)

        Z = complex(np.sum(a * a) + np.sum(b * b))
        H["chi"].append(round(chi, 8)); H["mid"].append(round(mid, 6))
        H["r"].append(round(r_now, 6)); H["v"].append(float(f"{v:.6g}"))
        H["omega"].append(round(omega, 8))
        H["cloRe"].append(float(f"{Z.real:.4g}"))
        H["cloIm"].append(float(f"{Z.imag:.4g}"))
        H["PA"].append(round(float(np.vdot(a, a).real), 8))
        H["PB"].append(round(float(np.vdot(b, b).real), 8))
        H["za"].append(round(za, 6)); H["zb"].append(round(zb, 6))
        H["ca"].append(round(ca, 5)); H["cb"].append(round(cb, 5))
        H["fa"].append(round(fa, 6)); H["fb"].append(round(fb, 6))

        a = K.k_translate_flat(a, -v, n_chi, n_eta)
        a, b, _ = _uni.collision_step_exact(a, b, sp)

    WA = np.array(WA); WB = np.array(WB); EA = np.array(EA); EB = np.array(EB)
    # A・B・全時刻で共通の一つの係数で正規化（A↔B のパワー移動を見えるようにする）
    scale = float(max(EA.max(), EB.max()))            # 包絡用（A・B・全時刻で共通）
    wscale = float(max(np.abs(WA).max(), np.abs(WB).max()))  # 波形用（同上）
    rnd = lambda X, s_: [[round(float(x), 4) for x in row] for row in (X / s_)]

    chi = np.array(H["chi"])
    out = {
        "experiment": "cr5_waveform_series_v1", "case": tag,
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "config": {"packet_a": list(cfg["a"]), "packet_b": list(cfg["b"]),
                   "label": cfg["label"], "note": cfg["note"],
                   "deg_a": DEG_A, "deg_b": DEG_B, "T": T_STEPS, "every": 1,
                   "omega0": OMEGA0, "kappa": "1-r",
                   "chi_keep": CHI_KEEP, "n_chi": n_chi, "n_eta": n_eta,
                   "theta0": theta0, "r0": H["r"][0],
                   "period_pred": float(np.pi / theta0) if theta0 > 1e-14 else None,
                   "norm_scale": scale, "wave_scale": wscale},
        "chi": H["chi"], "dchi": [float(f"{x:.6g}") for x in np.gradient(chi)],
        "mid": H["mid"], "r": H["r"], "v": H["v"], "omega": H["omega"],
        "cloRe": H["cloRe"], "cloIm": H["cloIm"],
        "PA": H["PA"], "PB": H["PB"], "za": H["za"], "zb": H["zb"],
        "carrierA": H["ca"], "carrierB": H["cb"],
        "fracA": H["fa"], "fracB": H["fb"],
        "waveA": rnd(WA, wscale), "waveB": rnd(WB, wscale),
        "envA": rnd(EA, scale), "envB": rnd(EB, scale),
    }
    p = HERE / f"cr5_waveform_{tag}_data_v1.json"
    p.write_text(json.dumps(out, ensure_ascii=False, separators=(",", ":")),
                 encoding="utf-8")
    return out, p


def main() -> None:
    tags = sys.argv[1:] or list(CASES)
    rows = []
    for tag in tags:
        t0 = time.time()
        print(f"走行中 {tag} ({CASES[tag]['label']}) …", end="", flush=True)
        out, p = run_case(tag)
        rows.append((tag, out, p))
        print(f" 完了 {time.time()-t0:.1f}s  {p.stat().st_size//1024}KB")

    print()
    print(f"{'ケース':11} {'r':>10} {'θ[°]':>8} {'予測周期':>9} "
          f"{'P_A 初→終':>22} {'|z_A| 初→終':>20}")
    print("-" * 88)
    for tag, out, _ in rows:
        c = out["config"]
        print(f"{tag:11} {c['r0']:10.6f} {np.degrees(c['theta0']):8.3f} "
              f"{c['period_pred']:9.4f} "
              f"{out['PA'][0]:10.6f}→{out['PA'][-1]:10.6f} "
              f"{out['za'][0]:9.6f}→{out['za'][-1]:9.6f} "
              f"m=+1占有 {min(out['fracA']):.4f}〜{max(out['fracA']):.4f}")
    print("\n※円の大きさは廃止した。波形は中心位相を 0° に揃えて記録している。"
          "\n  正規化は A・B・全時刻で共通なので、波形の高さの差はパワーの差そのもの。")


if __name__ == "__main__":
    main()
