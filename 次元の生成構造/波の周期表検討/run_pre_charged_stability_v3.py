#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""予備実験 v3: 帯電種の長時間安定性——周期表の第二元素の直接探索

背景: v2 で低N（N≤16）の海の安定種は 1/1 無質量基底種のみと判明。
質量・電荷を持つ側帯は全て過渡（共鳴）だった。安定な第二種の本命経路は
毛（巻き数）付き帯電シード——空間時間論文 §9.5 の帯電構造（種+相棒）を
長時間走らせ、安定性を直接測る。

系: 二体万能非弾性写像（閉形式厳密解・分布実験と同一構成）。
帯電状態: 単一巻きポンプ (k=30,32,34) + 単一巻き種 (k=21, 0.2S)
（run_distribution_readouts_v1.py 時期3と同一。決定論・シードなし）。

測定（J=4000衝突、窓40×100窓）:
  (1) 電荷持続: フェルミオン的マスク上の巻き数 q=+1 の重みの時間発展。
  (2) 種モードの回転数 ρ: 帯電スペクトル最大モードの位相前進/衝突
     （全系列アンラップ、δρ≈1/J）。部分区間4分割でドリフト測定。
  (3) 相棒 (q=+3相当帯) の重み持続（census型対構造の寿命）。

判定（事前固定）:
  H_second（第二元素）: (a) q=+1重み 最終/初期 ≥ 0.5 かつ
    (b) ρ の四分区間ドリフト < 2δρ_sub かつ最良有理番地（分母≤6優先、
    次いで特別族{31,62,124,248}）に |ρ−m/n|<δρ → 安定な帯電種と判定。
  H_resonance: 減衰する場合、q重みの指数フィットで寿命 τ を、
    ρドリフトで幅を測る（共鳴のカタログ値として記録）。
  参照予言: C論文のタング実測では分母3の幅が支配的（幅比>461）——
  ロックするなら分母3系（Z₄/Z₆結晶学的メニュー）が最有力。

使い方: python3 run_pre_charged_stability_v3.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_cs3", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
v1, toy, base = ex.v1, ex.toy, ex.base

S = 8.0
J_TOT = 4000
J_WIN = 40


def main() -> None:
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    shape = (n, ne)

    def single_winding(v):
        f = np.fft.fft(v.reshape(shape), axis=0, norm="ortho")
        f[0, :] = 0.0
        f[n // 2:, :] = 0.0
        return np.fft.ifft(f, axis=0, norm="ortho").reshape(v.shape)

    a = single_winding(v1.make_bundle(sp, (30, 32, 34), "A", scale=1.0)) * S
    b = single_winding(v1.make_bundle(sp, (30, 32, 34), "B", scale=1.0)) * S
    a = a + single_winding(v1.make_bundle(sp, (21,), "A", scale=1.0)) * (0.2 * S)

    ks = np.arange(n); kk = np.where(ks <= n // 2, ks, ks - n)
    ferm_k = (np.abs(kk) % 2 == 0) & (np.abs(kk) >= 4)
    ms = np.arange(ne); mm = np.where(ms <= ne // 2, ms, ms - ne)

    # 追跡モード: 初期帯電スペクトルの最大占有 (k,m)（m=+1帯から選ぶ）
    fa0 = np.fft.fft(np.fft.fft(a.reshape(shape), axis=0, norm="ortho"),
                     axis=1, norm="ortho")
    P0 = np.abs(fa0) ** 2
    mask_m1 = (mm[None, :] == 1)
    idx = np.unravel_index(np.argmax(P0 * mask_m1), P0.shape)
    k_tr, m_tr = int(idx[0]), int(idx[1])
    print(f"追跡モード: (k={kk[k_tr]}, m=+1)  初期占有={P0[idx]:.4f}")

    phases = np.zeros(J_TOT)
    windows = []
    q_series = {"+1": [], "+3": [], "0": [], "-1": []}
    Pwin = np.zeros(shape)
    for j in range(J_TOT):
        a, b, _ = ex.collision_step_exact(a, b, sp)
        fa = np.fft.fft(np.fft.fft(a.reshape(shape), axis=0, norm="ortho"),
                        axis=1, norm="ortho")
        fb = np.fft.fft(np.fft.fft(b.reshape(shape), axis=0, norm="ortho"),
                        axis=1, norm="ortho")
        phases[j] = float(np.angle(fa[k_tr, m_tr]))
        Pwin += (np.abs(fa) ** 2 + np.abs(fb) ** 2) / 2
        if (j + 1) % J_WIN == 0:
            P = Pwin / J_WIN
            occ = P > (P.max() * 1e-8)
            F = ferm_k[:, None] & occ
            for tag, m_ in (("+1", 1), ("+3", 3), ("0", 0), ("-1", -1)):
                q_series[tag].append(float(np.sum(P[F & (mm[None, :] == m_)])))
            windows.append(j + 1)
            Pwin = np.zeros(shape)

    # 回転数（全系列アンラップ）
    ph = np.unwrap(phases)
    rho_full = float((ph[-1] - ph[0]) / (2 * np.pi * (J_TOT - 1)))
    q4 = J_TOT // 4
    rhos = [float((ph[(k + 1) * q4 - 1] - ph[k * q4]) / (2 * np.pi * (q4 - 1)))
            for k in range(4)]
    drift = float(max(rhos) - min(rhos))
    drho = 1.0 / J_TOT
    drho_sub = 1.0 / q4
    stationary = bool(drift < 2 * drho_sub)

    rho_frac = rho_full % 1.0
    fr6 = Fraction(rho_frac).limit_denominator(6)
    dev6 = rho_frac - float(fr6)
    best_sp = None
    for nd in (31, 62, 124, 248):
        m_ = int(round(rho_frac * nd))
        d_ = rho_frac - m_ / nd
        if best_sp is None or abs(d_) < abs(best_sp[2]):
            best_sp = (m_, nd, d_)

    q1 = np.array(q_series["+1"])
    ratio_q1 = float(q1[-1] / max(q1[0], 1e-300))
    # 指数フィット（寿命）
    pos = q1 > 0
    tau = None
    if pos.all() and q1[0] > 0:
        w = np.array(windows, float)
        coef = np.polyfit(w, np.log(q1), 1)
        tau = float(-1.0 / coef[0]) if coef[0] < 0 else float("inf")

    h_second = bool(ratio_q1 >= 0.5 and stationary and
                    (abs(dev6) < drho or (best_sp and abs(best_sp[2]) < drho)))

    print(f"\n電荷 q=+1 重み: 初期窓={q1[0]:.4f} → 最終窓={q1[-1]:.4f} 比={ratio_q1:.3f}")
    print(f"  寿命τ（指数フィット, 衝突数）= {tau if tau is not None else '正値のため∞/未定義'}")
    print(f"相棒 q=+3: 初期={q_series['+3'][0]:.4e} → 最終={q_series['+3'][-1]:.4e}")
    print(f"回転数 ρ = {rho_full:.6f}（小数部 {rho_frac:.6f} ± {drho:.6f}）")
    print(f"  四分区間 ρ = {[round(r,5) for r in rhos]}  ドリフト={drift:.6f} "
          f"（2δρ_sub={2*drho_sub:.6f}）→ {'定常' if stationary else '漂'}")
    print(f"  番地: 分母≤6 → {fr6}（偏差{dev6:+.5f}）  "
          f"特別族 → {best_sp[0]}/{best_sp[1]}（偏差{best_sp[2]:+.5f}）")
    print(f"\n判定 H_second（安定な第二元素）= {h_second}")

    out = {"S": S, "J_TOT": J_TOT, "J_WIN": J_WIN,
           "track_mode": {"k": int(kk[k_tr]), "m": 1},
           "windows": windows, "q_series": q_series,
           "ratio_q1_final_initial": ratio_q1, "tau_collisions": tau,
           "rho_full": rho_full, "rho_frac": rho_frac, "drho": drho,
           "rho_quarters": rhos, "drift": drift, "stationary": stationary,
           "address_small6": [fr6.numerator, fr6.denominator, float(dev6)],
           "address_special": list(best_sp),
           "H_second": h_second, "runtime_sec": time.time() - t0}
    (HERE / "pre_charged_stability_result_v3.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s → pre_charged_stability_result_v3.json")


if __name__ == "__main__":
    main()
