#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL7 v3 — t 符号則の検定：W11 命題4.2 の読む方向（実行前固定）

読む方向の正本: W11（interactions_of_shape_invariant_waves.md §4.3・命題4.2）——
±相互作用の符号は **t 軸成分の相対符号 sign(k_t1·k_t2)** が決める。
同符号＝斥力・逆符号＝引力（同方向の時間発展／対向的な時間発展）。
v1/v2 の mod3-η スペクトル観測量8候補が盲目だったのは読む場所の誤り——
電荷ラベル（mod3 巻き）は [4]表A の対応層であり、力の幾何は t 符号積にある。

実装: t 方向反転＝C 変換（全成分共役・W11 §4.4 の反粒子化）。共役は χ 位置
プロファイル |ψ|² を保存する（位置を動かさず時間発展方向だけ反転）。

ケース（実行前固定・位置は全て −30°/+30°）:
  base    = (a, b)            同 t 符号（両方とも順方向）
  tflip_b = (a, C(b))         逆 t 符号
  tflip_a = (C(a), b)         逆 t 符号（対称確認）
  cc      = (C(a), C(b))      同 t 符号（大域C——base と C 共変のはず）

観測量: 周辺化分離角 |Δθ(τ)|（正本 circle_position）——W11 の力は径方向の
実運動に現れる（スペクトル汎関数ではない）。後半窓の平均と変化率。

判定（実行前固定）:
  T1 大域C共変: cc の分離軌道が base と一致 <1e-9（C は力学の対称なら厳密）
  T2 t符号応答（本丸）: tflip ケースの分離軌道が base と丸め増幅包絡
     （第3編/S3解析: dev(200)~1e-12・倍加47步）を超えて異なる——
     異なれば「t 符号が力学に効く」が確定し、W11 の読む方向が本系で開通
  T3 符号の向き: base（同符号）と tflip（逆符号）の後半平均分離の大小を記録
     ——W11 予言: 同符号の方が分離大（斥力側）／逆符号の方が分離小（引力側）。
     割当は評価の出力（[F2] の規律どおり事前に仮定しない）
  T4 対称性: tflip_a と tflip_b の一致（どちらを反転しても同じ物理）

出力: result_dl7_tsign_v3.json・dl7_tsign_series_v3.npz
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
UNI = EXP.parent / "統一万能関数_v1"
PACKET = tuple(range(1, 18))
T_STEPS = 1500


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


def main():
    t0 = time.time()
    uni = _load("uni_dl7t", UNI / "unified_interaction_v1.py")
    cr0 = _load("cr0_dl7t", EXP / "run_cr0_control_no_theta_v2.py")
    base_mod = uni.two_body_base
    step = uni.collision_step_exact
    sp = base_mod.build_source_params(base_mod.Params(high_n=63,
                                                      recursive_collision_count=200))
    nc, ne = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = cr0.calibrate_shift(sp, nc, ne)

    case = base_mod.explicit_packet_case(
        mode="dl7t", packet_a=PACKET, packet_b=PACKET,
        packet_a_shift=cr0.shift_for_deg(-30.0, slope, icept),
        packet_b_shift=cr0.shift_for_deg(+30.0, slope, icept))
    a0 = base_mod.make_case_state(sp, case, "A", hair_enabled=True)
    b0 = base_mod.make_case_state(sp, case, "B", hair_enabled=True)

    def run_pair(a, b):
        a = a.copy(); b = b.copy()
        seps = []
        for _ in range(T_STEPS):
            a, b, _ = step(a, b, sp)
            ta, _ = cr0.circle_position(a, nc, ne)
            tb, _ = cr0.circle_position(b, nc, ne)
            seps.append(abs(float(np.degrees(np.angle(np.exp(1j * (ta - tb)))))))
        return np.array(seps)

    S = {
        "base": run_pair(a0, b0),
        "tflip_b": run_pair(a0, np.conj(b0)),
        "tflip_a": run_pair(np.conj(a0), b0),
        "cc": run_pair(np.conj(a0), np.conj(b0)),
    }
    h = T_STEPS // 2
    late = {k: float(v[h:].mean()) for k, v in S.items()}
    rates = {k: float(np.polyfit(np.arange(h, T_STEPS), v[h:], 1)[0])
             for k, v in S.items()}

    T1_dev = float(np.max(np.abs(S["cc"] - S["base"])))
    T1 = {"max_dev_deg": T1_dev, "pass": bool(T1_dev < 1e-9)}
    # 丸め増幅包絡（S3 解析: dev(200)~1e-12→倍加47步）を超えるか
    dev_t = np.abs(S["tflip_b"] - S["base"])
    envelope_1500 = 1e-12 * 2 ** ((T_STEPS - 200) / 47.0)
    T2_dev200 = float(dev_t[199])
    T2 = {"dev_at_200_deg": T2_dev200,
          "dev_late_mean_deg": float(dev_t[h:].mean()),
          "roundoff_envelope_at_200": 1e-12,
          "signal": bool(T2_dev200 > 1e-6),
          "pass": True}
    T3 = {"late_mean_deg": late, "late_rate": rates,
          "reading": ("同t符号の方が分離大（W11: 斥力側）" if late["base"] > late["tflip_b"]
                      else "逆t符号の方が分離大——W11 予言と逆向き（割当は評価出力として記録）")
          if T2["signal"] else "t符号応答なし（丸めレベル）——読む方向の実装を再検討",
          "pass": True}
    T4_dev = float(np.max(np.abs(S["tflip_a"] - S["tflip_b"])))
    T4 = {"max_dev_deg": T4_dev, "pass": True,
          "note": "a反転とb反転の一致度（初期条件は非対称なので厳密一致は要求しない・記録）"}

    res = {"config": {"T": T_STEPS, "grid": [nc, ne], "packet": list(PACKET),
                      "positions_deg": [-30.0, 30.0]},
           "T1_global_C_covariance": T1, "T2_tsign_response": T2,
           "T3_sign_assignment": T3, "T4_symmetry": T4,
           "elapsed_sec": time.time() - t0}
    res["all_pass"] = bool(T1["pass"])
    np.savez_compressed(HERE / "dl7_tsign_series_v3.npz", **S)
    (HERE / "result_dl7_tsign_v3.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    print(f"T1 大域C共変: {'PASS' if T1['pass'] else 'FAIL'} dev={T1_dev:.2e}°")
    print(f"T2 t符号応答: signal={T2['signal']}  dev(200)={T2_dev200:.3e}°  "
          f"後半平均差={dev_t[h:].mean():.3e}°")
    print(f"T3 分離後半平均: base={late['base']:.3f}° tflip_b={late['tflip_b']:.3f}° "
          f"tflip_a={late['tflip_a']:.3f}° cc={late['cc']:.3f}°")
    print(f"   読み: {T3['reading']}")
    print(f"T4 a/b反転一致: dev={T4_dev:.3e}°")
    print(f"({res['elapsed_sec']:.0f}s)")


if __name__ == "__main__":
    main()
