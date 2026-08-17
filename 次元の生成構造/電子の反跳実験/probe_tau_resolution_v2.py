#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""τ 分解能 v2——衝突だけで Δθ が 44°/步 動く、その原因を切り分ける

v1 の結果
---------
  W1 読出しは連続。1/100 セル（0.00703°）の掃引に対し直線残差 8.5e-14°。
     **位置位相は量子化されていない。** 重心は第1フーリエ係数であって
     格子番号ではないので、512 セルは読出しを刻まない。
  W3 ω₀ を 1/4・1/16 にしても支配周期は 5.36 → 5.39 → 5.49 步のまま。
     **ズームできない。** 軌道も規格化時刻で一致しない（RMS差が基準の1.4倍）。
  W4 並進 v≡0（力学を全部止める）でも Δθ は全範囲 ±60° を掃き、
     1步あたり中央 44.5° 動く。

44.5°/步 には二つの説明があり、含意が正反対になる。

  仮説A「相互作用量子」: 1 步 = 1 衝突が本当に大角度の事象で、Δθ は
     実際にそれだけ動く。ならば細分は原理的に不可能で、CR3/CR4 の
     支配周期 5.36 步も衝突写像そのものの性質ということになる。
  仮説B「読出しの悪条件」: 混合で |z|（第1モーメントの大きさ）が小さくなり、
     重心角が不定になって数値が暴れているだけ。ならば Δθ の大変位は
     物理ではなく計器の破綻で、CR3/CR4 の周期も疑わしい。

検定
----
  X1 |z| の追跡: 衝突のみの走行で |z_A|,|z_B| を記録し、|dΔθ| との相関を見る。
     |z| が小さい步で大きく飛ぶなら仮説B。
  X2 衝突のみの走行のスペクトル: 支配周期が CR4 の 5.36 步と一致するか。
     一致するなら CR4 の「振動」は力学ではなく衝突写像の性質。
  X3 別計器での追試: 重心（第1モーメント）に依らない位置計器で同じことを
     測る。円周パワー分布の相互相関ピークで A・B の相対ずれを直接測る
     （|z| に依存しない）。二つの計器が一致すれば仮説A、割れれば仮説B。
  X4 単体の自由発展: B を置かず A だけを衝突させたときの θ_A の動き。
     二体相互作用ではなく衝突写像自体が位相を動かすかを見る。

使い方: python3 probe_tau_resolution_v2.py
出力  : result_tau_resolution_v2.json
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


_uni = _load("uni_tr2", UNI / "unified_interaction_v1.py")
K = _load("kin_tr2", UNI / "unified_kinetic_v1.py")
_cr0 = _load("cr0_tr2", HERE / "run_cr0_control_no_theta_v2.py")
_cr1 = _load("cr1_tr2", HERE / "run_cr1_kinetic_feedback_v1.py")

DEG_A, DEG_B = -30.0, +30.0
OMEGA0 = np.pi / 72.0
PACK_A, PACK_B = tuple(range(1, 18)), tuple(range(1, 4))


def make_ab(sp, slope, icept, tag):
    case = _uni.two_body_base.explicit_packet_case(
        mode=tag, packet_a=PACK_A, packet_b=PACK_B,
        packet_a_shift=_cr0.shift_for_deg(DEG_A, slope, icept),
        packet_b_shift=_cr0.shift_for_deg(DEG_B, slope, icept))
    a = _uni.two_body_base.make_case_state(sp, case, "A", hair_enabled=True)
    b = _uni.two_body_base.make_case_state(sp, case, "B", hair_enabled=True)
    return a, b


def xcorr_shift_deg(a, b, n_chi, n_eta):
    """|z| に依らない相対位置計器。

    円周パワー分布 P_A, P_B の巡回相互相関のピーク位置を、放物線内挿で
    セル以下まで求める。重心（第1モーメント）とは独立な計器。
    """
    PA = _cr0.power_chi(a, n_chi, n_eta)
    PB = _cr0.power_chi(b, n_chi, n_eta)
    c = np.fft.irfft(np.fft.rfft(PA) * np.conj(np.fft.rfft(PB)), n=n_chi).real
    k = int(np.argmax(c))
    y0, y1, y2 = c[(k - 1) % n_chi], c[k], c[(k + 1) % n_chi]
    den = (y0 - 2 * y1 + y2)
    frac = 0.5 * (y0 - y2) / den if den != 0 else 0.0
    pos = (k + frac) * 360.0 / n_chi
    return float(pos - 360.0 if pos > 180.0 else pos)


def spec_period(x, top=4):
    y = np.asarray(x, float); y = y - y.mean()
    Y = np.abs(np.fft.rfft(y * np.hanning(len(y)))); f = np.fft.rfftfreq(len(y))
    Y[0] = 0
    tot = float(np.sum(Y ** 2))
    idx = np.argsort(Y)[::-1]
    pk = [i for i in idx if 0 < i < len(Y) - 1 and Y[i] > Y[i - 1]
          and Y[i] > Y[i + 1]][:top]
    return [(float(1 / f[i]), float(100 * Y[i] ** 2 / tot)) for i in pk]


def main() -> None:
    t0 = time.time()
    sp = _uni.two_body_base.build_source_params(
        _uni.two_body_base.Params(high_n=63, recursive_collision_count=200))
    n_chi, n_eta = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = _cr0.calibrate_shift(sp, n_chi, n_eta)
    out = {}
    T = 1500

    # ---- X1/X2/X3 衝突のみ ------------------------------------------
    print("【X1-X3】並進 v≡0（衝突のみ）の走行 T=1500")
    a, b = make_ab(sp, slope, icept, "tr2_collonly")
    chi, zA, zB, xc = [], [], [], []
    for _ in range(T):
        pa, za = _cr0.circle_position(a, n_chi, n_eta)
        pb, zb = _cr0.circle_position(b, n_chi, n_eta)
        chi.append(float(np.degrees(np.angle(np.exp(1j * (pa - pb))))))
        zA.append(za); zB.append(zb)
        xc.append(xcorr_shift_deg(a, b, n_chi, n_eta))
        a, b, _ = _uni.collision_step_exact(a, b, sp)
    chi = np.array(chi); zA = np.array(zA); zB = np.array(zB)
    xc = np.array(xc)
    d = np.abs(np.diff(chi))
    dx = np.abs(np.diff(np.degrees(np.angle(np.exp(1j * np.radians(xc))))))

    print(f"  |z_A|: 初期 {zA[0]:.6f}  平均 {zA.mean():.6f}  "
          f"最小 {zA.min():.3e}  最大 {zA.max():.6f}")
    print(f"  |z_B|: 初期 {zB[0]:.6f}  平均 {zB.mean():.6f}  "
          f"最小 {zB.min():.3e}  最大 {zB.max():.6f}")
    zmin = np.minimum(zA, zB)[:-1]
    cc = float(np.corrcoef(np.log10(np.maximum(zmin, 1e-15)), d)[0, 1])
    print(f"  |dΔθ| と log|z|min の相関: {cc:+.4f}  "
          f"（強い負なら『|z| が小さい步で飛ぶ』＝仮説B）")
    q = np.quantile(zmin, [0.25, 0.75])
    print(f"  |z|min 下位25%の步の |dΔθ| 中央 "
          f"{np.median(d[zmin <= q[0]]):8.3f}°  /  "
          f"上位25% {np.median(d[zmin >= q[1]]):8.3f}°")
    out["X1"] = {"zA": [zA[0], float(zA.mean()), float(zA.min())],
                 "zB": [zB[0], float(zB.mean()), float(zB.min())],
                 "corr_d_logzmin": cc,
                 "d_med_lowz": float(np.median(d[zmin <= q[0]])),
                 "d_med_highz": float(np.median(d[zmin >= q[1]]))}

    print("\n  スペクトル（衝突のみ）と CR4 case17_3 の比較")
    sc = spec_period(chi)
    print("   衝突のみ 重心Δθ : " +
          "  ".join(f"{p:7.3f}步({s:5.2f}%)" for p, s in sc))
    sx = spec_period(xc)
    print("   衝突のみ 相互相関: " +
          "  ".join(f"{p:7.3f}步({s:5.2f}%)" for p, s in sx))
    p4 = HERE / "cr4_relative_case17_3_data_v1.json"
    if p4.exists():
        d4 = json.loads(p4.read_text(encoding="utf-8"))
        print("   CR4 case17_3     : " +
              "  ".join(f"{e['period']:7.3f}步({e['share']:5.2f}%)"
                        for e in d4["spec"][:4]))
        out["X2_cr4_spec"] = d4["spec"][:4]
    out["X2_collonly_spec"] = sc
    out["X2_collonly_xcorr_spec"] = sx

    print("\n  二つの計器の一致（重心 Δθ と相互相関）")
    print(f"   重心   : 範囲 [{chi.min():+8.3f},{chi.max():+8.3f}]°  "
          f"|dΔθ|中央 {np.median(d):8.3f}°")
    print(f"   相互相関: 範囲 [{xc.min():+8.3f},{xc.max():+8.3f}]°  "
          f"|dΔθ|中央 {np.median(dx):8.3f}°")
    agree = float(np.median(np.abs(np.degrees(np.angle(
        np.exp(1j * np.radians(chi - xc)))))))
    print(f"   両者の差の中央値: {agree:.3f}°  "
          f"（小さければ同じものを見ている＝仮説A）")
    out["X3"] = {"xcorr_range": [float(xc.min()), float(xc.max())],
                 "xcorr_dmed": float(np.median(dx)),
                 "centroid_dmed": float(np.median(d)),
                 "median_abs_diff": agree}

    # ---- X4 単体の自由発展 -------------------------------------------
    print("\n【X4】単体 A のみを衝突させる（相手なし＝自己衝突）")
    a1, b1 = make_ab(sp, slope, icept, "tr2_solo")
    th, z1 = [], []
    for _ in range(T):
        pa, za = _cr0.circle_position(a1, n_chi, n_eta)
        th.append(float(np.degrees(pa))); z1.append(za)
        a1, _, _ = _uni.collision_step_exact(a1, a1.copy(), sp)
    th = np.unwrap(np.radians(np.array(th)))
    th = np.degrees(th)
    dth = np.abs(np.diff(th))
    print(f"  θ_A 総回転 {(th[-1]-th[0])/360.0:+.4f} 周  "
          f"1步あたり |dθ| 中央 {np.median(dth):.4e}°  "
          f"|z_A| 平均 {np.mean(z1):.6f}")
    out["X4"] = {"turns": float((th[-1] - th[0]) / 360.0),
                 "dth_med": float(np.median(dth)),
                 "z_mean": float(np.mean(z1))}

    out["meta"] = {"experiment": "tau_resolution_probe_v2",
                   "date": time.strftime("%Y-%m-%d %H:%M:%S"),
                   "T": T, "n_chi": n_chi, "n_eta": n_eta}
    p = HERE / "result_tau_resolution_v2.json"
    p.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"\n保存 {p.name}  ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
