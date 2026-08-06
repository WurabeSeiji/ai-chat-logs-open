#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""粒子周期表・番地走査 v1: N=4..144 の分配観測量を収穫し U^n=I 番地格子に写す

目標（木原 2026-08-06）: 粒子の周期表。周期律の担い手は m/n 番地という見込み。
  - 各粒子種（ボゾン/フェルミオン/バリオン）＝番地（族）。分配値=種の指紋（反射率）
  - N は実質エネルギー。低N=低エネルギーでも第一周期の特性は現れるはず
  - 時計 ω≈π/72/step → 一周=144ステップ。N=144 まで走査しないと周期の閉じ
    （対称性）が見えない疑い。144の約数は結晶学的メニュー{3,4,6}を全て含む
  - 質量階層性の痕跡: 質量度と番地の相関

方法: 各Nで T=4000 走行（準安定窓[2000,4000]・400サンプル）。ヤコビアンは大Nで
不可能なため力学のみの観測量に限定:
  (1) 集団時計 ω_clock と π/72 比（時計番地）
  (2) 分配値5種（反射率型）:
      perp_ratio  = 親平面外パワー/全パワー
      axis_share  = sv3²/(sv1²+sv2²+sv3²)（占有3方向中の軸の取り分）
      plane_bal   = sv2²/sv1²（平面2方向の均衡）
      corr_zx/zy  = SVD方向復調読出しの軸-平面相関
  (3) 質量度 = 4detΓ/(trΓ)²（平面読出し対のGram・無次元非コヒーレンス）
  (4) 有意方向数 n_sig（>5%特異値）= 見かけ次元
  番地写像: 各分配値 v → θ/π=arcsin(√v)/π の連分数近似。分母階層 n≤16 と
  特別族 {5,10,25,31,62,72,124,144,248} での最良番地と偏差を記録。

読み方（事前登録・記述的）:
  P1 周期律: 同一番地が N を跨いで再帰するか（プラトー列）
  P2 種の閾値: 新しい番地が特定の N から出現するか
  P3 質量階層: 質量度が番地分母と相関するか
  P4 時計の閉じ: N→144 で時計番地に通約構造が現れるか

使い方: python3 run_periodic_address_scan_v1.py
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
SPACE = HERE.parent / "空間軸3軸と固有時間の創生_v1"
spec1 = importlib.util.spec_from_file_location(
    "pre1_pts", SPACE / "run_pre_2plus1_structure_v1.py")
pre1 = importlib.util.module_from_spec(spec1)
sys.modules[spec1.name] = pre1
spec1.loader.exec_module(pre1)
abl = pre1.abl

T_END = 4000
WIN = (2000, 4000)
SAMPLE_EVERY = 5
N_LIST = (4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
          20, 24, 28, 32, 36, 40, 48, 56, 64, 72, 80, 96, 112, 128, 144)
SPECIAL_DENOMS = (5, 10, 25, 31, 62, 72, 124, 144, 248)


def address_map(v: float) -> dict:
    """分配値 v∈(0,1) → 番地候補。小分母(≤16)最良と特別族最良を返す。"""
    v = min(max(v, 1e-12), 1 - 1e-12)
    x = float(np.arcsin(np.sqrt(v)) / np.pi)          # θ/π ∈ (0, 1/2)
    fr_small = Fraction(x).limit_denominator(16)
    dev_small = v - float(np.sin(np.pi * fr_small.numerator / fr_small.denominator) ** 2)
    best_sp = None
    for nd in SPECIAL_DENOMS:
        m = int(round(x * nd))
        if m < 1:
            continue
        pred = float(np.sin(np.pi * m / nd) ** 2)
        d = v - pred
        if best_sp is None or abs(d) < abs(best_sp[2]):
            best_sp = (m, nd, d)
    return {"theta_over_pi": x,
            "small": [fr_small.numerator, fr_small.denominator, dev_small],
            "special": list(best_sp) if best_sp else None}


def scan_one(n: int) -> dict:
    t0 = time.time()
    sys_lr, v, _, _, _, p, q, Z, wp = abl.build_init(n, True)
    M = sys_lr.m
    samples = []
    for t in range(T_END):
        Z, wp = abl.evolve(sys_lr, Z, wp)
        if WIN[0] <= t < WIN[1] and (t % SAMPLE_EVERY == 0):
            samples.append(Z.copy())
    S = np.array(samples)
    ns = S.shape[0]

    def perp(Zc):
        return Zc - p * (p @ Zc) - q * (q @ Zc)
    Sp = S - np.outer(S @ p, p) - np.outer(S @ q, q)
    pow_perp = float(np.mean(np.sum(np.abs(Sp) ** 2, axis=1)))
    pow_tot = float(np.mean(np.sum(np.abs(S) ** 2, axis=1)))
    perp_ratio = pow_perp / pow_tot

    X = np.hstack([Sp.real, Sp.imag])
    Xc = X - X.mean(axis=0)
    _, sv, Vt = np.linalg.svd(Xc, full_matrices=False)
    sv_rel = sv / sv[0]
    n_sig = int(np.sum(sv_rel > 0.05))
    p3 = sv[:3] ** 2
    axis_share = float(p3[2] / p3.sum())
    plane_bal = float(p3[1] / p3[0])

    # SVD上位3方向を実体方向とみなす復調読出し（大Nでも可能な代理）
    d1, d2, d3 = Vt[0], Vt[1], Vt[2]
    cp = S @ p
    cq = S @ q
    phi = np.unwrap(np.angle(cp + 1j * cq))
    omega_clock = float(np.polyfit(np.arange(ns), phi, 1)[0])

    def demod(d):
        dc = d[:M] + 1j * d[M:]
        dc = dc / np.linalg.norm(dc)
        c = Sp @ np.conj(dc)
        F = np.fft.fft(c - c.mean())
        pk = int(np.argmax(np.abs(F)))
        f = pk / ns if pk <= ns // 2 else (pk - ns) / ns
        n_k = round(2 * np.pi * f / omega_clock) if omega_clock != 0 else 0
        return c * np.exp(-1j * n_k * phi)

    r1, r2, r3 = demod(d1), demod(d2), demod(d3)
    def corr(a, b):
        return float(abs(np.mean(a * np.conj(b)))
                     / max(np.sqrt(np.mean(np.abs(a) ** 2) * np.mean(np.abs(b) ** 2)), 1e-300))
    corr_zx, corr_zy, corr_xy = corr(r3, r1), corr(r3, r2), corr(r1, r2)

    # 質量度: 平面読出し対のGram非コヒーレンス（無次元 4detΓ/trΓ²）
    G11 = np.mean(np.abs(r1) ** 2)
    G22 = np.mean(np.abs(r2) ** 2)
    G12 = np.mean(r1 * np.conj(r2))
    detG = G11 * G22 - abs(G12) ** 2
    mass_deg = float(4 * detG / (G11 + G22) ** 2)

    parts = {"perp_ratio": perp_ratio, "axis_share": axis_share,
             "plane_bal": plane_bal, "corr_zx": corr_zx, "corr_zy": corr_zy}
    addr = {k: address_map(vv) for k, vv in parts.items()}

    return {"N": n, "M": M, "omega_clock": omega_clock,
            "clock_over_pi72": float(omega_clock / (np.pi / 72)),
            "n_sig": n_sig, "sv_rel_top6": [float(x) for x in sv_rel[:6]],
            "partitions": parts, "addresses": addr,
            "corr_xy": corr_xy, "mass_deg": mass_deg,
            "runtime_sec": time.time() - t0}


def main() -> None:
    t0 = time.time()
    out = {"T_END": T_END, "WIN": list(WIN), "SAMPLE_EVERY": SAMPLE_EVERY,
           "N_LIST": list(N_LIST), "SPECIAL_DENOMS": list(SPECIAL_DENOMS),
           "scan": []}
    print(f"{'N':>4} {'ω/(π/72)':>9} {'次元':>4} {'perp':>7} {'axis':>7} "
          f"{'bal':>7} {'zx':>7} {'zy':>7} {'質量度':>8}  時間")
    for n in N_LIST:
        r = scan_one(n)
        out["scan"].append(r)
        pt = r["partitions"]
        print(f"{n:>4} {r['clock_over_pi72']:>9.4f} {r['n_sig']:>4} "
              f"{pt['perp_ratio']:>7.4f} {pt['axis_share']:>7.4f} "
              f"{pt['plane_bal']:>7.4f} {pt['corr_zx']:>7.4f} {pt['corr_zy']:>7.4f} "
              f"{r['mass_deg']:>8.5f}  {r['runtime_sec']:.0f}s")
        out["runtime_sec"] = time.time() - t0
        (HERE / "periodic_address_scan_result_v1.json").write_text(
            json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s → periodic_address_scan_result_v1.json")


if __name__ == "__main__":
    main()
