#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
コインと井戸型コップの交差項の図を、波から計算した R で描き直す (2026-11-25)。

fig_cup_well_cross.py の図では、円の半径 R を手で置いていた。
ここでは、円の半径を波の倍音の係数だけから計算する (パーセバルの等式)。
出どころ: 「閉じた弦の複素パルスにおける R の検討」(2026-10-29) §6.2、
          「同じ弦の上の二つの複素パルスにおける R の軌道」(2026-10-29) 要旨 2。
          どちらもその日付の時点の仮説。

  コインの R    : R1^2 = sum_{m>=1} A_m^2       (直流は含めない)
  コップの R    : R2^2 = sum_{m>=1} B_m^2
  和の R        : R^2(Δ) = R1^2 + R2^2 + 2 sum_m A_m B_m cos(m Δ)、Δ = 中心位相の差
  円の半径には和の R(Δ) を使う。コップを動かすと半径も変わる。

波形・係数・交差項は fig_cup_well_cross.py と同じ数のまま、何も変えていない。
変えたのは円の半径だけである。縮尺も前の図と同じ (単位の長さ 1 = 200 ピクセル) なので、
円は前の図の約 1/10 の大きさになる。

描き方
  コイン・コップ : 中心からの距離 r = R(Δ) + 値
  交差項         : r = R(Δ) + c。c は振幅の単位の二乗なので、単位の長さ 1 で一度割って長さにしている
                   (前の図の「R で一度割る」は、手で置いた R = 単位の長さ 1 で割ることだった)。

出力
  figures/well_cross_derived_R/fig_cup_well_cross_derivedR_phaseXXX.svg   15 度刻み、24 枚
  results/fig_cup_well_cross_derived_R.txt

実行: python3 fig_cup_well_cross_derived_R.py
"""
import math
import os

import fig_cup_coin_table as base
import fig_cup_well_cross as wc
from fig_cup_coin_table import (BLACK, RED, GREEN, BLUE, GREY, N_COIN, THETA_COIN,
                                f_coin, legend, save, svg_open, text)

OUT_DIR = os.path.join(base.FIG_DIR, "well_cross_derived_R")
RES_PATH = os.path.join(base.RES_DIR, "fig_cup_well_cross_derived_R.txt")

PX = 200.0          # 単位の長さ 1 あたりのピクセル。前の図の R_PX と同じ
SAMPLES = 3600

R1 = math.sqrt(sum(a * a for a in wc.A))
R2 = math.sqrt(sum(b * b for b in wc.B))


def cross_sum(theta):
    """2 sum_m A_m B_m cos(m Δ)。和の R^2 に入る交差項。"""
    d = math.radians(theta - THETA_COIN)
    return 2.0 * sum(wc.A[m] * wc.B[m] * math.cos(m * d) for m in range(1, N_COIN + 1))


def r_total(theta):
    return math.sqrt(R1 * R1 + R2 * R2 + cross_sum(theta))


def phase_figure(path, theta):
    w, h = 900, 760
    cx, cy = 450.0, 400.0
    rad = r_total(theta)

    def pt(r, deg):
        a = math.radians(deg)
        return cx + r * PX * math.cos(a), cy - r * PX * math.sin(a)

    def curve(func, color, width, dash=None):
        pts = []
        for i in range(SAMPLES + 1):
            deg = 360.0 * i / SAMPLES
            x, y = pt(rad + func(deg), deg)
            pts.append(f"{x:.2f},{y:.2f}")
        d = f' stroke-dasharray="{dash}"' if dash else ""
        return (f'<polyline points="{" ".join(pts)}" fill="none" stroke="{color}" '
                f'stroke-width="{width}" stroke-linejoin="round"{d}/>')

    s = svg_open(w, h)
    s.append(text(30, 34, f"コップ (井戸型) の中心位相 {theta:.0f}°　円の半径 = 波から計算した R",
                  size=18))
    s += legend(30, 66, [
        (BLACK, None, 4, f"半径 R(Δ) の円。R(Δ) = {rad:.4f} (倍音のノルム。直流は含めない)"),
        (RED, "2 3", 2, f"コイン (凸、中心 {THETA_COIN:.0f}°)"),
        (GREEN, "2 3", 2, f"コップ (井戸型、中心 {theta:.0f}°)"),
        (BLUE, None, 4, "交差項 2 × コイン × コップ"),
    ])
    # 位相の目盛り。円が小さいので、文字は波形の外に置く
    for deg in sorted({THETA_COIN, theta}):
        a, b = pt(0.0, deg), pt(rad + wc.H + 0.10, deg)
        s.append(f'<line x1="{a[0]:.1f}" y1="{a[1]:.1f}" x2="{b[0]:.1f}" y2="{b[1]:.1f}" '
                 f'stroke="{GREY}" stroke-width="0.8" stroke-dasharray="3 4"/>')
        tx, ty = pt(rad + wc.H + 0.17, deg)
        s.append(text(tx, ty + 5, f"{deg:.0f}°", size=14, fill=GREY, anchor="middle"))
    s.append(f'<circle cx="{cx}" cy="{cy}" r="2.5" fill="{GREY}"/>')
    s.append(f'<circle cx="{cx}" cy="{cy}" r="{rad * PX:.2f}" fill="none" stroke="{BLACK}" '
             f'stroke-width="2"/>')
    s.append(curve(lambda d: wc.well(d, theta), GREEN, 1.3, "2 3"))
    s.append(curve(f_coin, RED, 1.3, "2 3"))
    s.append(curve(lambda d: wc.cross(d, theta), BLUE, 2.2))

    y = 640
    s.append(text(30, y, f"コインの R₁ = √Σ A² = {R1:.4f}　　コップの R₂ = √Σ B² = {R2:.4f}　　"
                         f"比 R₁/R₂ = {R1 / R2:.4f}", size=15))
    s.append(text(30, y + 26, f"和の R² = R₁² + R₂² + 2 Σ A B cos(mΔ)。Δ = {theta - THETA_COIN:+.0f}° で "
                              f"2 Σ A B cos(mΔ) = {cross_sum(theta):+.6f}、R(Δ) = {rad:.4f}",
                  size=15))
    s.append(text(30, y + 52, f"コインの中心 (90°) での交差項 c = {wc.cross(THETA_COIN, theta):+.4f}　"
                              f"(コインの山 {f_coin(THETA_COIN):+.4f}、そこでのコップの値 "
                              f"{wc.well(THETA_COIN, theta):+.4f})", size=15, fill=BLUE))
    s.append(text(30, y + 78, "描画: r = R(Δ) + 値。波形・交差項の数と縮尺は前の図と同じ。"
                              "変えたのは円の半径だけ。",
                  size=14, fill=GREY))
    save(path, s)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    lines = [
        "コインと井戸型コップ: 波から計算した R (倍音のノルム。直流は含めない)",
        f"コインの R1 = {R1:.6f}",
        f"コップの R2 = {R2:.6f}",
        f"比 R1/R2 = {R1 / R2:.6f}",
        f"ピタゴラス和 sqrt(R1^2 + R2^2) = {math.sqrt(R1 * R1 + R2 * R2):.6f}",
        f"コップの直流 c0 = {wc.DC:.6f}、c0/R2 = {wc.DC / R2:.6f}",
        f"コインの山 / R1 = {f_coin(THETA_COIN) / R1:.6f}、井戸の深さ h / R2 = {wc.H / R2:.6f}",
        "",
        "コップの中心位相ごとの値",
        "  位相    Δ      2ΣABcos(mΔ)   和の R     R/ピタゴラス和   c(90度)     C/π",
    ]
    pyt = math.sqrt(R1 * R1 + R2 * R2)
    for theta in range(0, 360, wc.STEP_FIG):
        th = float(theta)
        phase_figure(os.path.join(OUT_DIR, f"fig_cup_well_cross_derivedR_phase{theta:03d}.svg"), th)
        lines.append(f"  {theta:3d} 度  {th - THETA_COIN:+5.0f}  {cross_sum(th):+.6f}     "
                     f"{r_total(th):.6f}   {r_total(th) / pyt:.6f}       "
                     f"{wc.cross(THETA_COIN, th):+.6f}   {wc.cross_total(th) / math.pi:+.6f}")
    out = "\n".join(lines) + "\n"
    with open(RES_PATH, "w", encoding="utf-8") as f:
        f.write(out)
    print(out)


if __name__ == "__main__":
    main()
