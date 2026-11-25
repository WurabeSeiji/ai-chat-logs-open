#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
コップを井戸型にしたときの交差項の図を SVG で出力する (2026-11-25)。

ねらい: 交差項 (2 × コイン × コップ = 情報の AND) だけで「隠れる」を表す。
  交差項は、コインとコップの両方に振幅がある場所にだけ立つ。
  だからコップを「覆う場所で振幅が 0、ほかの場所で振幅がある」波形にすれば、
  覆われたコインは交差項に出ず、覆われていないコインは出る。
  その波形が井戸型 (一つだけの矩形パルスを反転したもの) である。

波形
  コイン  : f(psi) = a R * sum_{m=1}^{35} w_m cos( m (psi - 90度) )      きのうと同じ (凸)
  コップ  : g(psi) = h (1 - W/2π) - sum_{m=1}^{N} w_m (2h / mπ) sin(mW/2) cos( m (psi - theta) )
            幅 W = 30 度、深さ h の井戸をフーリエ級数で合成したもの。外の高さ h、底 0。
  交差項  : c(psi) = 2 f(psi) g(psi)   (R^2)
  和は使わない。

決めた値と、その理由
  W = 30 度      きのうのコップの幅。
  h = 23/72 R    きのうのコップの山の高さ。
  N = 71         フェイエール重みでは壁の幅が 720/(N+1) = 10 度になる。
                 平らな底の幅 = 30 - 10 = 20 度 = コインの幅。コインがちょうど収まる最小の N。
  重み           フェイエール w_m = 1 - m/(N+1)。きのうと同じ。

出力
  figures/well_cross/fig_cup_well_cross_phaseXXX.svg   コップを 0 度から 15 度刻みで一周 (24 枚)
  figures/well_cross/fig_cup_well_cross_summary.svg    コップの位置と交差項の関係
  results/fig_cup_well_cross.txt                       数値

色: テーブル = 黒、コイン = 赤の点線、コップ = 緑の点線、交差項 = 青。
交差項は R^2 の次元なので、R で一度割って長さにして描く (きのうの交差項の図と同じ)。
それ以外の割り算と正規化は入れていない。

実行: python3 fig_cup_well_cross.py
"""
import math
import os

import fig_cup_coin_table as base
from fig_cup_coin_table import (BLACK, RED, GREEN, BLUE, ORANGE, GREY, N_COIN, A_COIN,
                                A_COIN_DEN, THETA_COIN, Circle, f_coin, legend, save,
                                svg_open, text, weight)

OUT_DIR = os.path.join(base.FIG_DIR, "well_cross")
RES_PATH = os.path.join(base.RES_DIR, "fig_cup_well_cross.txt")

W_DEG = 30.0
H_NUM, H_DEN = 23, 72
H = H_NUM / H_DEN
N_WELL = 71
STEP_FIG = 15

_W = math.radians(W_DEG)
DC = H * (1.0 - _W / (2.0 * math.pi))
B = [0.0] + [-(1.0 - m / (N_WELL + 1.0)) * 2.0 * H * math.sin(m * _W / 2.0) / (m * math.pi)
             for m in range(1, N_WELL + 1)]          # B[m]: 井戸の m 次の係数
A = [0.0] + [A_COIN * weight(m, N_COIN) for m in range(1, N_COIN + 1)]   # A[m]: コインの係数


def well(deg, theta):
    x = math.radians(deg - theta)
    return DC + sum(B[m] * math.cos(m * x) for m in range(1, N_WELL + 1))


def cross(deg, theta):
    return 2.0 * f_coin(deg) * well(deg, theta)


def cross_total(theta):
    """C = ∫ c dpsi を倍音の和で (R^2)。コインに DC 成分が無いので井戸の DC は効かない。"""
    d = math.radians(theta - THETA_COIN)
    return 2.0 * math.pi * sum(A[m] * B[m] * math.cos(m * d) for m in range(1, N_COIN + 1))


def cross_total_numeric(theta, n=7200):
    return sum(cross(360.0 * i / n, theta) for i in range(n)) * (2.0 * math.pi / n)


def cross_total_coin_dc(theta):
    """参考: コインに m=0 の項 (a R / 2) を足して、コインをどこでも 0 以上にした場合の C。"""
    return cross_total(theta) + 2.0 * (A_COIN / 2.0) * DC * 2.0 * math.pi


def scan(theta):
    """交差項の波形の最大と最小 (0.05 度刻み)。"""
    vals = [(cross(i * 0.05, theta), i * 0.05) for i in range(7200)]
    return max(vals), min(vals)


def phase_figure(path, theta):
    w, h = 900, 880
    c = Circle(450, 490)
    (cmx, cmx_deg), (cmn, cmn_deg) = scan(theta)
    s = svg_open(w, h)
    s.append(text(30, 34, f"コップ (井戸型) の中心位相 {theta:.0f}°", size=18))
    s += legend(30, 66, [
        (BLACK, None, 5, "テーブル (円周、半径 R)"),
        (RED, "2 3", 2, f"コイン (凸。倍音 1〜{N_COIN}、フェイエール重み、一本あたりの振幅 "
                        f"R/{A_COIN_DEN}、中心 {THETA_COIN:.0f}°)"),
        (GREEN, "2 3", 2, f"コップ (井戸型。幅 {W_DEG:.0f}°、深さ h = {H_NUM}/{H_DEN} R、"
                          f"倍音 0〜{N_WELL}、フェイエール重み、中心 {theta:.0f}°)"),
        (BLUE, None, 5, "交差項 2 × コイン × コップ"),
    ])
    s += c.base({THETA_COIN, theta})
    s.append(c.curve(lambda d: well(d, theta), GREEN, 1.3, "2 3"))
    s.append(c.curve(f_coin, RED, 1.3, "2 3"))
    s.append(c.curve(lambda d: cross(d, theta), BLUE, 2.4))
    y = 786
    s.append(text(30, y, "コップ: g(ψ) = h (1 − W/2π) − Σ (1 − m/(N+1)) (2h / mπ) sin(mW/2) "
                         "cos m(ψ − θ)　(井戸の幅 W、深さ h)", size=15))
    s.append(text(30, y + 26, "交差項: c(ψ) = 2 f(コイン) g(コップ)。R² の次元なので、"
                              "R で一度割って長さにして描く: r(ψ) = R + c(ψ)/R", size=15))
    s.append(text(30, y + 52, f"コインの中心 (90°) での値 c = {cross(THETA_COIN, theta):+.4f} R²　　"
                              f"波形の最大 {cmx:+.4f} R² ({cmx_deg:.1f}°)　"
                              f"最小 {cmn:+.4f} R² ({cmn_deg:.1f}°)", size=15, fill=BLUE))
    s.append(text(30, y + 78, f"コインの中心でのコップの値 g = {well(THETA_COIN, theta):+.4f} R　　"
                              f"コインの山 f = {f_coin(THETA_COIN):+.4f} R", size=15))
    save(path, s)


def summary_figure(path, rows):
    """rows = [(theta, c(90度), 波形の最大, C, コインに DC を足した場合の C), ...] (1 度刻み)。"""
    w, h = 940, 900
    left, right = 100, 900

    def X(deg):
        return left + deg / 360.0 * (right - left)

    def panel(y_top, ph, lo, hi, title, series, ticks, fmt):
        out = [text(left, y_top - 12, title, size=16)]

        def Y(v):
            return y_top + (hi - v) / (hi - lo) * ph
        out.append(f'<rect x="{left}" y="{y_top}" width="{right - left}" height="{ph}" '
                   f'fill="none" stroke="{GREY}" stroke-width="1"/>')
        for v in ticks:
            out.append(f'<line x1="{left}" y1="{Y(v):.1f}" x2="{right}" y2="{Y(v):.1f}" '
                       f'stroke="{BLACK if v == 0 else "#dddddd"}" '
                       f'stroke-width="{1.5 if v == 0 else 1}"/>')
            out.append(text(left - 8, Y(v) + 5, "0" if v == 0 else fmt.format(v), size=13,
                            fill=GREY, anchor="end"))
        for deg in range(0, 361, 30):
            out.append(f'<line x1="{X(deg):.1f}" y1="{y_top + ph}" x2="{X(deg):.1f}" '
                       f'y2="{y_top + ph + 6}" stroke="{GREY}" stroke-width="1.2"/>')
            out.append(text(X(deg), y_top + ph + 21, f"{deg}°", size=13, fill=GREY,
                            anchor="middle"))
        out.append(f'<line x1="{X(THETA_COIN):.1f}" y1="{y_top}" x2="{X(THETA_COIN):.1f}" '
                   f'y2="{y_top + ph}" stroke="{RED}" stroke-width="1" stroke-dasharray="4 4"/>')
        for idx, color, width_px, dash in series:
            pts = " ".join(f"{X(r[0]):.2f},{Y(r[idx]):.2f}" for r in rows)
            d = f' stroke-dasharray="{dash}"' if dash else ""
            out.append(f'<polyline points="{pts}" fill="none" stroke="{color}" '
                       f'stroke-width="{width_px}" stroke-linejoin="round"{d}/>')
        return out

    s = svg_open(w, h)
    s.append(text(30, 34, "コップ (井戸型) の位置と交差項の関係", size=18))
    s.append(text(30, 62, f"横はコップの中心位相。コインは {THETA_COIN:.0f}° に固定 (赤の縦の破線)。"
                          f"井戸の幅 {W_DEG:.0f}°、深さ h = {H:.4f} R、倍音 0〜{N_WELL}。",
                  size=14, fill=GREY))

    s += legend(30, 96, [
        (BLUE, None, 3, "コインの中心 (90°) での交差項 c(90°) = 2 f g"),
        (BLUE, "5 4", 1.6, "交差項の波形の最大"),
    ], size=14)
    s += panel(170, 250, -0.02, 0.18,
               "1. 交差項の波形 (その場所での値。単位 R²)",
               [(2, BLUE, 1.6, "5 4"), (1, BLUE, 2.6, None)],
               [0, 0.05, 0.10, 0.15], "{:+.2f}")

    s += legend(30, 496, [
        (ORANGE, None, 3, "C = ∫ c dψ (コインはきのうのまま。DC 成分なし)"),
        (ORANGE, "5 4", 1.6, "参考: コインに m = 0 の項 (a R / 2) を足した場合の C"),
    ], size=14)
    s += panel(570, 250, -0.03, 0.035,
               "2. 交差項を円周全体で足した値 C (単位 R²)",
               [(4, ORANGE, 1.6, "5 4"), (3, ORANGE, 2.6, None)],
               [-0.02, -0.01, 0, 0.01, 0.02, 0.03], "{:+.2f}")
    s.append(text(30, 872, "コインに DC 成分が無いと、コインを一周で足した値が 0 なので、"
                           "覆われていないコインは C に出ない (波形には出る)。", size=14, fill=GREY))
    save(path, s)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(base.RES_DIR, exist_ok=True)

    for theta in range(0, 360, STEP_FIG):
        phase_figure(os.path.join(OUT_DIR, f"fig_cup_well_cross_phase{theta:03d}.svg"),
                     float(theta))

    rows = []
    for theta in range(0, 361):
        # 波形の最大はコインの山の中 (70〜110 度) にあるので、そこだけ 0.1 度刻みで探す
        cmx = max(cross(70.0 + 0.1 * i, float(theta)) for i in range(401))
        rows.append([float(theta), cross(THETA_COIN, float(theta)), cmx,
                     cross_total(float(theta)), cross_total_coin_dc(float(theta))])
    summary_figure(os.path.join(OUT_DIR, "fig_cup_well_cross_summary.svg"), rows)

    free = 2.0 * f_coin(THETA_COIN) * H
    lines = [
        "コップを井戸型にしたときの交差項 c = 2 × コイン × コップ",
        f"コイン: 倍音 1〜{N_COIN}, a = 1/{A_COIN_DEN}, フェイエール重み, 中心 {THETA_COIN:.0f} 度, "
        f"山 {f_coin(THETA_COIN):+.6f} R",
        f"コップ: 井戸の幅 {W_DEG:.0f} 度, 深さ h = {H_NUM}/{H_DEN} R = {H:.6f} R, 倍音 0〜{N_WELL}, "
        f"フェイエール重み, DC 成分 {DC:.6f} R",
        f"  井戸の底 (中心) {well(0.0, 0.0):+.6f} R、中心から 10 度 {well(10.0, 0.0):+.6f} R、"
        f"15 度 {well(15.0, 0.0):+.6f} R、20 度 {well(20.0, 0.0):+.6f} R",
        f"覆わないときの目安 2 × コインの山 × h = {free:.6f} R^2",
        f"C の検算 (コップ 90 度): 倍音の和 {cross_total(90.0):+.8f}、数値積分 "
        f"{cross_total_numeric(90.0):+.8f}",
        "",
        "コップの中心位相ごとの値 (単位 R^2)",
        "  位相    c(90度)   波形の最大 (位置)      波形の最小 (位置)      C (DCなし)   C (コインにDC)",
    ]
    for theta in range(0, 360, STEP_FIG):
        (cmx, cmx_deg), (cmn, cmn_deg) = scan(float(theta))
        lines.append(f"  {theta:3d} 度  {cross(THETA_COIN, float(theta)):+.6f}  "
                     f"{cmx:+.6f} ({cmx_deg:6.2f} 度)  {cmn:+.6f} ({cmn_deg:6.2f} 度)  "
                     f"{cross_total(float(theta)):+.6f}   {cross_total_coin_dc(float(theta)):+.6f}")
    lines += ["", "コインの中心での交差項が、覆わないときの目安の半分になるコップの位置"]
    prev = None
    for r in rows:
        if prev is not None and (prev[1] - free / 2) * (r[1] - free / 2) < 0:
            t = (free / 2 - prev[1]) / (r[1] - prev[1])
            lines.append(f"  {prev[0] + t:.2f} 度")
        prev = r
    out = "\n".join(lines) + "\n"
    with open(RES_PATH, "w", encoding="utf-8") as f:
        f.write(out)
    print(out)


if __name__ == "__main__":
    main()
