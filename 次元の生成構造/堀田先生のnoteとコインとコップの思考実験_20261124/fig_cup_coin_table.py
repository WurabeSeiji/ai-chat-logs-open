#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
コインとコップの思考実験の図を SVG で出力する。

  模式図 (横から見た断面)
    figures/fig_cup_coin_schematic_covered.svg  コップをコインの上に被せた状態
    figures/fig_cup_coin_schematic_beside.svg   コップをコインの隣に置いた状態

  S^1 上の精密図 (テーブル = 半径 R の円周)
    figures/fig_cup_coin_s1_phase90.svg         コップの中心位相 90 度 (コインと重なる)
    figures/fig_cup_coin_s1_phase60.svg         コップの中心位相 60 度 (ずらした状態)

  交差項の図 (コップを 0 度から 15 度刻みで一周させる。24 枚)
    figures/cross_term/fig_cup_coin_s1_cross_phaseXXX.svg

色: テーブル = 黒、コイン = 赤、コップ = 緑、交差項 = 青。

量はすべて R を単位にする (位相の幅も振幅も、全体 R に対する比で意味を持つ)。

波形
  f(psi) = a R * sum_{m=1}^{N} w_m cos( m (psi - theta) )
    コイン: N = 35, a = 1/72, theta = 90 度
    コップ: N = 23, a = 1/36, theta = 90 度 または 60 度 (交差項の図では 0〜345 度)
  a は倍音一本あたりの振幅を R で測った値。コインはコップの半分
  (コインの山がコップの山からはみ出さないため)。
  倍音ごとの重み w_m は WEIGHT で選ぶ。
    "fejer" : w_m = 1 - m/(N+1)   (フェイエール。高い倍音ほど振幅を下げる。山の脇の波が小さくなる)
    "equal" : w_m = 1             (等振幅。山の脇に波が残る)
  フェイエールでは山の幅 (底から底) が 720/(N+1) 度になり、等振幅の 2 倍に広がる。
  幅をコイン 20 度、コップ 30 度に保つため、倍音の上限を 35 と 23 にしてある。
  (等振幅で同じ幅にするなら 18 と 12。)
  山の高さはフェイエールで a R N/2 (コイン 0.243 R、コップ 0.319 R)。
  定数項 (m = 0) は入れていない。このためフェイエールでは、山の外の波形は
  円周より a R/2 だけ内側に平らに沿う (コイン 0.0069 R、コップ 0.0139 R)。

描画
  波形は、中心からの距離 r(psi) = R + f(psi) として描く。
  f は R を単位にした長さなので、表示のための縮尺は入れていない。

交差項
  c(psi) = 2 * f_coin(psi) * f_cup(psi)     ( (a + b)^2 = a^2 + b^2 + 2ab の 2ab )
  c は R^2 の次元を持つ。波形と同じ図に描くには長さにする必要があるので、
  R で一度だけ割って  r(psi) = R + c(psi)/R  として描く。
  この割り算は次元を合わせるためのもので、図ごとに高さを揃える操作はしていない。
  各図に、円周全体で足した値  C = ∫_0^{2π} c(psi) dpsi  (psi はラジアン、単位 R^2) を書き入れる。
  C は倍音の和  2π R^2 * Σ_m a_coin a_cup w_m^coin w_m^cup cos(m (θ_cup − θ_coin))  とも
  書けるので、数値積分と照らし合わせて results に出す。

交差項のまとめの図
    figures/fig_cup_coin_cross_summary.svg
  横軸にコップの中心位相 (0〜360 度、1 度刻み)、縦軸に C (単位 R^2) を取る。
  15 度刻みの点は、交差項の図 24 枚に対応する。
  灰色の帯は振幅の分解能 ±EPS_C (単位 R^2)。帯の内側は読めない範囲。
  EPS_C は次の仮説の式で決める (導いたものではない。図にも仮説と書く):
    振幅の分解能 = k × (最大の振幅) ÷ N
      k は重みで決まる無次元の係数 (等振幅 1、フェイエール 2。山の幅が広がる倍率と同じ)
      N は倍音の数
    波形に当てると、k × A ÷ N は倍音一本あたりの振幅 a R に等しい (results に出す)。
    交差項 C に当てるときは、最大の振幅に C の最大 (重なったときの値) を、
    N に C に含まれる倍音の数 min(N_coin, N_cup) を使う。

偶数倍音と奇数倍音の図 (窓はコップ一つ、中心 90 度)
    figures/parity/fig_cup_s1_parity_even.svg   偶数倍音だけ (m = 2, 4, ..., 22)
    figures/parity/fig_cup_s1_parity_odd.svg    奇数倍音だけ (m = 1, 3, ..., 23)
    figures/parity/fig_cup_s1_parity_both.svg   両方 (m = 1, ..., 23)
  重み w_m = 1 - m/(N+1) (N = 23) と一本あたりの振幅 R/36 は三枚とも同じで、
  足す倍音を選ぶだけ。高さを揃える操作はしていない。
  「両方」の波形は「偶数だけ」と「奇数だけ」の和に等しい。

開いた系にできない理由の図 (窓はコップ一つ、中心 90 度)
    figures/open_system/fig_cup_s1_all_harmonics.svg   円周の上に、倍音 1〜23 を一本ずつ点線で描く
    figures/open_system/fig_cup_cut_open_270.svg       270 度で切り開いて直線にした図
  倍音 m 番の波形は  h_m(psi) = a R w_m cos( m (psi - theta) )。
  奇数倍音は紫、偶数倍音は橙の点線、合成波形は緑の実線。
  円周の図では、倍音一本ずつは小さくて見えないので、倍音の点線だけ半径方向に
  HARM_MAG 倍して描く (合成波形は等倍)。倍率は図に書く。
  切り開いた図は、縦軸に値 (R 単位) の目盛りを付けた二段のグラフで、倍率は使わない。
  切り口 (270 度) での各倍音の値と傾きを results に出す。

両端を固定した直線の上の定常波
    figures/open_system/fig_line_fixed_ends_equal16.svg
  長さ L = 2πR の直線。両端は固定端 (節)。k 番の定常波は sin(k π x / L)  (k = 1, 2, ...)。
  k = 1〜N_LINE を同じ振幅 a R (a = 1/A_LINE_DEN) で足す:
    g(x) = a R * sum_{k=1}^{N_LINE} sin(k π x / L)
  横軸は x/R (0 〜 2π)。縦軸は R を単位にした値。倍率は使わない。
    figures/open_system/fig_line_fixed_ends_fejer16.svg
  同じ直線で、フェイエール重み w_k = 1 - k/(N_LINE+1) を掛けた場合。
  比較のため、等振幅の波形 (灰色の破線) と、
  裾の曲線 (a R / 2) cot(π x / (2 L)) (黒の点線) を重ねる。
  裾の曲線は、定常波を何番まで足しても、また重みで波打ちを消しても残る部分である。

外部ライブラリは使わない。実行: python3 fig_cup_coin_table.py
"""

import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(HERE, "figures")
RES_DIR = os.path.join(HERE, "results")

BLACK = "#000000"
RED = "#d62728"
GREEN = "#2ca02c"
BLUE = "#1f77b4"
PURPLE = "#9467bd"
ORANGE = "#ff7f0e"
GREY = "#666666"
FONT = "'Hiragino Sans','Yu Gothic','Noto Sans CJK JP','Meiryo',sans-serif"

# ---- 設定 (長さはすべて R を単位にする) ------------------------------------
WEIGHT = "fejer"     # 倍音の重み。"fejer" または "equal"
N_COIN = 35          # コインの倍音の上限 (フェイエールで山の幅 720/36 = 20 度)
N_CUP = 23           # コップの倍音の上限 (フェイエールで山の幅 720/24 = 30 度)
A_COIN_DEN = 72      # コインの倍音一本あたりの振幅 = R/72
A_CUP_DEN = 36       # コップの倍音一本あたりの振幅 = R/36
A_COIN = 1.0 / A_COIN_DEN
A_CUP = 1.0 / A_CUP_DEN
THETA_COIN = 90.0    # コインの中心位相 (度)。北極
THETA_CUP_LIST = [90.0, 60.0]   # 精密図でのコップの中心位相 (度)
CROSS_STEP = 15      # 交差項の図でコップを動かす刻み (度)
KAPPA = {"fejer": 2.0, "equal": 1.0}[WEIGHT]   # 振幅の分解能の式の係数 k (山の幅の倍率と同じ)
N_LINE = 16          # 両端固定の直線に載せる定常波の上限
A_LINE_DEN = 36      # その一本あたりの振幅 = R/36 (コップと同じ)
HARM_MAG = 5         # 円周の図で倍音の点線を半径方向に拡大する倍率 (合成波形は等倍)
SAMPLES = 3600       # 円周の分割数 (0.1 度刻み)
R_PX = 200.0         # 図の上での R の長さ (ピクセル)

WNAME = {"fejer": "フェイエール重み", "equal": "等振幅"}[WEIGHT]
WFORM = "(1 − m/(N+1)) " if WEIGHT == "fejer" else ""
SUM_SUB = '<tspan dy="5" font-size="11">m=1…N</tspan><tspan dy="-5"> </tspan>'


def weight(m, n_max):
    if WEIGHT == "fejer":
        return 1.0 - m / (n_max + 1.0)
    if WEIGHT == "equal":
        return 1.0
    raise ValueError(WEIGHT)


def synth(psi_deg, n_max, theta_deg, amp):
    """倍音 1..n_max の合成波形。値は R を単位にした長さ。"""
    x = math.radians(psi_deg - theta_deg)
    return amp * sum(weight(m, n_max) * math.cos(m * x) for m in range(1, n_max + 1))


def f_coin(deg):
    return synth(deg, N_COIN, THETA_COIN, A_COIN)


def f_cup(deg, theta_cup):
    return synth(deg, N_CUP, theta_cup, A_CUP)


def profile(n_max, amp):
    """山の高さ、山の半幅 (中心から最初の極小まで)、山の外側の最大と最小。"""
    step = 0.005
    vals = [synth(i * step, n_max, 0.0, amp) for i in range(int(180.0 / step) + 1)]
    i = 1
    while i < len(vals) and vals[i] <= vals[i - 1]:
        i += 1
    i0 = i - 1
    outside = vals[i0:]
    return vals[0], i0 * step, max(outside), min(outside)


# ---- SVG の部品 ------------------------------------------------------------
def svg_open(w, h):
    return [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
        f'viewBox="0 0 {w} {h}">',
        f'<rect x="0" y="0" width="{w}" height="{h}" fill="#ffffff"/>',
    ]


def text(x, y, s, size=16, fill=BLACK, anchor="start"):
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" font-size="{size}" '
            f'fill="{fill}" text-anchor="{anchor}">{s}</text>')


def legend(x, y, rows, size=15):
    """rows = [(色, 点線の指定または None, 線の太さ, 文字列), ...] を縦に並べる。"""
    out = []
    for i, (color, dash, lw, label) in enumerate(rows):
        yy = y + i * (size + 9)
        d = f' stroke-dasharray="{dash}"' if dash else ""
        out.append(f'<line x1="{x}" y1="{yy - 5}" x2="{x + 34}" y2="{yy - 5}" '
                   f'stroke="{color}" stroke-width="{lw}"{d}/>')
        out.append(text(x + 44, yy, label, size=size))
    return out


def save(path, lines):
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines + ["</svg>", ""]))


# ---- 模式図 ----------------------------------------------------------------
def schematic(path, cup_center_x, title):
    w, h = 680, 330
    table_y = 250           # テーブルの上面
    table_x0, table_x1 = 40, 640
    table_t = 16
    coin_cx, coin_w, coin_h = 250, 80, 12
    cup_inner_w, cup_h, cup_stroke = 150, 130, 9

    s = svg_open(w, h)
    s.append(text(40, 34, title, size=18))
    s += legend(40, 66, [(BLACK, None, 5, "テーブル"), (RED, None, 5, "コイン"),
                         (GREEN, None, 5, "コップ")])

    # テーブル (黒)
    s.append(f'<rect x="{table_x0}" y="{table_y}" width="{table_x1 - table_x0}" '
             f'height="{table_t}" fill="{BLACK}"/>')
    # コイン (赤)。テーブルの上に置く
    s.append(f'<rect x="{coin_cx - coin_w / 2}" y="{table_y - coin_h}" width="{coin_w}" '
             f'height="{coin_h}" rx="3" fill="{RED}"/>')
    # コップ (緑)。伏せた断面を線で描く。中が見えるように塗らない
    x0 = cup_center_x - cup_inner_w / 2
    x1 = cup_center_x + cup_inner_w / 2
    top = table_y - cup_h
    base = table_y   # 縁はテーブルの上面に接する
    r = 14
    d = (f'M {x0:.1f} {base:.1f} L {x0:.1f} {top + r:.1f} '
         f'Q {x0:.1f} {top:.1f} {x0 + r:.1f} {top:.1f} '
         f'L {x1 - r:.1f} {top:.1f} '
         f'Q {x1:.1f} {top:.1f} {x1:.1f} {top + r:.1f} '
         f'L {x1:.1f} {base:.1f}')
    s.append(f'<path d="{d}" fill="none" stroke="{GREEN}" stroke-width="{cup_stroke}" '
             f'stroke-linecap="butt" stroke-linejoin="round"/>')
    save(path, s)


# ---- S^1 上の図の共通部品 --------------------------------------------------
class Circle:
    def __init__(self, cx, cy):
        self.cx, self.cy = cx, cy

    def pt(self, r_over_R, deg):
        """中心からの距離 (R 単位) と角度から画面上の点へ。90 度が上、反時計回り。"""
        a = math.radians(deg)
        r = r_over_R * R_PX
        return self.cx + r * math.cos(a), self.cy - r * math.sin(a)

    def curve(self, func, color, width, dash=None):
        """r(psi) = R + func(psi)。func は R を単位にした長さを返す。"""
        pts = []
        for i in range(SAMPLES + 1):
            deg = 360.0 * i / SAMPLES
            x, y = self.pt(1.0 + func(deg), deg)
            pts.append(f"{x:.2f},{y:.2f}")
        d = f' stroke-dasharray="{dash}"' if dash else ""
        return (f'<polyline points="{" ".join(pts)}" fill="none" stroke="{color}" '
                f'stroke-width="{width}" stroke-linejoin="round"{d}/>')

    def base(self, marks, label_r=0.85):
        """テーブルの円、半径 R の表示、位相の目盛り (marks は度の集まり)。"""
        s = []
        ex, ey = self.pt(1.0, 210.0)
        s.append(f'<line x1="{self.cx}" y1="{self.cy}" x2="{ex:.1f}" y2="{ey:.1f}" '
                 f'stroke="{GREY}" stroke-width="1.2" stroke-dasharray="6 5"/>')
        s.append(f'<circle cx="{self.cx}" cy="{self.cy}" r="3" fill="{GREY}"/>')
        mx, my = self.pt(0.5, 210.0)
        s.append(text(mx + 6, my - 8, "R", size=18, fill=GREY))
        for deg in sorted(set(marks)):
            a, b = self.pt(0.93, deg), self.pt(0.99, deg)
            s.append(f'<line x1="{a[0]:.1f}" y1="{a[1]:.1f}" x2="{b[0]:.1f}" y2="{b[1]:.1f}" '
                     f'stroke="{GREY}" stroke-width="1.5"/>')
            tx, ty = self.pt(label_r, deg)
            s.append(text(tx, ty + 6, f"{deg:.0f}°", size=15, fill=GREY, anchor="middle"))
        s.append(f'<circle cx="{self.cx}" cy="{self.cy}" r="{R_PX}" fill="none" '
                 f'stroke="{BLACK}" stroke-width="2.5"/>')
        return s


def coin_label(theta=THETA_COIN):
    return (f"コイン (倍音 1〜{N_COIN}、{WNAME}、一本あたりの振幅 R/{A_COIN_DEN}、"
            f"中心 {theta:.0f}°)")


def cup_label(theta):
    return (f"コップ (倍音 1〜{N_CUP}、{WNAME}、一本あたりの振幅 R/{A_CUP_DEN}、"
            f"中心 {theta:.0f}°)")


def wave_formula_lines(theta_cup, y0):
    hc, _, _, _ = profile(N_COIN, A_COIN)
    hk, _, _, _ = profile(N_CUP, A_CUP)
    return [
        text(30, y0, f"波形: f(ψ) = a R Σ{SUM_SUB}{WFORM}cos m(ψ − θ)　　"
                     f"描画: 中心からの距離 r(ψ) = R + f(ψ)", size=15),
        text(30, y0 + 26, f"コイン: N = {N_COIN}、a = 1/{A_COIN_DEN}、θ = {THETA_COIN:.0f}°、"
                          f"山の高さ {hc:.3f} R　　"
                          f"コップ: N = {N_CUP}、a = 1/{A_CUP_DEN}、θ = {theta_cup:.0f}°、"
                          f"山の高さ {hk:.3f} R", size=15),
    ]


# ---- S^1 上の精密図 --------------------------------------------------------
def s1_figure(path, theta_cup):
    w, h = 900, 720
    c = Circle(450, 390)
    s = svg_open(w, h)
    s.append(text(30, 34, f"コップの中心位相 {theta_cup:.0f}°", size=18))
    s += legend(30, 66, [
        (BLACK, None, 5, "テーブル (円周、半径 R)"),
        (RED, None, 5, coin_label()),
        (GREEN, None, 5, cup_label(theta_cup)),
    ])
    s += c.base({THETA_COIN, theta_cup})
    s.append(c.curve(lambda d: f_cup(d, theta_cup), GREEN, 2.2))
    s.append(c.curve(f_coin, RED, 2.2))
    s += wave_formula_lines(theta_cup, 650)
    save(path, s)


# ---- 交差項 ----------------------------------------------------------------
def cross(deg, theta_cup):
    """交差項 c(psi)。R^2 を単位にした値。"""
    return 2.0 * f_coin(deg) * f_cup(deg, theta_cup)


def cross_total_numeric(theta_cup):
    """C = ∫ c dpsi を数値で (等間隔の和。周期関数なので台形則と同じ)。単位 R^2。"""
    n = 7200
    tot = sum(cross(360.0 * i / n, theta_cup) for i in range(n))
    return tot * (2.0 * math.pi / n)


def cross_total_harmonic(theta_cup):
    """同じ C を倍音の和で。単位 R^2。"""
    d = math.radians(theta_cup - THETA_COIN)
    tot = 0.0
    for m in range(1, min(N_COIN, N_CUP) + 1):
        tot += A_COIN * A_CUP * weight(m, N_COIN) * weight(m, N_CUP) * math.cos(m * d)
    return 2.0 * math.pi * tot


def s1_cross_figure(path, theta_cup):
    w, h = 900, 800
    c = Circle(450, 440)
    s = svg_open(w, h)
    s.append(text(30, 34, f"コップの中心位相 {theta_cup:.0f}°", size=18))
    s += legend(30, 66, [
        (BLACK, None, 5, "テーブル (円周、半径 R)"),
        (RED, "2 3", 2, coin_label()),
        (GREEN, "2 3", 2, cup_label(theta_cup)),
        (BLUE, None, 5, "交差項 2 × コイン × コップ"),
    ])
    s += c.base({THETA_COIN, theta_cup})
    s.append(c.curve(lambda d: f_cup(d, theta_cup), GREEN, 1.3, "2 3"))
    s.append(c.curve(f_coin, RED, 1.3, "2 3"))
    # c は R^2 単位。R で一度割って長さ (R 単位) にして描く
    s.append(c.curve(lambda d: cross(d, theta_cup), BLUE, 2.4))

    c_num = cross_total_numeric(theta_cup)
    s += wave_formula_lines(theta_cup, 696)
    s.append(text(30, 748, "交差項: c(ψ) = 2 f(コイン) f(コップ)。R² の次元なので、"
                           "R で一度割って長さにして描く: r(ψ) = R + c(ψ)/R", size=15))
    s.append(text(30, 774, f"円周全体で足した値  C = ∫ c(ψ) dψ = {c_num:+.5f} R²",
                  size=15, fill=BLUE))
    save(path, s)
    return c_num


def cross_figures():
    out_dir = os.path.join(FIG_DIR, "cross_term")
    os.makedirs(out_dir, exist_ok=True)
    lines = [
        "交差項 c(psi) = 2 f_coin(psi) f_cup(psi)、C = ∫_0^{2π} c dpsi。単位はすべて R^2",
        f"重み: {WEIGHT}   コイン N = {N_COIN}, a = 1/{A_COIN_DEN}, 中心 {THETA_COIN:g} 度   "
        f"コップ N = {N_CUP}, a = 1/{A_CUP_DEN}",
        "図では c を R で一度割って長さにし、r = R + c/R として描いている",
        "",
        "コップ中心(度)  C(数値積分)     C(倍音の和)     差            c の最大      c の最小",
    ]
    for th in range(0, 360, CROSS_STEP):
        path = os.path.join(out_dir, f"fig_cup_coin_s1_cross_phase{th:03d}.svg")
        c_num = s1_cross_figure(path, float(th))
        c_har = cross_total_harmonic(float(th))
        vals = [cross(360.0 * i / SAMPLES, float(th)) for i in range(SAMPLES)]
        lines.append(f"{th:10d}  {c_num:+14.8f}  {c_har:+14.8f}  {c_num - c_har:12.2e}  "
                     f"{max(vals):+12.8f}  {min(vals):+12.8f}")
    lines += ["", f"離れたときに C が近づく値  −π a_coin a_cup = {-math.pi * A_COIN * A_CUP:+.8f}"]
    out = "\n".join(lines) + "\n"
    with open(os.path.join(RES_DIR, "fig_cup_coin_cross_term.txt"), "w", encoding="utf-8") as f:
        f.write(out)
    print(out)


def amp_resolution_wave(n_max, amp):
    """波形の振幅の分解能 (仮説): k × 山の高さ ÷ N。R 単位。"""
    peak, _, _, _ = profile(n_max, amp)
    return KAPPA * peak / n_max


def amp_resolution_cross():
    """交差項 C の振幅の分解能 (仮説): k × C の最大 ÷ min(N_coin, N_cup)。R^2 単位。"""
    return KAPPA * cross_total_harmonic(THETA_COIN) / min(N_COIN, N_CUP)


def cross_summary_figure(path):
    EPS_C = amp_resolution_cross()
    n_common = min(N_COIN, N_CUP)
    w, h = 900, 664
    x0, x1, y0, y1 = 100, 860, 80, 450      # 描画範囲 (画面上)
    vmin, vmax = -0.005, 0.025              # 縦軸の範囲 (R^2)

    def X(deg):
        return x0 + (x1 - x0) * deg / 360.0

    def Y(v):
        return y1 - (y1 - y0) * (v - vmin) / (vmax - vmin)

    s = svg_open(w, h)
    s.append(text(30, 34, "コップの中心位相と、円周全体で足した交差項 C", size=18))

    # 振幅の分解能の帯
    s.append(f'<rect x="{x0}" y="{Y(EPS_C):.1f}" width="{x1 - x0}" '
             f'height="{Y(-EPS_C) - Y(EPS_C):.1f}" fill="#d9d9d9"/>')
    # 目盛りと枠
    v = vmin
    while v <= vmax + 1e-12:
        yy = Y(v)
        s.append(f'<line x1="{x0}" y1="{yy:.1f}" x2="{x1}" y2="{yy:.1f}" stroke="#bbbbbb" '
                 f'stroke-width="{1.4 if abs(v) < 1e-12 else 0.6}"/>')
        s.append(text(x0 - 8, yy + 5, f"{v:+.3f}" if abs(v) > 1e-12 else "0", size=13,
                      fill=GREY, anchor="end"))
        v += 0.005
    for deg in range(0, 361, 45):
        xx = X(deg)
        s.append(f'<line x1="{xx:.1f}" y1="{y1}" x2="{xx:.1f}" y2="{y1 + 6}" stroke="{GREY}" '
                 f'stroke-width="1"/>')
        s.append(text(xx, y1 + 22, f"{deg}°", size=13, fill=GREY, anchor="middle"))
    s.append(f'<rect x="{x0}" y="{y0}" width="{x1 - x0}" height="{y1 - y0}" fill="none" '
             f'stroke="{BLACK}" stroke-width="1.2"/>')
    s.append(text((x0 + x1) / 2, y1 + 46, "コップの中心位相", size=15, anchor="middle"))
    s.append(f'<text x="34" y="{(y0 + y1) / 2:.1f}" font-family="{FONT}" font-size="15" '
             f'text-anchor="middle" transform="rotate(-90 34 {(y0 + y1) / 2:.1f})">'
             f'C (単位 R²)</text>')

    # コインの中心位相
    s.append(f'<line x1="{X(THETA_COIN):.1f}" y1="{y0}" x2="{X(THETA_COIN):.1f}" y2="{y1}" '
             f'stroke="{RED}" stroke-width="1.3" stroke-dasharray="2 3"/>')
    s.append(text(X(THETA_COIN) + 8, y0 + 18, f"コインの中心 {THETA_COIN:.0f}°", size=13, fill=RED))

    # 曲線 (1 度刻み) と 15 度刻みの点
    pts = [f"{X(d):.2f},{Y(cross_total_harmonic(float(d))):.2f}" for d in range(0, 361)]
    s.append(f'<polyline points="{" ".join(pts)}" fill="none" stroke="{BLUE}" '
             f'stroke-width="2.2" stroke-linejoin="round"/>')
    for d in range(0, 360, CROSS_STEP):
        s.append(f'<circle cx="{X(d):.2f}" cy="{Y(cross_total_harmonic(float(d))):.2f}" '
                 f'r="3.6" fill="{BLUE}"/>')

    # 帯の外に出る範囲 (読める範囲) を 0.01 度刻みで探す
    lo = hi = None
    for i in range(0, 36000):
        d = i / 100.0
        if cross_total_harmonic(d) > EPS_C:
            if lo is None:
                lo = d
            hi = d
    c_peak = cross_total_harmonic(THETA_COIN)
    c_far = -math.pi * A_COIN * A_CUP

    s.append(text(X(THETA_COIN) + 14, Y(c_peak) + 5, f"重なったとき {c_peak:+.4f} R²",
                  size=14, fill=BLUE))
    s.append(text(X(120), Y(-EPS_C) + 20,
                  f"離れたとき {c_far:+.4f} R² に近づく (0 ではない)", size=14, fill=BLUE))
    s.append(text(X(120), Y(EPS_C) - 8,
                  f"灰色の帯: 振幅の分解能 ±{EPS_C:.5f} R² (仮説)。内側は読めない", size=14))

    s += legend(30, 524, [
        (BLUE, None, 5, "C = ∫ c(ψ) dψ、c(ψ) = 2 f(コイン) f(コップ)。点は 15° 刻み (交差項の図 24 枚に対応)"),
    ])
    s.append(text(30, 552, f"コイン: 倍音 1〜{N_COIN}、{WNAME}、一本あたりの振幅 R/{A_COIN_DEN}、"
                           f"中心 {THETA_COIN:.0f}°", size=14))
    s.append(text(30, 578, f"コップ: 倍音 1〜{N_CUP}、{WNAME}、一本あたりの振幅 R/{A_CUP_DEN}",
                  size=14))
    if lo is not None:
        s.append(text(30, 606, f"C が帯の外に出るコップの位置: {lo:.2f}° 〜 {hi:.2f}° "
                               f"(コインの中心から ±{(hi - lo) / 2:.2f}°)", size=14))
    s.append(text(30, 634, f"帯の幅 (仮説): 振幅の分解能 = k × (C の最大) ÷ N = "
                           f"{KAPPA:g} × {c_peak:.5f} ÷ {n_common} = {EPS_C:.5f} R²。"
                           f"k は重みの係数、N は倍音の数", size=14, fill=GREY))
    save(path, s)

    lines = ["コップの中心位相 (度) と C (単位 R^2)。1 度刻み。C は倍音の和で計算",
             f"振幅の分解能 (仮説) EPS_C = k × C最大 ÷ N = {KAPPA:g} × {c_peak:.8f} ÷ {n_common} = {EPS_C:.8f} R^2",
             f"  離れたときの値 ÷ 分解能 = {c_far / EPS_C:+.4f}",
             f"重なったときの C = {c_peak:+.8f}   離れたときに近づく値 = {c_far:+.8f}"]
    if lo is not None:
        lines.append(f"C > EPS_C となる範囲: {lo:.2f} 度 〜 {hi:.2f} 度 (幅 {hi - lo:.2f} 度)")
    lines.append("")
    lines += [f"{d:4d}  {cross_total_harmonic(float(d)):+.8f}" for d in range(0, 360)]
    with open(os.path.join(RES_DIR, "fig_cup_coin_cross_summary.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines[:5]))


# ---- 偶数倍音と奇数倍音 ----------------------------------------------------
PARITY_THETA = 90.0
PARITY_CASES = [
    ("even", "偶数倍音だけ", lambda m: m % 2 == 0),
    ("odd", "奇数倍音だけ", lambda m: m % 2 == 1),
    ("both", "偶数倍音と奇数倍音の両方", lambda m: True),
]


def synth_sel(psi_deg, sel):
    """コップの波形から、sel(m) が真の倍音だけを足したもの。値は R を単位にした長さ。"""
    x = math.radians(psi_deg - PARITY_THETA)
    return A_CUP * sum(weight(m, N_CUP) * math.cos(m * x)
                       for m in range(1, N_CUP + 1) if sel(m))


def parity_figures():
    out_dir = os.path.join(FIG_DIR, "parity")
    os.makedirs(out_dir, exist_ok=True)
    opp = PARITY_THETA + 180.0
    lines = ["偶数倍音と奇数倍音。窓はコップ一つ、中心 90 度。長さは R を単位にした値",
             f"重み: {WEIGHT} (N = {N_CUP})、一本あたりの振幅 R/{A_CUP_DEN}。足す倍音を選ぶだけ", ""]
    vals = {}
    for key, title, sel in PARITY_CASES:
        ms = [m for m in range(1, N_CUP + 1) if sel(m)]
        v0, v1 = synth_sel(PARITY_THETA, sel), synth_sel(opp, sel)
        vals[key] = (v0, v1)
        if key == "both":
            mtxt = f"m = 1, 2, …, {ms[-1]}"
        else:
            mtxt = f"m = {ms[0]}, {ms[1]}, …, {ms[-1]}"

        w, h = 900, 720
        c = Circle(450, 380)
        s = svg_open(w, h)
        s.append(text(30, 34, f"{title} ({mtxt})", size=18))
        s += legend(30, 66, [
            (BLACK, None, 5, "テーブル (円周、半径 R)"),
            (GREEN, None, 5, f"コップ ({WNAME} N = {N_CUP}、一本あたりの振幅 R/{A_CUP_DEN}、"
                             f"中心 {PARITY_THETA:.0f}°)"),
        ])
        s += c.base({PARITY_THETA, opp}, label_r=0.74)   # 内向きの凹みと重ならない位置
        s.append(c.curve(lambda d, sel=sel: synth_sel(d, sel), GREEN, 2.2))
        s.append(text(30, 640, f"波形: f(ψ) = a R Σ {WFORM}cos m(ψ − θ)、Σ は {mtxt} について。"
                               f"a = 1/{A_CUP_DEN}、N = {N_CUP}、θ = {PARITY_THETA:.0f}°", size=15))
        s.append(text(30, 666, f"描画: 中心からの距離 r(ψ) = R + f(ψ)　　"
                               f"f({PARITY_THETA:.0f}°) = {v0:+.4f} R、"
                               f"f({opp:.0f}°) = {v1:+.4f} R", size=15))
        if key == "both":
            s.append(text(30, 692, "この波形は、偶数倍音だけの波形と奇数倍音だけの波形の和に等しい。"
                                   f"{opp:.0f}° では両者が打ち消し合う", size=15))
        elif key == "even":
            s.append(text(30, 692, f"{opp:.0f}° にも、{PARITY_THETA:.0f}° と同じ向きで同じ高さの山が立つ",
                          size=15))
        else:
            s.append(text(30, 692, f"{opp:.0f}° には、{PARITY_THETA:.0f}° と逆向きで同じ大きさの凹みが立つ",
                          size=15))
        save(os.path.join(out_dir, f"fig_cup_s1_parity_{key}.svg"), s)

        peak_out = max(synth_sel(i / 10.0, sel) for i in range(3600))
        peak_in = min(synth_sel(i / 10.0, sel) for i in range(3600))
        lines += [f"{title}  {mtxt}  ({len(ms)} 本)",
                  f"  f({PARITY_THETA:.0f} 度) = {v0:+.6f} R   f({opp:.0f} 度) = {v1:+.6f} R",
                  f"  円周上の最大 {peak_out:+.6f} R   最小 {peak_in:+.6f} R", ""]
    dmax = max(abs(synth_sel(i / 10.0, PARITY_CASES[0][2]) + synth_sel(i / 10.0, PARITY_CASES[1][2])
                   - synth_sel(i / 10.0, PARITY_CASES[2][2])) for i in range(3600))
    lines.append(f"偶数だけ + 奇数だけ − 両方 の絶対値の最大 (0.1 度刻み): {dmax:.2e}")
    out = "\n".join(lines) + "\n"
    with open(os.path.join(RES_DIR, "fig_cup_parity.txt"), "w", encoding="utf-8") as f:
        f.write(out)
    print(out)


# ---- 開いた系にできない理由 ------------------------------------------------
def harmonic(psi_deg, m, theta=90.0):
    """コップの m 番の倍音。値は R を単位にした長さ。"""
    return A_CUP * weight(m, N_CUP) * math.cos(m * math.radians(psi_deg - theta))


def open_system_figures():
    out_dir = os.path.join(FIG_DIR, "open_system")
    os.makedirs(out_dir, exist_ok=True)
    theta, cut = 90.0, 270.0
    col = lambda m: PURPLE if m % 2 == 1 else ORANGE
    total = lambda d: sum(harmonic(d, m) for m in range(1, N_CUP + 1))
    leg = [
        (BLACK, None, 5, "テーブル (円周、半径 R)"),
        (GREEN, None, 5, f"コップの合成波形 ({WNAME} N = {N_CUP}、一本あたりの振幅 R/{A_CUP_DEN}、中心 {theta:.0f}°)"),
        (PURPLE, "2 3", 2, f"奇数倍音 m = 1, 3, …, {N_CUP if N_CUP % 2 else N_CUP - 1} (一本ずつ)"),
        (ORANGE, "2 3", 2, f"偶数倍音 m = 2, 4, …, {N_CUP - 1 if N_CUP % 2 else N_CUP} (一本ずつ)"),
    ]

    # (1) 円周の上に全ての倍音
    w, h = 900, 800
    c = Circle(450, 440)
    s = svg_open(w, h)
    s.append(text(30, 34, f"コップの倍音を一本ずつ描いた図 (中心 {theta:.0f}°)", size=18))
    s += legend(30, 66, leg)
    s += c.base({theta, cut}, label_r=0.70)
    for m in range(1, N_CUP + 1):
        s.append(c.curve(lambda d, m=m: HARM_MAG * harmonic(d, m), col(m), 0.9, "2 3"))
    s.append(c.curve(total, GREEN, 2.4))
    s.append(text(30, 708, f"倍音: h(ψ) = a R {WFORM}cos m(ψ − θ)、a = 1/{A_CUP_DEN}、N = {N_CUP}、"
                           f"θ = {theta:.0f}°。合成波形 f(ψ) は m = 1…{N_CUP} の和", size=15))
    s.append(text(30, 734, f"描画: 合成波形は r = R + f(ψ)。倍音の点線は小さいので、"
                           f"半径方向に {HARM_MAG} 倍して r = R + {HARM_MAG} h(ψ) で描いた", size=15))
    s.append(text(30, 760, f"{cut:.0f}° では、どの倍音も振れが最大 (奇数は内向き、偶数は外向き)。"
                           f"節になる倍音は一つもない", size=15))
    save(os.path.join(out_dir, "fig_cup_s1_all_harmonics.svg"), s)

    # (2) 270 度で切り開いた直線の図
    w, h = 900, 800
    x0, x1 = 100, 860
    s = svg_open(w, h)
    s.append(text(30, 34, f"{cut:.0f}° で切り開いて直線にした図", size=18))
    s += legend(30, 66, leg[1:])

    def X(u):      # u = 切り口からの角度 (0〜360 度)
        return x0 + (x1 - x0) * u / 360.0

    def panel(y0, y1, vmin, vmax, step, ylabel):
        Y = lambda v: y1 - (y1 - y0) * (v - vmin) / (vmax - vmin)
        out = []
        v = vmin
        while v <= vmax + 1e-12:
            yy = Y(v)
            zero = abs(v) < 1e-12
            out.append(f'<line x1="{x0}" y1="{yy:.1f}" x2="{x1}" y2="{yy:.1f}" stroke="#bbbbbb" '
                       f'stroke-width="{1.4 if zero else 0.6}"/>')
            out.append(text(x0 - 8, yy + 5, "0" if zero else f"{v:+.2f}", size=13, fill=GREY,
                            anchor="end"))
            v += step
        for u in range(0, 361, 45):
            out.append(f'<line x1="{X(u):.1f}" y1="{y1}" x2="{X(u):.1f}" y2="{y1 + 6}" '
                       f'stroke="{GREY}" stroke-width="1"/>')
            out.append(text(X(u), y1 + 22, f"{int((cut + u) % 360)}°", size=13, fill=GREY,
                            anchor="middle"))
        out.append(f'<rect x="{x0}" y="{y0}" width="{x1 - x0}" height="{y1 - y0}" fill="none" '
                   f'stroke="{BLACK}" stroke-width="1.2"/>')
        # 切り口 (両端) を太く
        for xx in (x0, x1):
            out.append(f'<line x1="{xx}" y1="{y0}" x2="{xx}" y2="{y1}" stroke="{BLACK}" '
                       f'stroke-width="3"/>')
        out.append(f'<text x="34" y="{(y0 + y1) / 2:.1f}" font-family="{FONT}" font-size="14" '
                   f'text-anchor="middle" transform="rotate(-90 34 {(y0 + y1) / 2:.1f})">'
                   f'{ylabel}</text>')
        return out, Y

    def poly(func, Y, color, width, dash=None):
        pts = [f"{X(360.0 * i / SAMPLES):.2f},{Y(func(cut + 360.0 * i / SAMPLES)):.2f}"
               for i in range(SAMPLES + 1)]
        d = f' stroke-dasharray="{dash}"' if dash else ""
        return (f'<polyline points="{" ".join(pts)}" fill="none" stroke="{color}" '
                f'stroke-width="{width}" stroke-linejoin="round"{d}/>')

    # 上段: 合成波形
    out, Y = panel(160, 360, -0.05, 0.35, 0.05, "合成波形 (R 単位)")
    s += out
    s.append(poly(total, Y, GREEN, 2.4))
    s.append(text(x0 + 10, 178, "切り口", size=13))
    s.append(text(x1 - 10, 178, "切り口", size=13, anchor="end"))
    # 下段: 倍音を一本ずつ
    out, Y = panel(420, 660, -0.03, 0.03, 0.01, "倍音 一本ずつ (R 単位)")
    s += out
    for m in range(1, N_CUP + 1):
        s.append(poly(lambda d, m=m: harmonic(d, m), Y, col(m), 0.9, "2 3"))
    for m in range(1, N_CUP + 1):      # 切り口での値を点で示す
        for xx in (x0, x1):
            s.append(f'<circle cx="{xx}" cy="{Y(harmonic(cut, m)):.2f}" r="2.6" fill="{col(m)}"/>')
    s.append(text((x0 + x1) / 2, 700, "円周の上の位相 ψ (切り口から一周)", size=14, anchor="middle"))
    s.append(text(30, 734, f"切り口 ({cut:.0f}°) では、どの倍音も 0 にならない。値は奇数倍音が −a R w、"
                           f"偶数倍音が +a R w で、振れが最大の位置にある", size=15))
    s.append(text(30, 760, f"合成波形は切り口で {total(cut):+.4f} R と小さいが、"
                           f"それは奇数と偶数が打ち消し合った結果であって、節ではない", size=15))
    save(os.path.join(out_dir, "fig_cup_cut_open_270.svg"), s)

    # 数値
    lines = [f"コップの倍音 (中心 {theta:.0f} 度) の、切り口 {cut:.0f} 度での値と傾き。R を単位にした値",
             "  m   振幅 a w_m     切り口での値    値/振幅   切り口での傾き d/dpsi (ラジアンあたり)"]
    for m in range(1, N_CUP + 1):
        amp = A_CUP * weight(m, N_CUP)
        val = harmonic(cut, m)
        slope = -amp * m * math.sin(m * math.radians(cut - theta))
        lines.append(f"{m:3d}  {amp:.8f}   {val:+.8f}   {val / amp:+.3f}   {slope:+.2e}")
    lines += ["", f"合成波形の切り口での値: {total(cut):+.8f} R"]
    # 全ての倍音が同時に節になる位置があるか (0.01 度刻み)
    best = min(((max(abs(math.cos(m * math.radians(i / 100.0 - theta))) for m in range(1, N_CUP + 1)),
                 i / 100.0) for i in range(36000)), key=lambda t: t[0])
    lines += ["",
              "全ての倍音が同時に 0 になる位置を探した (0.01 度刻み)。",
              "各位置で max_m |cos m(psi − theta)| を取り、その最小を見る。0 なら共通の節がある。",
              f"  最小値 {best[0]:.6f}  (位置 {best[1]:.2f} 度)  → 共通の節は無い"]
    out = "\n".join(lines) + "\n"
    with open(os.path.join(RES_DIR, "fig_cup_open_system.txt"), "w", encoding="utf-8") as f:
        f.write(out)
    print(out)


# ---- 両端を固定した直線の上の定常波 ----------------------------------------
def line_mode(u, k, fejer=False):
    """u = x/L (0〜1)。k 番の定常波。値は R を単位にした長さ。"""
    wk = 1.0 - k / (N_LINE + 1.0) if fejer else 1.0
    return wk * math.sin(k * math.pi * u) / A_LINE_DEN


def line_tail(u):
    """裾の曲線 (a R / 2) cot(π x / (2 L))。u = x/L。"""
    return 0.5 / A_LINE_DEN / math.tan(0.5 * math.pi * u)


def line_fixed_ends_figure(fejer=False):
    out_dir = os.path.join(FIG_DIR, "open_system")
    os.makedirs(out_dir, exist_ok=True)
    total = lambda u: sum(line_mode(u, k, fejer) for k in range(1, N_LINE + 1))
    total_eq = lambda u: sum(line_mode(u, k, False) for k in range(1, N_LINE + 1))
    wname = "フェイエール重み" if fejer else "等振幅"
    wform = "(1 − k/(N+1)) " if fejer else ""
    tag = "fejer" if fejer else "equal"
    col = lambda k: PURPLE if k % 2 == 1 else ORANGE
    n = 7200
    vals = [total(i / n) for i in range(n + 1)]
    imax = max(range(n + 1), key=lambda i: vals[i])
    imin = min(range(n + 1), key=lambda i: vals[i])

    w, h = 900, 810
    x0, x1 = 100, 860
    X = lambda u: x0 + (x1 - x0) * u
    s = svg_open(w, h)
    s.append(text(30, 34, f"両端を固定した長さ L = 2πR の直線に、定常波 1〜{N_LINE} 番を足した波形 ({wname})",
                  size=18))
    leg = [(GREEN, None, 5, f"合成波形 (k = 1〜{N_LINE}、{wname}、一本あたりの振幅 R/{A_LINE_DEN})")]
    if fejer:
        leg += [("#999999", "7 4", 2, "等振幅で足した波形 (比較)"),
                (BLACK, "2 3", 2, "裾の曲線 (a R/2) cot(πx/(2L))")]
    leg += [(PURPLE, "2 3", 2, "奇数番の定常波 (一本ずつ)"),
            (ORANGE, "2 3", 2, "偶数番の定常波 (一本ずつ)")]
    s += legend(30, 62, leg, size=14)

    def panel(y0, y1, vmin, vmax, step, ylabel, fmt):
        Y = lambda v: y1 - (y1 - y0) * (v - vmin) / (vmax - vmin)
        out = []
        nstep = int(round((vmax - vmin) / step))
        for j in range(nstep + 1):
            v = vmin + j * step
            zero = abs(v) < 1e-12
            out.append(f'<line x1="{x0}" y1="{Y(v):.1f}" x2="{x1}" y2="{Y(v):.1f}" stroke="#bbbbbb" '
                       f'stroke-width="{1.4 if zero else 0.6}"/>')
            out.append(text(x0 - 8, Y(v) + 5, "0" if zero else fmt.format(v), size=13, fill=GREY,
                            anchor="end"))
        for j, lab in enumerate(["0", "π/2", "π", "3π/2", "2π"]):
            xx = X(j / 4.0)
            out.append(f'<line x1="{xx:.1f}" y1="{y1}" x2="{xx:.1f}" y2="{y1 + 6}" stroke="{GREY}" '
                       f'stroke-width="1"/>')
            out.append(text(xx, y1 + 22, lab, size=13, fill=GREY, anchor="middle"))
        out.append(f'<rect x="{x0}" y="{y0}" width="{x1 - x0}" height="{y1 - y0}" fill="none" '
                   f'stroke="{BLACK}" stroke-width="1.2"/>')
        for xx in (x0, x1):
            out.append(f'<line x1="{xx}" y1="{y0}" x2="{xx}" y2="{y1}" stroke="{BLACK}" '
                       f'stroke-width="3"/>')
        out.append(f'<text x="34" y="{(y0 + y1) / 2:.1f}" font-family="{FONT}" font-size="14" '
                   f'text-anchor="middle" transform="rotate(-90 34 {(y0 + y1) / 2:.1f})">'
                   f'{ylabel}</text>')
        return out, Y

    def poly(func, Y, color, width, dash=None):
        pts = [f"{X(i / SAMPLES):.2f},{Y(func(i / SAMPLES)):.2f}" for i in range(SAMPLES + 1)]
        d = f' stroke-dasharray="{dash}"' if dash else ""
        return (f'<polyline points="{" ".join(pts)}" fill="none" stroke="{color}" '
                f'stroke-width="{width}" stroke-linejoin="round"{d}/>')

    out, Y = panel(190, 400, -0.10, 0.35, 0.05, "合成波形 (R 単位)", "{:+.2f}")
    s += out
    if fejer:
        s.append(poly(total_eq, Y, "#999999", 1.4, "7 4"))
        # 裾の曲線は端で発散するので、枠の中に入る範囲だけ描く
        pts = [f"{X(i / SAMPLES):.2f},{Y(line_tail(i / SAMPLES)):.2f}"
               for i in range(1, SAMPLES + 1) if line_tail(i / SAMPLES) <= 0.35]
        s.append(f'<polyline points="{" ".join(pts)}" fill="none" stroke="{BLACK}" '
                 f'stroke-width="1.4" stroke-dasharray="2 3"/>')
    s.append(poly(total, Y, GREEN, 2.4))
    s.append(text(x0 + 10, 390, "固定端", size=13))
    s.append(text(x1 - 10, 390, "固定端", size=13, anchor="end"))
    out, Y = panel(450, 650, -0.03, 0.03, 0.01, "定常波 一本ずつ (R 単位)", "{:+.2f}")
    s += out
    for k in range(1, N_LINE + 1):
        s.append(poly(lambda u, k=k: line_mode(u, k, fejer), Y, col(k), 0.9, "2 3"))
    s.append(text((x0 + x1) / 2, 690, "直線の上の位置 x/R", size=14, anchor="middle"))

    s.append(text(30, 722, f"波形: g(x) = a R Σ{SUM_SUB.replace('m=1', 'k=1')}{wform}sin(k π x / L)、L = 2πR、a = 1/{A_LINE_DEN}、"
                           f"N = {N_LINE}。両端では全ての定常波が 0", size=15))
    s.append(text(30, 748, f"山は直線の端に寄る: 最大 {vals[imax]:+.4f} R が x/R = "
                           f"{2 * math.pi * imax / n:.4f} (端から全長の {100.0 * imax / n:.2f}%)", size=15))
    if fejer:
        s.append(text(30, 774, "重みを掛けると波打ちは消えるが、正の側に片寄った長い裾 (黒の点線) は残る。"
                               f"中央での値 {vals[n // 2]:+.4f} R", size=15))
    else:
        s.append(text(30, 774, f"山の先は、正の側に片寄ったまま波打って小さくなり、反対の端で 0 になる。"
                               f"最小は {vals[imin]:+.4f} R (x/R = {2 * math.pi * imin / n:.4f})", size=15))
    save(os.path.join(out_dir, f"fig_line_fixed_ends_{tag}{N_LINE}.svg"), s)

    lines = [f"両端固定、長さ L = 2πR の直線。g(x) = a R sum_{{k=1}}^{{{N_LINE}}} w_k sin(k π x / L)、a = 1/{A_LINE_DEN}",
             f"重み: {wname}",
             "値は R を単位にした長さ", "",
             f"最大 {vals[imax]:+.6f} R   位置 x/L = {imax / n:.5f}  (x/R = {2 * math.pi * imax / n:.5f})",
             f"最小 {vals[imin]:+.6f} R   位置 x/L = {imin / n:.5f}  (x/R = {2 * math.pi * imin / n:.5f})",
             f"両端での値 {vals[0]:+.2e}, {vals[n]:+.2e}",
             f"中央 x = L/2 での値 {vals[n // 2]:+.6f} R",
             "", "x/L      g/R        裾の曲線 (a/2)cot(πx/2L)   差"]
    lines += [f"{i / 40:.3f}  {total(i / 40):+.6f}   {line_tail(i / 40):+.6f}   "
              f"{total(i / 40) - line_tail(i / 40):+.6f}" for i in range(1, 41)]
    out = "\n".join(lines) + "\n"
    with open(os.path.join(RES_DIR, f"fig_line_fixed_ends_{tag}.txt"), "w", encoding="utf-8") as f:
        f.write(out)
    print(out)


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    os.makedirs(RES_DIR, exist_ok=True)

    schematic(os.path.join(FIG_DIR, "fig_cup_coin_schematic_covered.svg"),
              cup_center_x=250, title="コップをコインの上に被せた状態")
    schematic(os.path.join(FIG_DIR, "fig_cup_coin_schematic_beside.svg"),
              cup_center_x=460, title="コップをコインの隣に置いた状態")
    for th in THETA_CUP_LIST:
        s1_figure(os.path.join(FIG_DIR, f"fig_cup_coin_s1_phase{th:.0f}.svg"), th)

    lines = ["コインとコップの合成波形 f(psi) = a R * sum_{m=1}^{N} w_m cos(m (psi - theta))",
             f"重み: {WEIGHT}" + ("  w_m = 1 - m/(N+1)" if WEIGHT == "fejer" else "  w_m = 1"),
             "長さはすべて R を単位にした値", ""]
    for name, n, amp, den in (("コイン", N_COIN, A_COIN, A_COIN_DEN),
                              ("コップ", N_CUP, A_CUP, A_CUP_DEN)):
        peak, half, omax, omin = profile(n, amp)
        lines += [
            f"{name}  N = {n}, a = 1/{den}",
            f"  山の高さ                        : {peak:.6f} R",
            f"  山の半幅 (中心から最初の極小)   : {half:.3f} 度",
            f"  山の幅 (底から底)               : {2 * half:.3f} 度",
            f"  山の外側の最大と最小            : {omax:+.6f} R, {omin:+.6f} R",
            f"  山の外側の振れ幅 ÷ 山の高さ     : {(omax - omin) / peak:.4f}",
            f"  正反対 (中心から 180 度) の値   : {synth(180.0, n, 0.0, amp):+.6f} R",
            "",
        ]
    lines += ["振幅の分解能 (仮説): k × 山の高さ ÷ N   (k = %g)" % KAPPA]
    for name, n, amp, den in (("コイン", N_COIN, A_COIN, A_COIN_DEN),
                              ("コップ", N_CUP, A_CUP, A_CUP_DEN)):
        peak, _, omax, omin = profile(n, amp)
        eps = amp_resolution_wave(n, amp)
        lines.append(f"  {name}: {eps:.6f} R  (一本あたりの振幅 a R = {amp:.6f} R)   "
                     f"山の外の最大・最小 ÷ 分解能 = {omax / eps:+.3f}, {omin / eps:+.3f}")
    lines.append("")
    lines += ["コップの中心位相ごとの、コインの中心 (90 度) でのコップの値"]
    for th in THETA_CUP_LIST:
        lines.append(f"  コップ中心 {th:.0f} 度 : {f_cup(THETA_COIN, th):+.6f} R")
    out = "\n".join(lines) + "\n"
    with open(os.path.join(RES_DIR, "fig_cup_coin_table.txt"), "w", encoding="utf-8") as f:
        f.write(out)
    print(out)
    cross_figures()
    cross_summary_figure(os.path.join(FIG_DIR, "fig_cup_coin_cross_summary.svg"))
    parity_figures()
    open_system_figures()
    line_fixed_ends_figure(fejer=False)
    line_fixed_ends_figure(fejer=True)


if __name__ == "__main__":
    main()
