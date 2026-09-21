#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_coulomb_readout_figures_v1.py

第二思考実験の図 fig01〜fig05 を、日本語版・英語版、SVG・PNG で作る。

  * 図は一度だけ記述し（Canvas に線・点・文字を置く）、同じ記述から SVG と PNG を出す。
    SVG は依存なしの簡易出力（文字はテキストのまま）。PNG は matplotlib がある場合だけ。
  * fig01〜fig03 は数値実験の内容（軌道は決定的に再計算。収束の表は
    coulomb_readout_experiments_results_v1.json から読む）。fig04・fig05 は概念図。
  * 出力はこのスクリプトと同じフォルダ。名前は figNN_coulomb_readout_<内容>_<ja|en>_v1.<svg|png>

依存: numpy。PNG も作る場合は matplotlib。  実行: python3 make_coulomb_readout_figures_v1.py
"""
import json
import math
import os
from pathlib import Path

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager as fm
    HAVE_MPL = True
except ImportError:
    HAVE_MPL = False

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "coulomb_readout_experiments_results_v1.json"
C = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#7f7f7f"]
JP_FONTS = ["/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc", "/usr/share/fonts/truetype/noto/NotoSansCJK-Regular.ttc",
            "/System/Library/Fonts/ヒラギノ角ゴシック W3.ttc", "/Library/Fonts/Arial Unicode.ttf",
            "C:/Windows/Fonts/YuGothR.ttc", "C:/Windows/Fonts/meiryo.ttc", "C:/Windows/Fonts/msgothic.ttc"]


def _n(v):
    s = f"{v:.1f}"
    return s[:-2] if s.endswith(".0") else s


def _esc(t):
    return t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


class Canvas:
    """図の記述。座標はピクセル（左上が原点）。"""

    def __init__(self, w, h):
        self.w, self.h, self.ops = w, h, []

    def line(self, pts, color="#000", width=1.5, dash=None):
        self.ops.append(("line", [(float(x), float(y)) for x, y in pts], color, width, dash))

    def segs(self, segments, color="#ddd", width=1):
        self.ops.append(("segs", [tuple(map(float, sg)) for sg in segments], color, width))

    def dots(self, pts, color, r=2.5):
        self.ops.append(("dots", [(float(x), float(y)) for x, y in pts], color, r))

    def poly(self, pts, fill, opacity=0.25):
        self.ops.append(("poly", [(float(x), float(y)) for x, y in pts], fill, opacity))

    def rect(self, x, y, w, h, fill="none", stroke="#000", opacity=1.0):
        self.ops.append(("rect", x, y, w, h, fill, stroke, opacity))

    def text(self, x, y, s, size=12, anchor="start", rot=0, color="#000"):
        self.ops.append(("text", float(x), float(y), s, size, anchor, rot, color))

    def arrow(self, x0, y0, x1, y1, color="#000", width=1.5, head=8):
        a = math.atan2(y1 - y0, x1 - x0)
        self.line([(x0, y0), (x1, y1)], color, width)
        self.poly([(x1, y1), (x1 - head * math.cos(a - 0.4), y1 - head * math.sin(a - 0.4)),
                   (x1 - head * math.cos(a + 0.4), y1 - head * math.sin(a + 0.4))], color, 1.0)

    # ---- SVG ----
    @staticmethod
    def _body(s, size):
        if "^" in s:                                         # "10^−5" → 右肩の指数
            b, e = s.split("^", 1)
            return f'{_esc(b)}<tspan dy="-5" font-size="{size - 3}">{_esc(e)}</tspan>'
        return _esc(s)

    def to_svg(self, path):
        o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {self.w} {self.h}" '
             f'font-family="Hiragino Sans, Yu Gothic, Noto Sans CJK JP, sans-serif">',
             f'<rect width="{self.w}" height="{self.h}" fill="#fff"/>']
        i = 0
        while i < len(self.ops):
            op = self.ops[i]
            i += 1
            k = op[0]
            if k == "text" and not op[6]:                    # 同じ書式の文字が続くときは <g> にまとめる
                run = [op]
                while i < len(self.ops) and self.ops[i][0] == "text" and not self.ops[i][6] and self.ops[i][4:6] == op[4:6] and self.ops[i][7] == op[7]:
                    run.append(self.ops[i])
                    i += 1
                if len(run) > 2:
                    fc = "" if op[7] == "#000" else f' fill="{op[7]}"'
                    o.append(f'<g font-size="{op[4]}" text-anchor="{op[5]}"{fc}>' + "".join(
                        f'<text x="{_n(r[1])}" y="{_n(r[2])}">{self._body(r[3], r[4])}</text>' for r in run) + "</g>")
                    continue
                i -= len(run) - 1
            if k == "line":
                d = " ".join(f"{_n(x)},{_n(y)}" for x, y in op[1])
                da = f' stroke-dasharray="{op[4]}"' if op[4] else ""
                o.append(f'<polyline points="{d}" fill="none" stroke="{op[2]}" stroke-width="{op[3]}"{da}/>')
            elif k == "segs":
                d = "".join(f"M{_n(a)} {_n(b)}L{_n(c)} {_n(e)}" for a, b, c, e in op[1])
                o.append(f'<path d="{d}" stroke="{op[2]}" stroke-width="{op[3]}" fill="none"/>')
            elif k == "dots":
                d = "".join(f"M{_n(x)} {_n(y)}h0" for x, y in op[1])
                o.append(f'<path d="{d}" stroke="{op[2]}" stroke-width="{2 * op[3]}" stroke-linecap="round" fill="none"/>')
            elif k == "poly":
                d = " ".join(f"{_n(x)},{_n(y)}" for x, y in op[1])
                o.append(f'<polygon points="{d}" fill="{op[2]}" fill-opacity="{op[3]}"/>')
            elif k == "rect":
                _, x, y, w, h, fill, stroke, op_ = op
                o.append(f'<rect x="{_n(x)}" y="{_n(y)}" width="{_n(w)}" height="{_n(h)}" fill="{fill}" '
                         f'fill-opacity="{op_}" stroke="{stroke}"/>')
            else:
                _, x, y, s, size, anchor, rot, color = op
                tr = f' transform="rotate({-rot} {_n(x)} {_n(y)})"' if rot else ""
                fc = "" if color == "#000" else f' fill="{color}"'
                o.append(f'<text x="{_n(x)}" y="{_n(y)}" font-size="{size}" text-anchor="{anchor}"{tr}{fc}>{self._body(s, size)}</text>')
        with open(path, "w", encoding="utf-8", newline="\n") as f:
            f.write("\n".join(o) + "\n</svg>\n")

    # ---- PNG（matplotlib を描画器として使う。配置は SVG と同じ） ----
    def to_png(self, path, fp):
        k = 0.72                                        # px → pt（100 px = 1 inch）
        fig = plt.figure(figsize=(self.w / 100, self.h / 100), dpi=200)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_xlim(0, self.w); ax.set_ylim(self.h, 0); ax.axis("off")
        for op in self.ops:
            t = op[0]
            if t == "line":
                xs, ys = zip(*op[1])
                ls = (0, tuple(float(v) * k for v in op[4].split())) if op[4] else "-"
                ax.plot(xs, ys, color=op[2], lw=op[3] * k, ls=ls, solid_capstyle="butt")
            elif t == "segs":
                for x0, y0, x1, y1 in op[1]:
                    ax.plot([x0, x1], [y0, y1], color=op[2], lw=op[3] * k)
            elif t == "dots":
                xs, ys = zip(*op[1])
                ax.plot(xs, ys, "o", color=op[2], ms=2 * op[3] * k, mew=0)
            elif t == "poly":
                xs, ys = zip(*op[1])
                ax.fill(xs, ys, color=op[2], alpha=op[3], lw=0)
            elif t == "rect":
                _, x, y, w, h, fill, stroke, op_ = op
                ax.add_patch(plt.Rectangle((x, y), w, h, facecolor=fill if fill != "none" else "none",
                                           edgecolor=stroke if stroke != "none" else "none", alpha=op_, lw=0.8 * k))
            else:
                _, x, y, s, size, anchor, rot, color = op
                if "^" in s:
                    b, e = s.split("^", 1)
                    s = f"${b}^{{{e.replace('−', '-')}}}$"
                ha = {"start": "left", "middle": "center", "end": "right"}[anchor]
                kw = {"fontproperties": fp} if fp is not None else {}
                ax.text(x, y, s, fontsize=size * k, ha=ha, va="baseline", rotation=rot, color=color,
                        rotation_mode="anchor", **kw)
        fig.savefig(path)
        plt.close(fig)


class Axes:
    def __init__(self, cv, x0, y0, w, h, xlim, ylim, xlog=False, ylog=False):
        self.cv, self.x0, self.y0, self.w, self.h, self.xlog, self.ylog = cv, x0, y0, w, h, xlog, ylog
        self.xl = [math.log10(v) for v in xlim] if xlog else list(xlim)
        self.yl = [math.log10(v) for v in ylim] if ylog else list(ylim)

    def X(self, v):
        v = math.log10(v) if self.xlog else v
        return self.x0 + (v - self.xl[0]) / (self.xl[1] - self.xl[0]) * self.w

    def Y(self, v):
        v = math.log10(v) if self.ylog else v
        return self.y0 + self.h - (v - self.yl[0]) / (self.yl[1] - self.yl[0]) * self.h

    def pts(self, xs, ys):
        return [(self.X(x), self.Y(y)) for x, y in zip(xs, ys)]

    def frame(self, xt, yt, title="", xlabel="", ylabel="", xfmt=None, yfmt=None, dx=40):
        cv = self.cv
        grid = [(self.X(v), self.y0, self.X(v), self.y0 + self.h) for v in xt] + \
               [(self.x0, self.Y(v), self.x0 + self.w, self.Y(v)) for v in yt]
        if grid:
            cv.segs(grid)
        for v in xt:
            cv.text(self.X(v), self.y0 + self.h + 15, (xfmt or _n)(v), 11, "middle")
        for v in yt:
            cv.text(self.x0 - 6, self.Y(v) + 4, (yfmt or _n)(v), 11, "end")
        cv.rect(self.x0, self.y0, self.w, self.h)
        if title:
            cv.text(self.x0 + self.w / 2, self.y0 - 9, title, 13, "middle")
        if xlabel:
            cv.text(self.x0 + self.w / 2, self.y0 + self.h + 33, xlabel, 12, "middle")
        if ylabel:
            cv.text(self.x0 - dx, self.y0 + self.h / 2, ylabel, 12, "middle", rot=90)

    def legend(self, entries, x, y, box_w):
        self.cv.rect(x - 6, y - 11, box_w, 16 * len(entries) + 6, "#fff", "#ccc", 0.85)
        for i, (label, color, dash) in enumerate(entries):
            self.cv.line([(x, y + 16 * i), (x + 22, y + 16 * i)], color, 2, dash)
            self.cv.text(x + 28, y + 16 * i + 4, label, 11)


def _p10(v):
    e = int(round(math.log10(v)))
    return "1" if e == 0 else ("10" if e == 1 else "10^" + str(e).replace("-", "−"))


# --------------------------------------------------------------------------- 軌道（決定的）
def elliptic_orbit(n, e, psi=0.5, u0=0.3):
    th = 2 * math.pi / n
    R = lambda a: np.array([[math.cos(a), -math.sin(a)], [math.sin(a), math.cos(a)]])
    T = R(psi) @ np.diag([math.sqrt(1 - e), math.sqrt(1 + e)])
    S = T @ R(th) @ np.linalg.inv(T)
    X = [T @ np.array([math.cos(u0), math.sin(u0)])]
    for _ in range(n + 1):
        X.append(S @ X[-1])
    u = np.linspace(0, 2 * math.pi, 97)
    return S, np.array(X), (T @ np.vstack([np.cos(u), np.sin(u)])).T, th


def hyperbolic_orbits(H=0.25, a1=math.radians(10), a2=math.radians(70)):
    vp, vm = np.array([math.cos(a1), math.sin(a1)]), np.array([math.cos(a2), math.sin(a2)])
    out = {}
    for name, (al, be) in {"rep": (0.6, 0.6), "att": (0.6, -0.6)}.items():
        out[name] = np.array([al * math.exp(k * H) * vp + be * math.exp(-k * H) * vm for k in range(-7, 8)])
    return vp, vm, out


def qmap(X):
    """実二次写像 Phi(a,b)=(a^2-b^2,2ab)。"""
    X = np.asarray(X, dtype=float)
    return np.column_stack((X[:, 0] ** 2 - X[:, 1] ** 2, 2 * X[:, 0] * X[:, 1]))


def sq(X):
    Y = qmap(X)
    return Y[:, 0], Y[:, 1]


def real_hyperbolic_eigenvectors(S):
    tau = float(np.trace(S))
    disc = math.sqrt(max(0.0, tau * tau - 4.0))
    vals = [(tau + disc) / 2.0, (tau - disc) / 2.0]
    vecs = []
    for lam in vals:
        if abs(S[0, 1]) + abs(lam - S[0, 0]) >= abs(lam - S[1, 1]) + abs(S[1, 0]):
            v = np.array([S[0, 1], lam - S[0, 0]], dtype=float)
        else:
            v = np.array([lam - S[1, 1], S[1, 0]], dtype=float)
        v /= np.linalg.norm(v)
        vecs.append(v)
    order = np.argsort([-abs(vals[0]), -abs(vals[1])])
    return vals[order[0]], vecs[order[0]], vals[order[1]], vecs[order[1]]


# --------------------------------------------------------------------------- 図
L = {
 "ja": {
  "1a": "(a) 二つの値の平面：線形の法則の軌道", "1b": "(b) 読み出し平面：同じ軌道。原点が焦点になる",
  "1c": "(c) 双曲型の法則：二つの区画の軌道", "1d": "(d) 読み出し平面：引力と斥力の双曲線",
  "focus": "原点＝焦点（力の中心）", "att": "引力（μ>0）", "rep": "斥力（μ<0）", "null": "零方向 v+, v−",
  "2a": "(a) 時計の刻み Δt_k と距離 r_k（一周 48 歩、e=0.6）", "2b": "(b) 差分加速度 × r² ÷ μ（1 なら逆二乗）",
  "2c": "(c) 逆二乗則からのずれ：歩幅 θ の 2 次で減る", "2d": "(d) 時計がどの初期位相でも単調になる範囲",
  "k": "歩数 k", "dt": "Δt_k（平均で規格化）", "r": "r_k（平均で規格化）", "n64": "一周 64 歩", "n512": "一周 512 歩",
  "theta": "歩幅 θ", "err": "最大の相対誤差", "slope2": "傾き 2", "n": "一周の歩数 n", "e": "離心率 e",
  "mono": "e < cos(2π/n)：単調", "alias": "角度の折り返しが起こりうる",
  "3t": "結合の符号は、配置が入っている区画の開き角 ψ で決まる：μ = −8αβ cos ψ", "3x": "配置が入っている区画（二つの零方向のあいだ）の開き角 ψ（度）",
  "3y": "μ / (8αβ)", "3c": "−cos ψ", "3s": "乱数で選んだ法則", "3n": "直角（厳密に無名な法則を含む）：中性、μ=0",
  "4t": "三つの面積読み出しと、モデルの流れ", "4ab": "二つの値の平面", "4x": "x = ω(X, PX) = a²−b²", "4y": "y = ω(X, DX) = 2ab",
  "4r": "r = ω(X, J0X) = a²+b²", "4P": "PX（入れ替え）", "4D": "DX（符号反転）", "4J": "J0X（四分の一回転）", "4z": "x² + y² − r² = 0（ゼロ閉塞）", "4area": "面積 = r",
  "4f": ["二つの値 X=(a,b)", "線形の法則  X_{k+1}+X_{k−1}=τX_k", "位置＝三つの面積読み出し (x, y, r)", "時計  Δt_k ∝ X_k·X_{k+1}（面積の保存から）", "逆二乗の力・焦点が原点の円錐曲線"],
  "5t": "二つの値で読めるもの、三体で初めて読めるもの", "5a": "二体：関係性は一つ", "5b": "三体：三角形",
  "5a1": "読める：法則の数 τ、離心率 e、", "5a2": "結合 μ の符号、角度差と距離の比", "5a3": "読めない：質量比、電荷比、地図",
  "5b1": "形の自由度 2N−4：N=2 で 0、N=3 で 2", "5b2": "s12·s23·s31 = +1 ⇔ sij = σi·σj", "5b3": "（固有の電荷は関係性の整合として立つ）"},
 "en": {
  "1a": "(a) Two-value plane: orbit of the linear law", "1b": "(b) Readout plane: same orbit, focus at the origin",
  "1c": "(c) Hyperbolic law: orbits in two sectors", "1d": "(d) Readout plane: attractive and repulsive",
  "focus": "origin = focus (center of force)", "att": "attractive (μ>0)", "rep": "repulsive (μ<0)", "null": "null directions v+, v−",
  "2a": "(a) Clock tick and distance (48 steps per turn, e=0.6)", "2b": "(b) Difference acceleration × r² ÷ μ",
  "2c": "(c) Deviation from inverse-square: order 2 in θ", "2d": "(d) Where the clock is monotone for every phase",
  "k": "step k", "dt": "Δt_k (normalized by the mean)", "r": "r_k (normalized by the mean)", "n64": "64 steps per turn", "n512": "512 steps per turn",
  "theta": "step θ", "err": "maximum relative error", "slope2": "slope 2", "n": "steps per turn n", "e": "eccentricity e",
  "mono": "e < cos(2π/n): monotone", "alias": "angle wrap-around can occur",
  "3t": "Sign of the coupling: μ = −8αβ cos ψ (ψ: opening angle of the sector containing the configuration)",
  "3x": "opening angle ψ of the sector containing the configuration (degrees)",
  "3y": "μ / (8αβ)", "3c": "−cos ψ", "3s": "randomly chosen laws", "3n": "right angle (includes the strictly nameless law): neutral, μ=0",
  "4t": "The three area readouts and the flow of the model", "4ab": "plane of the two values", "4x": "x = ω(X, PX) = a²−b²", "4y": "y = ω(X, DX) = 2ab",
  "4r": "r = ω(X, J0X) = a²+b²", "4P": "PX (swap)", "4D": "DX (sign flip)", "4J": "J0X (quarter turn)", "4z": "x² + y² − r² = 0 (zero closure)", "4area": "area = r",
  "4f": ["two values X=(a,b)", "linear law  X_{k+1}+X_{k−1}=τX_k", "position = three area readouts (x, y, r)", "clock  Δt_k ∝ X_k·X_{k+1} (from area conservation)", "inverse-square force; conics with a focus at the origin"],
  "5t": "What two values can read, and what first becomes readable with three bodies", "5a": "two bodies: one relation", "5b": "three bodies: a triangle",
  "5a1": "readable: the law's number τ, eccentricity e,", "5a2": "sign of the coupling μ, angle steps and distance ratios", "5a3": "not readable: mass ratio, charge ratio, a map",
  "5b1": "shape degrees of freedom 2N−4: 0 for N=2, 2 for N=3", "5b2": "s12·s23·s31 = +1 ⇔ sij = σi·σj", "5b3": "(intrinsic charge arises as consistency of relations)"},
}


def fig1(t):
    cv = Canvas(920, 840)
    S, X, ell, th = elliptic_orbit(24, 0.6)
    a = Axes(cv, 60, 40, 380, 320, (-1.6, 1.6), (-1.35, 1.35))
    a.frame([-1, 0, 1], [-1, 0, 1], t["1a"], "a", "b")
    cv.line(a.pts(ell[::2, 0], ell[::2, 1]), C[1], 1)
    cv.dots(a.pts(X[:24, 0], X[:24, 1]), C[0])
    cv.dots([(a.X(0), a.Y(0))], "#000", 3)
    b = Axes(cv, 520, 40, 380, 320, (-2.0, 1.8), (-1.6, 1.6))
    b.frame([-2, -1, 0, 1], [-1, 0, 1], t["1b"], "x", "y")
    ex, ey = sq(ell[:49])                                # X の半周が読み出しの一周
    zx, zy = sq(X[:12])
    cv.line(b.pts(ex, ey), C[1], 1)
    cv.dots(b.pts(zx, zy), C[0])
    cv.dots([(b.X(0), b.Y(0))], "#000", 3.5)
    cv.text(b.X(0) + 8, b.Y(0) - 8, t["focus"], 11)
    vp, vm, orb = hyperbolic_orbits()
    c = Axes(cv, 60, 470, 380, 320, (-4.2, 4.2), (-3.5, 3.5))
    c.frame([-4, -2, 0, 2, 4], [-2, 0, 2], t["1c"], "a", "b")
    for v in (vp, vm):
        m_ = min(4.2 / abs(v[0]), 3.5 / abs(v[1]))          # 枠の中に収める
        cv.line(c.pts([-m_ * v[0], m_ * v[0]], [-m_ * v[1], m_ * v[1]]), "#000", 1, "4 3")
    for name, col in (("att", C[2]), ("rep", C[3])):
        cv.line(c.pts(orb[name][:, 0], orb[name][:, 1]), col, 1)
        cv.dots(c.pts(orb[name][:, 0], orb[name][:, 1]), col)
    c.legend([(t["att"], C[2], None), (t["rep"], C[3], None), (t["null"], "#000", "4 3")], 72, 486, 150)
    d = Axes(cv, 520, 470, 380, 320, (-9, 9), (-7.5, 7.5))
    d.frame([-8, -4, 0, 4, 8], [-4, 0, 4], t["1d"], "x", "y")
    for name, col in (("att", C[2]), ("rep", C[3])):
        zx, zy = sq(orb[name])
        keep = (np.abs(zx) < 9) & (np.abs(zy) < 7.5)
        cv.line(d.pts(zx[keep], zy[keep]), col, 1)
        cv.dots(d.pts(zx[keep], zy[keep]), col)
    cv.dots([(d.X(0), d.Y(0))], "#000", 3.5)
    cv.text(d.X(0) + 10, d.Y(0) + 34, t["focus"], 11)
    return cv


def clock_series(n, e):
    S, X, _, th = elliptic_orbit(n, e)
    dt = np.array([(th / math.cos(th)) * float(X[k] @ X[k + 1]) for k in range(n + 1)])
    Y = qmap(X)
    A = (S - np.linalg.inv(S)) / (2 * math.sin(th))
    dX = A @ X[0]
    mu = 2 * float(dX @ dX) + 2 * float(X[0] @ X[0])
    acc = [2 * ((Y[k + 1] - Y[k]) / dt[k] - (Y[k] - Y[k - 1]) / dt[k - 1]) / (dt[k] + dt[k - 1]) for k in range(1, n)]
    radius = np.linalg.norm(Y, axis=1)
    ratio = [np.linalg.norm(acc[k - 1]) * radius[k] ** 2 / mu for k in range(1, n)]
    return dt[:n], radius[:n], np.array(ratio)


def fig2(t):
    cv = Canvas(920, 840)
    dt, r, _ = clock_series(48, 0.6)
    a = Axes(cv, 60, 40, 380, 320, (0, 48), (0, 2.0))
    a.frame([0, 12, 24, 36, 48], [0, 0.5, 1, 1.5, 2], t["2a"], t["k"])
    cv.line(a.pts(range(48), dt / dt.mean()), C[0], 1.6)
    cv.line(a.pts(range(48), r / r.mean()), C[1], 1.6, "5 3")
    a.legend([(t["dt"], C[0], None), (t["r"], C[1], "5 3")], 72, 56, 200)
    b = Axes(cv, 520, 40, 380, 320, (0, 1), (0.96, 1.04))
    b.frame([0, 0.5, 1], [0.96, 0.98, 1.0, 1.02, 1.04], t["2b"], "k / n", yfmt=lambda v: f"{v:.2f}")
    for n, col, key in ((64, C[0], "n64"), (512, C[3], "n512")):
        rt = clock_series(n, 0.6)[2]
        kk = np.arange(1, n) / n
        step = max(1, n // 32)
        cv.line(b.pts(kk[::step], rt[::step]), col, 1.4)
    b.legend([(t["n64"], C[0], None), (t["n512"], C[3], None)], 532, 56, 150)
    res = json.load(open(RESULTS, encoding="utf-8"))["E4"]
    ths = [2 * math.pi / n for n in res["E4d_n"]]
    c = Axes(cv, 60, 470, 380, 320, (4e-3, 0.4), (3e-6, 1e-1), True, True)
    c.frame([1e-2, 1e-1], [1e-5, 1e-4, 1e-3, 1e-2, 1e-1], t["2c"], t["theta"], t["err"], _p10, _p10, 46)
    cv.line(c.pts(ths, res["E4d_errors"]), C[0], 1.6)
    cv.dots(c.pts(ths, res["E4d_errors"]), C[0], 3.5)
    cv.line(c.pts([6e-3, 0.3], [4e-6, 4e-6 * (0.3 / 6e-3) ** 2]), "#000", 1, "5 4")
    cv.text(c.X(0.05), c.Y(2e-4), t["slope2"], 11)
    d = Axes(cv, 520, 470, 380, 320, (4, 40), (0, 1))
    d.frame([5, 10, 20, 30, 40], [0, 0.25, 0.5, 0.75, 1], t["2d"], t["n"], t["e"], yfmt=lambda v: f"{v:g}")
    ns = np.linspace(4, 40, 37)
    ec = np.clip(np.cos(2 * math.pi / ns), 0, 1)
    cv.poly(d.pts(list(ns) + [40, 4], list(ec) + [0, 0]), C[0], 0.2)
    cv.line(d.pts(ns, ec), C[0], 1.8)
    cv.text(d.X(22), d.Y(0.35), t["mono"], 12, "middle")
    cv.text(d.X(9), d.Y(0.9), t["alias"], 11)
    return cv


def fig3(t):
    cv = Canvas(860, 520)
    a = Axes(cv, 80, 50, 740, 390, (0, 180), (-1.15, 1.15))
    a.frame([0, 30, 60, 90, 120, 150, 180], [-1, -0.5, 0, 0.5, 1], t["3t"], t["3x"], t["3y"], dx=48)
    cv.poly([(80, 50), (820, 50), (820, a.Y(0)), (80, a.Y(0))], C[2], 0.08)
    cv.poly([(80, a.Y(0)), (820, a.Y(0)), (820, 440), (80, 440)], C[3], 0.08)
    cv.text(812, 68, t["att"], 12, "end", color=C[2])
    cv.text(812, 430, t["rep"], 12, "end", color=C[3])
    ps = np.linspace(0, 180, 37)
    cv.line(a.pts(ps, -np.cos(np.radians(ps))), C[0], 1.6)
    rng = np.random.default_rng(20260921)
    px, py = [], []
    while len(px) < 40:
        A = rng.standard_normal((2, 2))
        dA = np.linalg.det(A)
        if abs(dA) < 1e-2:
            continue
        if dA < 0:
            A = A[::-1].copy(); dA = -dA
        S = A / math.sqrt(dA)
        tau = np.trace(S)
        if not 2.2 < tau < 8:
            continue
        _, vp, _, vm = real_hyperbolic_eigenvectors(S)
        X0 = rng.standard_normal(2)
        al, be = np.linalg.solve(np.c_[vp, vm], X0)
        if al < 0:
            vp, al = -vp, -al                              # 配置が入っている区画の二辺にそろえる（α, β > 0）
        if be < 0:
            vm, be = -vm, -be
        B = (S - np.linalg.inv(S)) / (2 * math.sinh(math.acosh(tau / 2)))
        mu = 2 * float((B @ X0) @ (B @ X0)) - 2 * float(X0 @ X0)
        px.append(math.degrees(math.acos(max(-1.0, min(1.0, float(vp @ vm))))))
        py.append(mu / (8 * al * be))
    cv.dots(a.pts(px, py), "#000", 2.6)
    cv.dots([(a.X(90), a.Y(0))], C[3], 5)
    cv.text(a.X(90) + 10, a.Y(0) + 20, t["3n"], 11)
    a.legend([(t["3c"], C[0], None), (t["3s"], "#000", "1 5")], 96, 70, 200)
    return cv


def fig4(t):
    cv = Canvas(920, 480)
    cv.text(460, 28, t["4t"], 14, "middle")
    a = Axes(cv, 40, 60, 380, 380, (-1.5, 1.5), (-1.5, 1.5))
    a.frame([], [], "", "", "")
    cv.text(230, 458, t["4ab"], 12, "middle")
    cv.line(a.pts([-1.5, 1.5], [0, 0]), "#999", 1); cv.line(a.pts([0, 0], [-1.5, 1.5]), "#999", 1)
    cv.line(a.pts([-1.4, 1.4], [-1.4, 1.4]), "#999", 1, "4 3")
    cv.text(a.X(1.42), a.Y(0) - 6, "a", 12, "end"); cv.text(a.X(0) + 6, a.Y(1.4), "b", 12)
    cv.text(a.X(1.25), a.Y(1.3), "a=b", 11, "end", color="#777")
    X = np.array([1.0, 0.4])
    imgs = {"X": (X, "#000"), "4P": (np.array([X[1], X[0]]), C[0]), "4D": (np.array([-X[0], X[1]]), C[1]), "4J": (np.array([-X[1], X[0]]), C[2])}
    J = imgs["4J"][0]
    cv.poly(a.pts([0, X[0], X[0] + J[0], J[0]], [0, X[1], X[1] + J[1], J[1]]), C[2], 0.15)
    for key, (v, col) in imgs.items():
        cv.arrow(a.X(0), a.Y(0), a.X(v[0]), a.Y(v[1]), col, 2)
        off = {"X": (8, 4), "4P": (8, -6), "4D": (-4, 20), "4J": (-10, -10)}[key]
        cv.text(a.X(v[0]) + off[0], a.Y(v[1]) + off[1], "X" if key == "X" else t[key], 12, "start" if key in ("X", "4P", "4D") else "end", color=col)
    cv.text(a.X(0.78), a.Y(1.18), t["4area"], 12, "start", color=C[2])
    for i, (key, col) in enumerate((("4x", C[0]), ("4y", C[1]), ("4r", C[2]))):
        cv.text(470, 92 + 24 * i, t[key], 13, color=col)
    cv.text(470, 92 + 24 * 3 + 4, t["4z"], 13)
    for i, s in enumerate(t["4f"]):
        y = 196 + 52 * i
        cv.rect(470, y, 420, 34, "#f4f4f4", "#888")
        cv.text(680, y + 22, s, 12, "middle")
        if i < 4:
            cv.arrow(680, y + 34, 680, y + 51, "#555", 1.5, 7)
    return cv


def fig5(t):
    cv = Canvas(920, 420)
    cv.text(460, 28, t["5t"], 14, "middle")
    cv.text(210, 70, t["5a"], 13, "middle")
    cv.line([(110, 150), (310, 150)], "#000", 2)
    cv.dots([(110, 150), (310, 150)], C[0], 9)
    cv.text(210, 138, "μ12 : + / 0 / −", 13, "middle")
    for i, k in enumerate(("5a1", "5a2", "5a3")):
        cv.text(40, 232 + 24 * i + (14 if i == 2 else 0), t[k], 12, color="#000" if i < 2 else C[3])
    cv.line([(440, 50), (440, 390)], "#ccc", 1)
    cv.text(690, 70, t["5b"], 13, "middle")
    tri = [(580, 250), (800, 250), (690, 100)]
    cv.line(tri + [tri[0]], "#000", 2)
    cv.dots(tri, C[0], 9)
    for (x, y), s in zip(tri, ("σ1", "σ2", "σ3")):
        cv.text(x + (-26 if x < 690 else 14), y + 6, s, 13)
    cv.text(690, 272, "s12", 13, "middle"); cv.text(762, 172, "s23", 13); cv.text(618, 172, "s31", 13, "end")
    for i, k in enumerate(("5b1", "5b2", "5b3")):
        cv.text(690, 318 + 24 * i, t[k], 12, "middle")
    return cv


FIGS = [("fig01_coulomb_readout_two_readouts", fig1), ("fig02_coulomb_readout_clock_and_force", fig2),
        ("fig03_coulomb_readout_coupling_sign", fig3), ("fig04_coulomb_readout_three_readouts", fig4),
        ("fig05_coulomb_readout_two_values_limit", fig5)]


def main():
    fp = None
    if HAVE_MPL:
        for c in JP_FONTS:
            if os.path.exists(c):
                fp = fm.FontProperties(fname=c)
                break
    for lang in ("ja", "en"):
        for name, fn in FIGS:
            cv = fn(L[lang])
            cv.to_svg(HERE / f"{name}_{lang}_v1.svg")
            if HAVE_MPL and (lang == "en" or fp is not None):
                cv.to_png(HERE / f"{name}_{lang}_v1.png", fp if lang == "ja" else None)
            print("written:", f"{name}_{lang}_v1", "(svg" + (", png)" if HAVE_MPL and (lang == "en" or fp is not None) else ")"))


if __name__ == "__main__":
    main()
