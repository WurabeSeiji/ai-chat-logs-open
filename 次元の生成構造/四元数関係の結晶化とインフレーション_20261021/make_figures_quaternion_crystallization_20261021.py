# -*- coding: utf-8 -*-
"""
四元数関係の結晶化とインフレーション（仮説論文）図版生成プログラム  2026-10-21
- すべての図を SVG で直接生成する（文字は <text> のまま残し、編集可能にする）
- PNG は cairosvg で SVG から変換する（任意）
- 乱数を使う図（fig05, fig06）は固定シードで再現可能
実行:  python make_figures_quaternion_crystallization_20261021.py [出力フォルダ]
依存:  numpy, scipy（ConvexHull, Delaunay）, cairosvg（PNG 変換のみ）
記法の宣言（検討ルール R2・R3・R4）:
  相対位相 u = q1^-1 q2（分母 = 観測者 q1、左除算）、結合係数 k = (q1+q2)^-1 (q1-q2)
  （分母 = 合成 q1+q2、読み出しの最後に一度だけ割る）。和と差の原点は q1, q2 共通の原点 O。
"""
import math, os, re, sys
import numpy as np
from scipy.spatial import ConvexHull, Delaunay

OUT = sys.argv[1] if len(sys.argv) > 1 else "out"
os.makedirs(OUT, exist_ok=True)
FONT = "'Noto Sans CJK JP','Yu Gothic','Meiryo','Hiragino Sans',sans-serif"
C = dict(ink="#222", gray="#888", light="#ccc", blue="#1f5fa8", red="#b03030",
         green="#2e7d32", fillb="#dce8f6", fillr="#f6dcdc", fillg="#e3f1e4")


def f(v):
    return f"{v:.1f}".rstrip("0").rstrip(".")


class SVG:
    def __init__(s, w, h, title):
        s.w, s.h, s.el, s.marks = w, h, [], set()
        s.el.append(f'<rect width="{w}" height="{h}" fill="#fff"/>')
        s.text(w / 2, 28, title, 18, "middle", weight="bold")

    def esc(s, t):
        return t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

    def mark(s, col):
        s.marks.add(col)
        return f'url(#a{col[1:]})'

    def line(s, x1, y1, x2, y2, col=C["ink"], sw=1.5, dash=None, arrow=False):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        a = f' marker-end="{s.mark(col)}"' if arrow else ""
        s.el.append(f'<line x1="{f(x1)}" y1="{f(y1)}" x2="{f(x2)}" y2="{f(y2)}" '
                    f'stroke="{col}" stroke-width="{sw}"{d}{a}/>')

    def segs(s, pairs, col=C["ink"], sw=1.0, op=1.0, dash=None):
        d = "".join(f"M{f(a[0])} {f(a[1])}L{f(b[0])} {f(b[1])}" for a, b in pairs)
        da = f' stroke-dasharray="{dash}"' if dash else ""
        s.el.append(f'<path d="{d}" stroke="{col}" stroke-width="{sw}" '
                    f'stroke-opacity="{op}" fill="none"{da}/>')

    def poly(s, pts, col=C["ink"], fill="none", sw=1.5, closed=True, dash=None):
        p = " ".join(f"{f(x)},{f(y)}" for x, y in pts)
        tag = "polygon" if closed else "polyline"
        da = f' stroke-dasharray="{dash}"' if dash else ""
        s.el.append(f'<{tag} points="{p}" stroke="{col}" fill="{fill}" stroke-width="{sw}"{da}/>')

    def circle(s, cx, cy, r, col=C["ink"], fill="none", sw=1.5, dash=None):
        da = f' stroke-dasharray="{dash}"' if dash else ""
        s.el.append(f'<circle cx="{f(cx)}" cy="{f(cy)}" r="{f(r)}" stroke="{col}" '
                    f'fill="{fill}" stroke-width="{sw}"{da}/>')

    def dot(s, x, y, col=C["ink"], r=4):
        s.el.append(f'<circle cx="{f(x)}" cy="{f(y)}" r="{r}" fill="{col}"/>')

    def ellipse(s, cx, cy, rx, ry, col=C["ink"], sw=1.2, dash=None, fill="none"):
        da = f' stroke-dasharray="{dash}"' if dash else ""
        s.el.append(f'<ellipse cx="{f(cx)}" cy="{f(cy)}" rx="{f(rx)}" ry="{f(ry)}" '
                    f'stroke="{col}" fill="{fill}" stroke-width="{sw}"{da}/>')

    def path(s, d, col=C["ink"], fill="none", sw=1.5, dash=None):
        da = f' stroke-dasharray="{dash}"' if dash else ""
        s.el.append(f'<path d="{d}" stroke="{col}" fill="{fill}" stroke-width="{sw}"{da}/>')

    def text(s, x, y, t, size=14, anchor="start", col=C["ink"], weight="normal"):
        # 添字は _{..}、上付きは ^{..} と書く（Unicode の添字字形に依存しない）
        t = s.esc(t)
        sm, d, cur, out = f(size * 0.72), size * 0.32, 0.0, []
        for tok in re.split(r"([_^]\{[^}]*\})", t):
            if not tok:
                continue
            if tok[0] in "_^" and tok[1:2] == "{":
                tgt, body, fs = (d if tok[0] == "_" else -1.3 * d), tok[2:-1], f' font-size="{sm}"'
            else:
                tgt, body, fs = 0.0, tok, ""
            dy = f' dy="{f(tgt - cur)}"' if abs(tgt - cur) > 1e-9 else ""
            out.append(f"<tspan{dy}{fs}>{body}</tspan>" if (dy or fs) else body)
            cur = tgt
        t = "".join(out)
        s.el.append(f'<text x="{f(x)}" y="{f(y)}" font-size="{size}" text-anchor="{anchor}" '
                    f'fill="{col}" font-weight="{weight}">{t}</text>')

    def lines(s, x, y, ts, size=14, anchor="start", col=C["ink"], dy=None):
        for i, t in enumerate(ts):
            s.text(x, y + i * (dy or size * 1.45), t, size, anchor, col)

    def box(s, x, y, w, h, col=C["ink"], fill="#fff", sw=1.5):
        s.el.append(f'<rect x="{f(x)}" y="{f(y)}" width="{f(w)}" height="{f(h)}" rx="8" '
                    f'stroke="{col}" fill="{fill}" stroke-width="{sw}"/>')

    def arc(s, cx, cy, r, a0, a1, col=C["ink"], sw=1.5):
        # 数学の角（反時計回り、ラジアン）。SVG は y 下向き
        x0, y0 = cx + r * math.cos(a0), cy - r * math.sin(a0)
        x1, y1 = cx + r * math.cos(a1), cy - r * math.sin(a1)
        large = 1 if (a1 - a0) % (2 * math.pi) > math.pi else 0
        s.path(f"M{f(x0)} {f(y0)}A{f(r)} {f(r)} 0 {large} 0 {f(x1)} {f(y1)}", col, sw=sw)

    def save(s, name):
        defs = "".join(
            f'<marker id="a{c[1:]}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" '
            f'markerHeight="7" orient="auto-start-reverse"><path d="M0 0L10 5L0 10z" fill="{c}"/></marker>'
            for c in sorted(s.marks))
        body = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{s.w}" height="{s.h}" '
                f'viewBox="0 0 {s.w} {s.h}" font-family="{FONT}">\n<defs>{defs}</defs>\n'
                + "\n".join(s.el) + "\n</svg>\n")
        p = os.path.join(OUT, name + ".svg")
        with open(p, "w", encoding="utf-8") as fh:
            fh.write(body)
        try:
            import cairosvg
            cairosvg.svg2png(bytestring=body.encode("utf-8"),
                             write_to=os.path.join(OUT, name + ".png"), output_width=s.w * 2)
        except ImportError:
            pass
        print("wrote", p, len(body), "bytes")


# ---------------------------------------------------------------- fig01
def fig01():
    s = SVG(920, 820, "図1　研究の経緯：行き詰まりが次の段階を強制する連鎖")
    steps = [
        ("① 中心投影（グノモン射影）と Σx_{n}² = R²", "曲率半径 R が光円錐と接続される可能性"),
        ("② 無名複素ゼロ閉塞 Σz_{n}² = 0", "複素ヒルベルト空間（清水『新版 量子論の基礎』）"),
        ("③ x²+y²+z²+t² = R²+Q_{1}²+Q_{2}²+Q_{3}² と無名性", "ラベルは読み出し後に付けた名前"),
        ("④ 完全ネットワーク（複素シンプレックス）", "関係の数 M = N(N−1)/2"),
        ("⑤ M² の自己相互作用", "高対称な集合で幾何級数的（インフレーション的）発展"),
        ("⑥ 四元数（単なる 4 次元ベクトル）", "無名性を破らず、非可換な積という構造をもつ"),
        ("⑦ 過冷却 ＋ シード → 結晶化（部分投射）", "本稿の主題：外から関係を与えない自律的な構造化"),
    ]
    why = [
        ["Σx_{n}² = 0 は実数では自明解のみ", "→ 複素数 z_{n} が必要"],
        ["四次元時空が背景として外から置かれ、", "進行位相 ωt が暗黙に持ち込まれている"],
        ["{z_{n}} そのものには構造がなく、", "四次元時空を自律的に読み出せない"],
        ["関係（どの頂点とどの頂点か）を", "外から与える → 無名性に違反"],
        ["発展は得られたが、", "四次元時空が読み出せない"],
        ["構造を自律的に生む仕組みが必要", "→ 物性物理の結晶化に類似"],
    ]
    x, w, h, gap, y0 = 40, 500, 74, 34, 56
    for i, (t1, t2) in enumerate(steps):
        y = y0 + i * (h + gap)
        last = i == len(steps) - 1
        s.box(x, y, w, h, C["blue"] if last else C["ink"], C["fillb"] if last else "#fff", 2 if last else 1.5)
        s.text(x + 16, y + 30, t1, 16, weight="bold")
        s.text(x + 16, y + 56, t2, 13, col="#444")
        if i < len(why):
            s.line(x + w / 2, y + h, x + w / 2, y + h + gap - 2, C["red"], 2, arrow=True)
            s.lines(x + w + 24, y + h + 4, why[i], 13, col=C["red"])
    s.text(x + w + 24, y0 + 6 * (h + gap) + 30, "（赤字＝その段の行き詰まり）", 12, col=C["gray"])
    s.save("fig01_研究経緯_必然の連鎖")


# ---------------------------------------------------------------- fig02
def fig02():
    s = SVG(920, 540, "図2　二等辺三角形と視野角距離：結合係数 g と距離 D")
    ox, oy, r = 190, 400, 150
    t1, t2 = math.radians(15), math.radians(95)
    z1 = (ox + r * math.cos(t1), oy - r * math.sin(t1))
    z2 = (ox + r * math.cos(t2), oy - r * math.sin(t2))
    zs = (z1[0] + z2[0] - ox, z1[1] + z2[1] - oy)
    s.arc(ox, oy, r, -0.25, math.pi + 0.25, C["light"], 1)
    s.line(ox, oy, *z1, C["ink"], 1.8)
    s.line(ox, oy, *z2, C["ink"], 1.8)
    s.poly([z1, zs, z2], C["gray"], sw=1, closed=False, dash="3 3")
    s.line(*z1, *z2, C["red"], 2.5)
    s.line(ox, oy, *zs, C["blue"], 2.5)
    s.arc(ox, oy, 38, t1, t2, C["ink"], 1.3)
    s.text(ox + 34, oy - 40, "Δθ", 15)
    for p, lab, dx, dy in [((ox, oy), "O（観測者・原点）", -10, 22), (z1, "z_{1}", 8, 6),
                           (z2, "z_{2}", -26, -6), (zs, "z_{1}+z_{2}", 8, -4)]:
        s.dot(*p)
        s.text(p[0] + dx, p[1] + dy, lab, 14)
    s.text(20, 72, "和（青）：|z_{1}+z_{2}| = 2r cos(Δθ/2)", 14, col=C["blue"])
    s.text(20, 94, "差（赤）：|z_{1}−z_{2}| = 2r sin(Δθ/2)", 14, col=C["red"])
    s.text(20, 520, "等振幅 |z_{1}| = |z_{2}| = r。和と差の原点は O（R3）", 12, col=C["gray"])
    # 右：グラフ
    gx, gy, gw, gh = 500, 410, 340, 300
    s.line(gx, gy, gx + gw, gy, arrow=True, col=C["ink"])
    s.line(gx, gy, gx, gy - gh, arrow=True, col=C["ink"])
    ymax = 4.0
    X = lambda a: gx + gw * 0.95 * a / math.pi
    Y = lambda v: gy - gh * 0.9 * min(v, ymax) / ymax
    for v in (1, 2, 3, 4):
        s.text(gx - 8, Y(v) + 5, str(v), 12, "end", C["gray"])
        s.line(gx, Y(v), gx + gw * 0.95, Y(v), C["light"], 0.6)
    for a, lab in ((0, "0"), (math.pi / 2, "π/2"), (math.pi, "π")):
        s.text(X(a), gy + 20, lab, 13, "middle")
    s.text(gx + gw + 8, gy + 5, "Δθ", 14)
    A = np.linspace(0.02, math.pi - 0.02, 90)
    s.poly([(X(a), Y(math.tan(a / 2))) for a in A if math.tan(a / 2) <= ymax],
           C["blue"], sw=2.4, closed=False)
    s.poly([(X(a), Y(1 / math.tan(a / 2))) for a in A if 1 / math.tan(a / 2) <= ymax],
           C["red"], sw=2.4, closed=False)
    s.dot(X(math.pi / 2), Y(1), C["ink"], 4)
    s.text(X(math.pi / 2) + 8, Y(1) + 18, "直交：g = D = 1", 12)
    s.text(gx + 20, 64, "g = |z_{1}−z_{2}| / |z_{1}+z_{2}| = tan(Δθ/2)　結合係数＝近さ", 13, col=C["blue"])
    s.text(gx + 20, 86, "D = 1/g = cot(Δθ/2)　視野角距離", 13, col=C["red"])
    s.lines(gx - 20, 470, ["同位相 Δθ→0：g→0, D→∞（作用なし・無限遠）",
                           "逆位相 Δθ=π：g→∞, D=0（全反射・接触）",
                           "分母は合成 |z_{1}+z_{2}|。読み出しの最後に一度だけ割る（R2）"], 12, dy=19)
    s.save("fig02_視野角距離と結合係数")


# ---------------------------------------------------------------- fig03
def fig03():
    s = SVG(920, 560, "図3　相対位相の円と二枚の地図：結合の地図 k と距離の地図 x")
    cx, cy, r = 280, 330, 140
    a = math.radians(64)
    U = lambda ang: (cx + r * math.cos(ang), cy - r * math.sin(ang))
    s.line(cx - r - 60, cy, cx + r + 60, cy, C["light"], 1)
    s.line(cx, cy + r + 40, cx, 60, C["ink"], 1.2, arrow=True)
    s.text(cx + 8, 70, "射影先の直線（Re u = 0）", 12, col=C["gray"])
    s.circle(cx, cy, r, C["ink"], sw=1.6)
    u = U(a)
    m, p = U(math.pi), U(0)
    tk, tx = math.tan(a / 2), 1 / math.tan(a / 2)
    K, Xp = (cx, cy - r * tk), (cx, cy - r * tx)
    s.line(*m, *K, C["blue"], 2)
    s.line(*K, *u, C["blue"], 2, dash="5 3")
    s.line(*p, *Xp, C["red"], 2)
    s.dot(*u, C["ink"], 5); s.text(u[0] + 8, u[1] - 6, "u = q_{1}^{−1}q_{2} = e^{nα}", 14)
    s.dot(*p, C["green"], 6); s.text(p[0] + 10, p[1] + 20, "u = +1：自己（観測者）", 13, col=C["green"])
    s.dot(*m, C["ink"], 5); s.text(m[0] - 10, m[1] + 22, "u = −1：逆位相", 13, anchor="end")
    s.dot(*K, C["blue"], 5); s.text(K[0] - 10, K[1] + 4, "|k| = tan(α/2)", 14, "end", C["blue"])
    s.dot(*Xp, C["red"], 5); s.text(Xp[0] + 10, Xp[1] + 4, "|x| = cot(α/2)", 14, col=C["red"])
    s.arc(cx, cy, 30, 0, a, C["ink"], 1.2); s.text(cx + 34, cy - 14, "α", 15)
    tx0 = 520
    s.lines(tx0, 110, ["結合の地図（−1 から射影）", "k = (1+u)^{−1}(1−u) = −n tan(α/2)",
                       "観測者 u = +1 は原点 k = 0", "同位相ほど k は小さい（作用が弱い）"], 14, col=C["blue"])
    s.lines(tx0, 250, ["距離の地図（+1 から射影）", "x = (1−u)^{−1}(1+u) = 1/k",
                       "観測者 u = +1 は無限遠（自己距離 = ∞）", "逆位相 u = −1 は x = 0（接触）"], 14, col=C["red"])
    s.lines(tx0, 390, ["二枚は反転 x = 1/k で結ばれ、どちらも角度を保つ（共形）",
                       "→ 視野角だけで測る局在観測者には区別できない",
                       "分母：相対位相は観測者 q_{1} で左除算、k は合成 q_{1}+q_{2} で除算",
                       "（R2・R4。複素数では n = i）"], 13)
    s.save("fig03_相対位相の円と二枚の地図")


# ---------------------------------------------------------------- fig04
def fig04():
    s = SVG(920, 520, "図4　四元数への拡張：S³ の相対位相を R³ へ射影する")
    cx, cy, R, e = 230, 290, 150, math.radians(16)
    P = lambda X, Y, Z: (cx + R * X, cy - R * (Z * math.cos(e) + Y * math.sin(e)))
    s.circle(cx, cy, R, C["ink"], sw=1.5)
    s.ellipse(cx, cy, R, R * math.sin(e), C["gray"], 1, dash="4 3", fill="#f3f6fa")
    N, S0 = P(0, 0, 1), P(0, 0, -1)
    al, ph = math.radians(58), math.radians(35)
    X, Y, Z = math.sin(al) * math.cos(ph), math.sin(al) * math.sin(ph), math.cos(al)
    u = P(X, Y, Z)
    k = P(X / (1 + Z), Y / (1 + Z), 0)
    O = P(0, 0, 0)
    s.line(*S0, *u, C["blue"], 1.6, dash="5 3")
    s.line(*O, *k, C["blue"], 3, arrow=True)
    for pt, lab, col, dx, dy in [(N, "u = +1（観測者）", C["green"], 10, -6),
                                 (S0, "u = −1（射影の極）", C["ink"], 10, 18),
                                 (u, "u = q_{1}^{−1}q_{2}", C["ink"], 8, -8),
                                 (O, "k = 0", C["green"], -50, 16)]:
        s.dot(*pt, col, 5)
        s.text(pt[0] + dx, pt[1] + dy, lab, 13, col=col)
    s.text(k[0] + 8, k[1] + 18, "k = −n tan(α/2)", 14, col=C["blue"])
    s.text(cx - R, cy + R + 36, "S³ を 1 次元省略して S² で表示。赤道面 = R³（純虚四元数）", 12, col=C["gray"])
    # 右：左右の枠
    bx, by = 600, 200
    s.circle(bx, by, 100, C["light"], dash="3 3", sw=1)
    aL, aR = math.radians(20), math.radians(78)
    eL = (bx + 100 * math.cos(aL), by - 100 * math.sin(aL))
    eR = (bx + 100 * math.cos(aR), by - 100 * math.sin(aR))
    s.line(bx, by, *eL, C["blue"], 3, arrow=True)
    s.line(bx, by, *eR, C["red"], 3, arrow=True)
    s.dot(bx, by)
    s.text(eL[0] + 6, eL[1] + 4, "n（左：q_{1}^{−1}q_{2}）", 13, col=C["blue"])
    s.text(eR[0] + 6, eR[1] - 4, "q_{1} n q_{1}^{−1}（右：q_{2}q_{1}^{−1}）", 13, col=C["red"])
    s.lines(500, 335, ["長さ |k| = tan(α/2) は共通（観測枠によらない）",
                       "方向 n は観測枠に依存（q_{1} による共役で回転）"], 13)
    s.box(480, 390, 420, 90, C["ink"], "#fafafa", 1)
    s.lines(496, 416, ["振幅が異なるとき（u = ρ e^{nα}）",
                       "k = [(1−ρ²) − 2 Im u] / |1+u|²",
                       "スカラー部：振幅の不整合　ベクトル部：向きをもつ位相結合"], 13)
    s.save("fig04_四元数への拡張_S3からR3")


# ---------------------------------------------------------------- fig05
def fig05():
    rng = np.random.default_rng(20261021)
    n = 30
    V = rng.normal(size=(n, 3)); V /= np.linalg.norm(V, axis=1)[:, None]
    hull = ConvexHull(V)
    s = SVG(920, 500, "図5　無名性：名前のない点集合の凸包が接続（トラス）を与える")
    cx, cy, R, e = 220, 270, 170, math.radians(20)
    rot = lambda p: (cx + R * p[0], cy - R * (p[2] * math.cos(e) + p[1] * math.sin(e)),
                     p[1] * math.cos(e) - p[2] * math.sin(e))   # 第 3 成分 < 0 が手前
    edges = set()
    for sim in hull.simplices:
        for i in range(3):
            a, b = sorted((sim[i], sim[(i + 1) % 3])); edges.add((a, b))
    front, back = [], []
    for a, b in edges:
        pa, pb = rot(V[a]), rot(V[b])
        (front if (pa[2] + pb[2]) < 0 else back).append(((pa[0], pa[1]), (pb[0], pb[1])))
    s.circle(cx, cy, R, C["light"], sw=1)
    s.segs(back, C["gray"], 0.8, 0.35)
    s.segs(front, C["blue"], 1.3)
    for p in V:
        q = rot(p); s.dot(q[0], q[1], C["ink"] if q[2] < 0 else C["gray"], 3)
    s.text(cx, cy + R + 34, f"点 {n} 個（ラベルなし）→ 凸包の面 {len(hull.simplices)} 枚・辺 {len(edges)} 本",
           13, "middle")
    s.text(cx, cy + R + 52, "S³ の類比として S² で表示（3 次元では面は四面体）", 12, "middle", C["gray"])
    # 右：観測者中心の立体射影（k の地図）
    o = V[0]; pole = -o
    t = np.cross(o, [0, 0, 1.0]); e1 = t / np.linalg.norm(t); e2 = np.cross(o, e1)
    proj = lambda x: np.array([x @ e1, x @ e2]) / (1 + x @ o)
    Q = np.array([proj(x) for x in V])
    good = []
    for sim, eq in zip(hull.simplices, hull.equations):
        nrm, d = eq[:3], -eq[3]
        if nrm @ pole < d - 1e-12:        # 極が面のキャップの外 → ドロネー三角形
            good.append(tuple(sorted(sim)))
    dl = {tuple(sorted(x)) for x in Delaunay(Q).simplices}
    print("fig05: hull->Delaunay match:", set(good) == dl, len(good), len(dl))
    px, py, sc, W = 680, 270, 95, 210
    inside = lambda q: abs(q[0]) * sc < W and abs(q[1]) * sc < W - 20
    pe = set()
    for sim in good:
        for i in range(3):
            pe.add(tuple(sorted((sim[i], sim[(i + 1) % 3]))))
    segs = [((px + sc * Q[a][0], py - sc * Q[a][1]), (px + sc * Q[b][0], py - sc * Q[b][1]))
            for a, b in pe if inside(Q[a]) and inside(Q[b])]
    s.segs(segs, C["blue"], 1.2)
    # 空の外接円の例：窓に収まる外接円のうち最大のもの
    def circ(sm):
        A, B, Cc = Q[list(sm)]
        d = 2 * (A[0] * (B[1] - Cc[1]) + B[0] * (Cc[1] - A[1]) + Cc[0] * (A[1] - B[1]))
        ux = ((A @ A) * (B[1] - Cc[1]) + (B @ B) * (Cc[1] - A[1]) + (Cc @ Cc) * (A[1] - B[1])) / d
        uy = ((A @ A) * (Cc[0] - B[0]) + (B @ B) * (A[0] - Cc[0]) + (Cc @ Cc) * (B[0] - A[0])) / d
        return A, B, Cc, ux, uy, float(np.linalg.norm(A - [ux, uy]))
    cand = [circ(sm) for sm in good if 0 not in sm]
    cand = [c for c in cand if (abs(c[3]) + c[5]) * sc < W - 30 and (abs(c[4]) + c[5]) * sc < W - 50]
    A, B, Cc, ux, uy, rr = max(cand, key=lambda c: c[5])
    s.poly([(px + sc * p[0], py - sc * p[1]) for p in (A, B, Cc)], C["red"], C["fillr"], 1.8)
    s.circle(px + sc * ux, py - sc * uy, sc * rr, C["red"], sw=1.4, dash="5 3")
    for q in Q:
        if inside(q): s.dot(px + sc * q[0], py - sc * q[1], C["ink"], 3)
    s.dot(px, py, C["green"], 6); s.text(px - 10, py + 22, "観測者（k = 0）", 12, "end", C["green"])
    s.text(px, cy + R + 34, "観測者中心の立体射影 → ドロネー分割（一致を数値確認）", 13, "middle")
    s.text(px, cy + R + 52, "赤：外接円の内部に他の点がない（共形変換で不変な判定）", 12, "middle", C["red"])
    s.save("fig05_無名点集合の凸包とドロネー分割")


# ---------------------------------------------------------------- fig06
def fig06():
    rng = np.random.default_rng(7)
    n = 13
    P = rng.uniform(-1, 1, size=(n, 2))
    s = SVG(920, 440, "図6　投射の三相（2 次元の類比：剛性の数え方 F − S = 2V − 3 − E）")
    order = np.argsort(np.linalg.norm(P - P.mean(0), axis=1))
    placed, E_part = [order[0], order[1]], [(order[0], order[1])]
    for v in order[2:]:            # ヘネバーグ操作：新しい点を既存の最近接 2 点に結ぶ
        dd = sorted(placed, key=lambda w: np.linalg.norm(P[v] - P[w]))[:2]
        E_part += [(v, dd[0]), (v, dd[1])]; placed.append(v)
    E_full = [(i, j) for i in range(n) for j in range(i + 1, n)]
    panels = [("ゼロ投射", [], "構造ができない"),
              ("部分投射", E_part, "剛性ぎりぎり・組み替えの余地"),
              ("全投射", E_full, "剛体結晶・発展しない")]
    for k, (name, E, note) in enumerate(panels):
        ox = 160 + k * 300; oy = 210; sc = 105
        X = lambda p: (ox + sc * p[0], oy - sc * p[1])
        s.box(ox - 140, 54, 280, 310, C["blue"] if k == 1 else C["light"],
              C["fillb"] if k == 1 else "#fff", 2 if k == 1 else 1)
        s.text(ox, 80, name, 16, "middle", weight="bold")
        s.segs([(X(P[a]), X(P[b])) for a, b in E], C["ink"], 0.8 if k == 2 else 1.4, 0.6 if k == 2 else 1)
        for i, p in enumerate(P):
            s.dot(*X(p), C["red"] if (k == 1 and i == order[0]) else C["ink"], 4)
        Ecount = len(E); val = 2 * n - 3 - Ecount
        Fm = max(val, 0); Sm = max(-val, 0)
        s.text(ox, 340, f"E = {Ecount}　F = {Fm}　S = {Sm}", 13, "middle")
        s.text(ox, 390, note, 13, "middle", C["blue"] if k == 1 else C["ink"])
    s.text(460, 418, f"V = {n}、2V − 3 = {2 * n - 3}。部分投射はシード（赤）からヘネバーグ操作で成長させた例。"
           "3 次元では 3V − 6", 12, "middle", C["gray"])
    s.save("fig06_投射の三相と剛性の数え方")


# ---------------------------------------------------------------- fig07
def fig07():
    s = SVG(920, 460, "図7　四面体のフラストレーション：平坦な空間では閉じず、S³ では閉じる")
    dih = math.degrees(math.acos(1 / 3))
    for k, (cx, ang, title) in enumerate([(230, dih, "平坦な R³：二面角 70.53° × 5 = 352.64°"),
                                          (680, 72.0, "S³（600 胞体）：二面角 72° × 5 = 360°")]):
        cy, r = 250, 150
        a0 = 90
        for i in range(5):
            t0, t1 = math.radians(a0 + i * ang), math.radians(a0 + (i + 1) * ang)
            pts = [(cx, cy), (cx + r * math.cos(t0), cy - r * math.sin(t0)),
                   (cx + r * math.cos(t1), cy - r * math.sin(t1))]
            s.poly(pts, C["blue"], C["fillb"] if i % 2 == 0 else "#eef4fb", 1.4)
        if k == 0:
            t0, t1 = math.radians(a0 + 5 * ang), math.radians(a0 + 360)
            s.poly([(cx, cy), (cx + r * math.cos(t0), cy - r * math.sin(t0)),
                    (cx + r * math.cos(t1), cy - r * math.sin(t1))], C["red"], C["fillr"], 1.6)
            s.text(cx + 22, cy - r - 10, f"隙間 {360 - 5 * dih:.2f}°", 14, col=C["red"])
        s.dot(cx, cy)
        s.text(cx, 64, title, 15, "middle", weight="bold")
    s.lines(230, 420, ["共通の辺を真上から見た図（各扇形が正四面体 1 個）",
                       "隙間＝閉じきれない分＝欠陥線（不足角）＝曲率の担い手"], 12, "middle", C["gray"], dy=18)
    s.lines(680, 420, ["頂点 120 個 = 単位四元数（二項正二十面体群）",
                       "正四面体 600 個で S³ をちょうど埋め尽くす"], 12, "middle", C["gray"], dy=18)
    s.save("fig07_四面体フラストレーションと600胞体")


# ---------------------------------------------------------------- fig08
def fig08():
    s = SVG(920, 510, "図8　シードからの幾何級数成長：面の数がプロセスごとに r 倍になる")
    # 左：面の木
    x0, y0 = 40, 70
    s.text(x0, y0, "四面体を 1 個付けると、面 1 枚を塞ぎ新しい面 3 枚を露出する", 13)
    lv = [1, 3, 9, 27]
    for L, cnt in enumerate(lv):
        y = y0 + 50 + L * 80
        xs = [60 + (330 / (cnt + 1)) * (i + 1) for i in range(cnt)]
        for xx in xs: s.dot(xx, y, C["blue"], 4 if cnt < 27 else 3)
        if L > 0:
            prev = [60 + (330 / (lv[L - 1] + 1)) * (i + 1) for i in range(lv[L - 1])]
            s.segs([((prev[i // 3], y - 80), (xs[i], y)) for i in range(cnt)], C["gray"], 0.9)
        s.text(410, y + 5, f"n = {L}：F = {cnt}", 13)
    s.lines(x0, 420, ["隣どうしが面を共有して閉じると倍率 r < 3", "閉じが支配すると r → 1（インフレーションの終わり）"],
            13, col=C["red"])
    # 右：グラフ
    gx, gy, gw, gh = 560, 400, 320, 300
    s.line(gx, gy, gx + gw, gy, arrow=True); s.line(gx, gy, gx, gy - gh, arrow=True)
    nmax, ymax = 200, 100
    X = lambda n: gx + gw * 0.95 * n / nmax
    Y = lambda v: gy - gh * 0.92 * v / ymax
    for v in (0, 25, 50, 75, 100):
        s.text(gx - 6, Y(v) + 4, str(v), 11, "end", C["gray"])
        s.line(gx, Y(v), gx + gw * 0.95, Y(v), C["light"], 0.6)
    for n in (0, 50, 100, 150, 200): s.text(X(n), gy + 18, str(n), 11, "middle")
    s.text(gx + gw, gy + 34, "プロセス回数 n", 12, "end")
    s.text(gx - 40, gy - gh - 8, "log_{10}（ブロック数 V）", 12)
    for r, col in ((3, C["blue"]), (2.5, C["green"]), (2, C["gray"])):
        pts = [(X(n), Y(n * math.log10(r))) for n in range(0, nmax + 1, 5) if n * math.log10(r) <= ymax]
        s.poly(pts, col, sw=2, closed=False)
        s.text(pts[-1][0] - 4, pts[-1][1] - 6, f"r = {r}", 12, "end", col)
    # 閉じのモデル（例示）：r_n = 1 + 2(1 − V_n / V_max)
    lv_, lmax, pts = 0.0, 78.0, []
    for n in range(0, nmax + 1):
        if n % 4 == 0 or 150 <= n <= 175:
            pts.append((X(n), Y(lv_)))
        r_n = 1 + 2 * max(0.0, 1 - 10 ** (lv_ - lmax))
        lv_ += math.log10(r_n)
    s.poly(pts, C["red"], sw=2.2, closed=False, dash="6 3")
    s.text(X(8), Y(96), "破線：閉じを含む例示モデル", 12, col=C["red"])
    s.text(X(8), Y(96) + 17, "r_{n} = 1 + 2(1 − V_{n}/V_{max})、V_{max} = 10^{78}", 12, col=C["red"])
    n60 = 180 / math.log(3)
    s.line(X(n60), gy, X(n60), Y(n60 * math.log10(3)), C["ink"], 1, dash="3 3")
    s.text(X(n60) + 4, gy - 8, f"n ≈ {n60:.0f}", 11)
    s.lines(gx - 40, 462, ["e-fold 数 N_{e} = (n/3) ln r（長さの尺度 a ∝ V^{1/3}）",
                           f"r = 3 で 60 e-fold ⇔ n ≈ {n60:.0f}（V ≈ 10^{{{n60 * math.log10(3):.0f}}}）"], 12, dy=19)
    s.save("fig08_シードからの幾何級数成長とe-fold")


# ---------------------------------------------------------------- fig09
def fig09():
    s = SVG(920, 500, "図9　時空の式は分解型八元数のヌル条件 ＝ 等振幅の四元数の組")
    gx, gy, L = 90, 450, 320
    s.line(gx, gy, gx + L + 30, gy, arrow=True); s.line(gx, gy, gx, gy - L - 30, arrow=True)
    s.text(gx + L + 30, gy + 24, "|a|　a = (x, y, z, t)：時空側", 13, "end")
    s.text(gx - 10, gy - L - 40, "|b|　b = (R, Q_{1}, Q_{2}, Q_{3})：内部側", 13)
    s.poly([(gx, gy), (gx + L, gy), (gx + L, gy - L)], "none", "#f4f8fc", 0)
    s.poly([(gx, gy), (gx, gy - L), (gx + L, gy - L)], "none", "#fcf4f4", 0)
    s.line(gx, gy, gx + L, gy - L, C["blue"], 3)
    s.text(gx + L + 8, gy - L + 4, "N = |a|² − |b|² = 0", 14, col=C["blue"])
    s.text(gx + L * 0.72, gy - L * 0.25, "N > 0", 14, "middle", C["gray"])
    s.text(gx + L * 0.2, gy - L * 0.9, "N < 0", 14, "middle", C["gray"])
    on = (gx + L * 0.55, gy - L * 0.55); off = (gx + L * 0.55, gy - L * 0.78)
    s.dot(*on, C["blue"], 6); s.text(on[0] + 10, on[1] + 18, "等振幅：u = a^{−1}b ∈ S³", 13, col=C["blue"])
    s.dot(*off, C["red"], 6); s.text(off[0] - 12, off[1] + 5, "|a| ≠ |b|：k にスカラー部", 13, "end", C["red"])
    s.line(*off, on[0], on[1] - 8, C["red"], 1.2, dash="3 3", arrow=True)
    tx = 530
    s.lines(tx, 110, ["ケーリー＝ディクソン構成：o = (a, b)（a, b は四元数）",
                     "分解型八元数のノルム　N(o) = |a|² − |b|²（符号 (4,4)）"], 14)
    s.lines(tx, 190, ["x²+y²+z²+t² − (R²+Q_{1}²+Q_{2}²+Q_{3}²) = 0", "　⇔　|a| = |b|",
                      "　⇔　相対位相 u = a^{−1}b が単位四元数（S³ 上）"], 14, col=C["blue"])
    s.lines(tx, 305, ["→ 本稿の「等振幅の四元数二つの関係読み出し」と同じもの",
                      "→ 振幅の不整合（ヌル錐からのずれ）が重力歪みの候補"], 13)
    s.lines(tx, 395, ["宣言（R2・R3・R4・R7）：分母は時空側 a（左除算）、",
                      "二乗の和の原点は a, b 共通の原点 O。",
                      "a^{−1}b と b a^{−1} のどちらで読むかは外から入る選択"], 12, col=C["gray"])
    s.save("fig09_分解型八元数のヌル条件と等振幅")


if __name__ == "__main__":
    for fn in (fig01, fig02, fig03, fig04, fig05, fig06, fig07, fig08, fig09):
        fn()
