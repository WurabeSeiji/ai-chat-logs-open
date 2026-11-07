#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""HTML の対照（2026-11-07、4 パルス版）。原本 index_v8q3.html（3 パルス）と index_v9q4.html（4 パルス、A₄ = 0・n₄ = k）を同じ 1〜3 の配置で比べる。
(1) h1：v9q4 からパルス 4 の項（4: θ = …、Θ の 4 つ目、軌道の 4: (…)、平行移動の 4:、クリップの 4:）を除くと v8q3 の h1 と完全一致（和の量・1〜3 の読み・軌道 1〜3・枠）。
(2) h2：対 12・23・13 の「対ごとの読み」、合成則〜鏡映対の条件、和の尾（c₀/R … 巻き数）が v8q3 と一致。
(3) v9q4 自身：座標の軌道 4 本が n_i ∘ (Θ₁, Θ₂, Θ₃)、Θ_i = θ_i(初期) + ティック × φ_i（i = 1..4）、Θ_i ≡ θ_i (mod 360)、右の枠 = max(360, max|点|)、点数 = ティック + 1。
    既定（newdefault）の初期値（1〜3 と、4: θ₄ = 300, φ₄ = 1, δ₄ = −10, M₄ = 100, A₄ = 1, n₄ = (0.577, 0.577, 0.577)）。h1・h2 に NaN・undefined がない（A₄ = 0 の対照では 4 を含む対の読みは R₄ = 0 で割るので対照外、4 の項を除いた h1 と 12・23・13・和の尾で調べる）。"""
import re, os, sys, html
base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
NUM = r"(-?\d+\.\d)"
ORB = r"　軌道（i_i Θ₁, j_i Θ₂, k_i Θ₃）：Θ = \(%s, %s, %s, %s\)、1: \(%s, %s, %s\)、2: \(%s, %s, %s\)、3: \(%s, %s, %s\)、4: \(%s, %s, %s\)、各 (\d+) 点、右の枠 ±(\d+)°" % ((NUM,) * 16)
TH0 = {"default": [180, 180, 180, 300], "tri": [180, 150, 100, 300], "new": [97, 194, 58, 300]}
def head(path):
    if not os.path.exists(path): return None
    s = open(path, encoding="utf-8").read()
    h1 = re.search(r'id="h1">(.*?)</div>', s, re.S); h2 = re.search(r'id="h2"[^>]*>(.*?)</div>', s, re.S)
    return None if not h1 or not h2 else (html.unescape(h1.group(1)), html.unescape(h2.group(1)))
def strip4(h1):
    """v9q4 の h1 からパルス 4 の項を除く"""
    h = re.sub(r"　4: θ = .*? R4 = \d+\.\d+", "", h1, count=1)
    h = re.sub(r"(Θ = \(-?\d+\.\d, -?\d+\.\d, -?\d+\.\d), -?\d+\.\d\)", r"\1)", h, count=1)
    h = re.sub(r"、4: \(-?\d+\.\d, -?\d+\.\d, -?\d+\.\d\)(、各 \d+ 点)", r"\1", h, count=1)
    h = re.sub(r"(平行移動 1: \([^)]*\)、2: \([^)]*\)、3: \([^)]*\))、4: \([^)]*\)", r"\1", h, count=1)
    h = re.sub(r"(1: -?[\d.]+°〜\+?-?[\d.]+°、2: -?[\d.]+°〜\+?-?[\d.]+°、3: -?[\d.]+°〜\+?-?[\d.]+°)、4: -?[\d.]+°〜\+?-?[\d.]+°", r"\1", h, count=1)
    return h
def seg_pairs(h2):
    return {lab: re.search(lab + r"：Δ = .*?vec V/\(R_iR_j\) = \([^)]*\)", h2).group(0) if re.search(lab + r"：Δ = ", h2) else None for lab in ["12", "23", "13"]}
def seg_rules(h2):
    m = re.search(r"合成則 .*?鏡映対の条件 δ₂＋δ₃ = -?\d+", h2); return m.group(0) if m else None
def seg_tail(h2):
    i = h2.find("｜　和：c₀/R"); return h2[i:] if i >= 0 else None
def check_orb(b, tag, th0):
    m = re.search(ORB, b[0])
    if not m: print(f"  {tag}軌道の項目がない"); return 1
    g = [float(m.group(i)) for i in range(1, 17)]; Th = g[0:4]; pts = [g[4:7], g[7:10], g[10:13], g[13:16]]; npts = int(m.group(17)); L2 = int(m.group(18))
    ths = [float(t) for t in re.findall(r"\d: θ = (-?\d+)°", b[0])]; phis = [float(t) for t in re.findall(r"\d: θ = -?\d+°, φ = (\d+)°", b[0])]
    nvs = [[float(x) for x in t] for t in re.findall(r"n = \((-?\d+\.\d+), (-?\d+\.\d+), (-?\d+\.\d+)\)", b[0])]
    t = re.search(r"ティック (\d+)", b[1]); tk = int(t.group(1)) if t else None
    ok1 = len(nvs) == 4 and all(abs(pts[i][a] - nvs[i][a] * Th[a]) < 0.3 for i in range(4) for a in range(3))   # n は 3 桁表示なので 0.577×300 で 0.15 程度の丸め
    ok2 = tk is not None and len(phis) == 4 and all(abs(Th[i] - (th0[i] + tk * phis[i])) < 0.06 for i in range(4))
    ok3 = len(ths) == 4 and all(min(abs((Th[i] - ths[i]) % 360), abs(360 - (Th[i] - ths[i]) % 360)) < 0.06 for i in range(4))
    ok4 = L2 == round(max(360.0, max(abs(v) for q in pts for v in q)))
    okn = tk is not None and npts == tk + 1
    print(f"  {tag}Θ = ({Th[0]:.1f}, {Th[1]:.1f}, {Th[2]:.1f}, {Th[3]:.1f})　点 " + "、".join(f"{i + 1}: ({pts[i][0]:.1f}, {pts[i][1]:.1f}, {pts[i][2]:.1f})" for i in range(4)) + f" = n_i ∘ (Θ₁, Θ₂, Θ₃)：{'一致' if ok1 else '不一致'}　Θ_i = θ_i(初期) + ティック × φ_i（4 本）：{'一致' if ok2 else '不一致'}　Θ_i ≡ θ_i (mod 360)：{'一致' if ok3 else '不一致'}　右の枠 ±{L2}°：{'一致' if ok4 else '不一致'}　点数 {npts}（ティック {tk if tk is not None else '?'} + 1）{'一致' if okn else '不一致'}")
    return (not ok1) + (not ok2) + (not ok3) + (not ok4) + (not okn)
def nonan(b, tag):
    bad = any(w in b[0] or w in b[1] for w in ["NaN", "undefined", "Infinity°"])
    print(f"  {tag}NaN・undefined：{'なし' if not bad else 'あり'}"); return int(bad)
def main():
    bad = 0
    for name in ["default", "tri"]:
        print(f"== {name}（A₄ = 0、n₄ = k、1〜3 は n 全部 k）==")
        for pw in (0, 1, 2):
            for rc in (0, 1):
                a = head(os.path.join(base, f"dom_v8q3_pow{pw}_rc{rc}_{name}.html")); b = head(os.path.join(base, f"dom_v9q4_pow{pw}_rc{rc}_{name}.html"))
                if a is None or b is None: print(f"  pow={pw} rc={rc} 読めない"); bad += 1; continue
                ok1 = strip4(b[0]) == a[0]
                pa, pb = seg_pairs(a[1]), seg_pairs(b[1]); ok2 = all(pa[l] is not None and pa[l] == pb[l] for l in pa)
                ok3 = seg_rules(a[1]) is not None and seg_rules(a[1]) == seg_rules(b[1]); ok4 = seg_tail(a[1]) is not None and seg_tail(a[1]) == seg_tail(b[1])
                print(f"  表示 {pw}、{'座標' if rc else '円周上'}：h1（4 の項を除く）{'一致' if ok1 else '不一致'}　h2 対 12・23・13 {'一致' if ok2 else '不一致'}　合成則〜鏡映対 {'一致' if ok3 else '不一致'}　和の尾 {'一致' if ok4 else '不一致'}")
                bad += (not ok1) + (not ok2) + (not ok3) + (not ok4) + nonan((strip4(b[0]), (seg_rules(b[1]) or "") + "".join(v or "" for v in seg_pairs(b[1]).values()) + (seg_tail(b[1]) or "")), "  ")   # A₄ = 0 では 4 を含む対の読み（R₄ = 0 で割る tan α・cos Θ）が NaN になるので対照外
                if rc == 1: bad += check_orb(b, "  ", TH0[name])
    print("== newdefault（rc=1）=="); b = head(os.path.join(base, "dom_v9q4_newdefault_rc1.html"))
    if b is None: print("  読めない"); bad += 1
    else:
        ok = all(w in b[0] for w in ["1: θ = 97°, φ = 1°, δ = 20°/t, M = 100", "2: θ = 194°, φ = 1°, δ = 0°/t, M = 100", "3: θ = 58°, φ = 1°, δ = 10°/t, M = 100", "4: θ = 300°, φ = 1°, δ = -10°/t, M = 100（奇 50・偶 50）A = 1.0 n = (0.577, 0.577, 0.577)", "n = (1.000, 0.000, 0.000) R1", "n = (0.000, 1.000, 0.000) R2", "n = (0.000, 0.000, 1.000) R3"])
        print(f"  初期値（1〜3 と 4: θ₄ = 300, φ₄ = 1, δ₄ = −10, M₄ = 100, A₄ = 1, n₄ = (1,1,1)/√3）{'一致' if ok else '不一致'}"); bad += (not ok) + nonan(b, "") + check_orb(b, "", TH0["new"])
    print("== newdefault（rc=0）=="); b = head(os.path.join(base, "dom_v9q4_newdefault_rc0.html"))
    if b is None: print("  読めない"); bad += 1
    else:
        ok = "辺 6 本" in b[1] and all(f"{l}：Δ = " in b[1] for l in ["12", "23", "13", "14", "24", "34"]); print(f"  辺 6 本・対 6 組 {'一致' if ok else '不一致'}"); bad += (not ok) + nonan(b, "")
    print("== auto6s（rc=1）=="); b = head(os.path.join(base, "dom_v9q4_auto6s_rc1.html"))
    if b is None: print("  読めない"); bad += 1
    else: bad += nonan(b, "") + check_orb(b, "", TH0["new"])
    print("== auto6s（rc=0）=="); b = head(os.path.join(base, "dom_v9q4_auto6s_rc0.html"))
    if b is None: print("  読めない"); bad += 1
    else: bad += nonan(b, "")
    print("対照 " + ("一致" if bad == 0 else f"不一致 {bad}"))
    return 0 if bad == 0 else 2
if __name__ == "__main__": sys.exit(main())
