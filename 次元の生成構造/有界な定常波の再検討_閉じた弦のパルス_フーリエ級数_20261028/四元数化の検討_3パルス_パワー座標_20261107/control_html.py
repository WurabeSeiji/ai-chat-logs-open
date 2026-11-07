#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""HTML の対照（2026-11-07）：(1) 右のパワー＝円周上（rc=0）では、index_v8q3.html の h1・h2 が原本 index_v7q3.html と完全一致。
(2) 座標（rc=1）では h2 が一致し、h1 は「　軌道：」の項目を除いて一致。波ごとの点 (i_i Θ₁, j_i Θ₂, k_i Θ₃) を h1 の n_i と Θ から検算し、Θ_i = θ_i(初期) + ティック × φ_i、Θ_i ≡ θ_i (mod 360)、右の枠 = max(360, max|点|) を検算。
(3) Auto を仮想時間 6 秒走らせた後、軌道の点数がティック数 + 1 になっている（ティックごとに 1 点）。"""
import re, os, sys, html, math
base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
ORB = r"　軌道（i_i Θ₁, j_i Θ₂, k_i Θ₃）：Θ = \((-?\d+\.\d), (-?\d+\.\d), (-?\d+\.\d)\)、1: \((-?\d+\.\d), (-?\d+\.\d), (-?\d+\.\d)\)、2: \((-?\d+\.\d), (-?\d+\.\d), (-?\d+\.\d)\)、3: \((-?\d+\.\d), (-?\d+\.\d), (-?\d+\.\d)\)、各 (\d+) 点、右の枠 ±(\d+)°"
def head(path):
    s = open(path, encoding="utf-8").read()
    h1 = re.search(r'id="h1">(.*?)</div>', s, re.S); h2 = re.search(r'id="h2"[^>]*>(.*?)</div>', s, re.S)
    return None if not h1 or not h2 else (html.unescape(h1.group(1)), html.unescape(h2.group(1)))
def strip_orb(h1):
    m = re.search(ORB, h1); return (h1[:m.start()] + h1[m.end():], m) if m else (h1, None)
def check_orb(b, tag, th0):
    h1s, m = strip_orb(b[0])
    if not m: print(f"  {tag}軌道の項目がない"); return 1
    ths = [float(t) for t in re.findall(r"\d: θ = (-?\d+)°", b[0])]; phis = [float(t) for t in re.findall(r"\d: θ = -?\d+°, φ = (\d+)°", b[0])]
    nvs = [[float(x) for x in t] for t in re.findall(r"n = \((-?\d+\.\d+), (-?\d+\.\d+), (-?\d+\.\d+)\)", b[0])] if "q" == "q" else [[1.0, 0.0, 0.0]] * 3
    Th = [float(m.group(a + 1)) for a in range(3)]; pts = [[float(m.group(3 + 3 * i + a + 1)) for a in range(3)] for i in range(3)]; npts = int(m.group(13)); L2 = int(m.group(14))
    t = re.search(r"ティック (\d+)", b[1]); tk = int(t.group(1)) if t else None
    ok1 = len(nvs) == 3 and all(abs(pts[i][a] - nvs[i][a] * Th[a]) < 0.06 for i in range(3) for a in range(3))
    ok2 = tk is not None and len(phis) == 3 and all(abs(Th[i] - (th0[i] + tk * phis[i])) < 0.06 for i in range(3))
    ok3 = len(ths) == 3 and all(min(abs((Th[i] - ths[i]) % 360), abs(360 - (Th[i] - ths[i]) % 360)) < 0.06 for i in range(3))
    ok4 = L2 == round(max(360.0, max(abs(v) for q in pts for v in q)))
    okn = tk is not None and npts == tk + 1
    print(f"  {tag}Θ = ({Th[0]:.1f}, {Th[1]:.1f}, {Th[2]:.1f})　波ごとの点 " + "、".join(f"{i + 1}: ({pts[i][0]:.1f}, {pts[i][1]:.1f}, {pts[i][2]:.1f})" for i in range(3)) + f" = n_i ∘ Θ：{'一致' if ok1 else '不一致'}　Θ_i = θ_i(初期) + ティック × φ_i：{'一致' if ok2 else '不一致'}　Θ_i ≡ θ_i (mod 360)：{'一致' if ok3 else '不一致'}　右の枠 ±{L2}° = max(360, max|点|)：{'一致' if ok4 else '不一致'}　点数 {npts}（ティック {tk if tk is not None else '?'} + 1）{'一致' if okn else '不一致'}")
    return (not ok1) + (not ok2) + (not ok3) + (not ok4) + (not okn)
TH0 = {"default": [180, 180, 180], "tri": [180, 150, 100], "new": [97, 194, 58]}   # 各配置の θ_i の初期値（shoot.sh の URL と既定値）
def main():
    bad = 0
    for name in ["default", "tri"]:
        print(f"== {name} ==")
        for pw in (0, 1, 2):
            for rc in (0, 1):
                a = head(os.path.join(base, f"dom_v7q3_pow{pw}_rc{rc}_{name}.html")); b = head(os.path.join(base, f"dom_v8q3_pow{pw}_rc{rc}_{name}.html"))
                if a is None or b is None: print(f"  pow={pw} rc={rc} 読めない"); bad += 1; continue
                if rc == 0:
                    ok = (a == b); print(f"  表示 {pw}、円周上：原本と h1・h2 {'一致' if ok else '不一致'}"); bad += (not ok)
                else:
                    h1s, m = strip_orb(b[0]); ok = (h1s == a[0] and b[1] == a[1]); print(f"  表示 {pw}、座標：h1（軌道の項目を除く）・h2 {'一致' if ok else '不一致'}"); bad += (not ok)
                    bad += check_orb(b, "", TH0[name])
    print("== auto6s（rc=1、既定 軸 i, j, k）=="); b = head(os.path.join(base, "dom_v8q3_auto6s_rc1.html"))
    if b is None: print("  読めない"); bad += 1
    else: bad += check_orb(b, "", TH0["new"])
    print("== auto6s_kkk（rc=1、軸 全部 k）=="); b = head(os.path.join(base, "dom_v8q3_auto6s_kkk_rc1.html"))
    if b is None: print("  読めない"); bad += 1
    else: bad += check_orb(b, "", TH0["new"])
    print("== newdefault（rc=1）=="); b = head(os.path.join(base, "dom_v8q3_newdefault_rc1.html"))
    ok = b is not None and "1: θ = 97°, φ = 1°, δ = 20°/t, M = 100" in b[0] and "A = 1.0 n = (1.000, 0.000, 0.000) R1" in b[0] and "A = 1.0 n = (0.000, 1.000, 0.000) R2" in b[0] and "A = 1.0 n = (0.000, 0.000, 1.000) R3" in b[0]; print(f"  初期値（θ・φ・δ・M、n = i, j, k）{'一致' if ok else '不一致'}"); bad += (not ok)
    if b is not None: bad += check_orb(b, "", TH0["new"])
    print("対照 " + ("一致" if bad == 0 else f"不一致 {bad}"))
    return 0 if bad == 0 else 2
if __name__ == "__main__": sys.exit(main())
