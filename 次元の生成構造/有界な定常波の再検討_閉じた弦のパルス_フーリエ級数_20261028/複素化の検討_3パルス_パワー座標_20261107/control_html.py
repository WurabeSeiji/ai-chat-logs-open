#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""HTML の対照（2026-11-07）：(1) 右のパワー＝円周上（rc=0）では、index_v8c3.html の h1・h2 が原本 index_v7c3.html と完全一致。
(2) 座標（rc=1）では h2 が一致し、h1 は「　軌道：」の項目を除いて一致。軌道の点 (kR sin θ₁, kR sin θ₂, kR sin θ₃) を h1 の kR と θ_i から検算。
(3) Auto を仮想時間 6 秒走らせた後、軌道の点数がティック数 + 1 になっている（ティックごとに 1 点）。"""
import re, os, sys, html, math
base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
ORB = r"　軌道：\(kR sin θ₁, kR sin θ₂, kR sin θ₃\) = \((-?\d+\.\d+), (-?\d+\.\d+), (-?\d+\.\d+)\)、(\d+) 点"
def head(path):
    s = open(path, encoding="utf-8").read()
    h1 = re.search(r'id="h1">(.*?)</div>', s, re.S); h2 = re.search(r'id="h2"[^>]*>(.*?)</div>', s, re.S)
    return None if not h1 or not h2 else (html.unescape(h1.group(1)), html.unescape(h2.group(1)))
def strip_orb(h1):
    m = re.search(ORB, h1); return (h1[:m.start()] + h1[m.end():], m) if m else (h1, None)
def check_orb(b, tag):
    h1s, m = strip_orb(b[0])
    if not m: print(f"  {tag}軌道の項目がない"); return 1
    ths = [float(t) for t in re.findall(r"\d: θ = (-?\d+)°", b[0])]; mk = re.search(r"表示 ×\d+ = (\d+\.\d+)", b[0])
    pt = [float(m.group(i)) for i in (1, 2, 3)]; npts = int(m.group(4))
    ok = bool(mk) and len(ths) == 3 and all(abs(pt[i] - float(mk.group(1)) * math.sin(math.radians(ths[i]))) < 2e-3 for i in range(3))
    t = re.search(r"ティック (\d+)", b[1]); okn = t is not None and npts == int(t.group(1)) + 1
    print(f"  {tag}軌道の点 ({pt[0]:.4f}, {pt[1]:.4f}, {pt[2]:.4f}) = kR sin θ：{'一致' if ok else '不一致'}　点数 {npts}（ティック {t.group(1) if t else '?'} + 1）{'一致' if okn else '不一致'}")
    return (not ok) + (not okn)
def main():
    bad = 0
    for name in ["default", "tri"]:
        print(f"== {name} ==")
        for pw in (0, 1, 2):
            for rc in (0, 1):
                a = head(os.path.join(base, f"dom_v7c3_pow{pw}_rc{rc}_{name}.html")); b = head(os.path.join(base, f"dom_v8c3_pow{pw}_rc{rc}_{name}.html"))
                if a is None or b is None: print(f"  pow={pw} rc={rc} 読めない"); bad += 1; continue
                if rc == 0:
                    ok = (a == b); print(f"  表示 {pw}、円周上：原本と h1・h2 {'一致' if ok else '不一致'}"); bad += (not ok)
                else:
                    h1s, m = strip_orb(b[0]); ok = (h1s == a[0] and b[1] == a[1]); print(f"  表示 {pw}、座標：h1（軌道の項目を除く）・h2 {'一致' if ok else '不一致'}"); bad += (not ok)
                    bad += check_orb(b, "")
    print("== auto6s（rc=1）=="); b = head(os.path.join(base, "dom_v8c3_auto6s_rc1.html"))
    if b is None: print("  読めない"); bad += 1
    else: bad += check_orb(b, "")
    print("== newdefault（rc=1）=="); b = head(os.path.join(base, "dom_v8c3_newdefault_rc1.html"))
    ok = b is not None and "1: θ = 97°, φ = 1°, δ = 20°/t, M = 100" in b[0]; print(f"  初期値 {'一致' if ok else '不一致'}"); bad += (not ok)
    if b is not None: bad += check_orb(b, "")
    print("対照 " + ("一致" if bad == 0 else f"不一致 {bad}"))
    return 0 if bad == 0 else 2
if __name__ == "__main__": sys.exit(main())
