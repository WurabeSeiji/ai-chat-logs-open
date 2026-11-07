#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""HTML の対照（2026-11-07）：表示 3 択のそれぞれ（＋観測者 1 の例、Auto 40 ティック）で、index_v6c3.html の h1 は「　右：」より前が原本 index_v5c3.html と
完全一致、h2 は完全一致（左パネルと計算は変えていない）。右の平行移動 (dx, dy) は −kR(cos θ_o, sin θ_o)：h1 の kR（表示 ×k = …）と観測者の θ から検算。"""
import re, os, sys, html, math
base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")

def head(path):
    s = open(path, encoding="utf-8").read()
    h1 = re.search(r'id="h1">(.*?)</div>', s, re.S); h2 = re.search(r'id="h2"[^>]*>(.*?)</div>', s, re.S)
    return None if not h1 or not h2 else (html.unescape(h1.group(1)), html.unescape(h2.group(1)))

def check(name, pw):
    a = head(os.path.join(base, f"dom_v5c3_pow{pw}_{name}.html")); b = head(os.path.join(base, f"dom_v6c3_pow{pw}_{name}.html"))
    if a is None or b is None: print(f"  pow={pw} 読めない"); return 1
    k = b[0].find("　右："); ok1 = (k > 0 and b[0][:k] == a[0]); ok2 = (b[1] == a[1])
    m = re.search(r"右：観測者 (\d) の中心位相 θ = (-?\d+)° を原点に（平行移動 \((-?\d+\.\d+), (-?\d+\.\d+)\)", b[0]); mk = re.search(r"表示 ×\d+ = (\d+\.\d+)", b[0])
    okd = False; txt = ""
    if m and mk:
        th = float(m.group(2)); Rd = float(mk.group(1)); dx, dy = float(m.group(3)), float(m.group(4))
        okd = abs(dx + Rd * math.cos(math.radians(th))) < 2e-3 and abs(dy + Rd * math.sin(math.radians(th))) < 2e-3
        txt = f"観測者 {m.group(1)}、θ = {th:.0f}°、kR = {Rd:.3f}、平行移動 ({dx:.3f}, {dy:.3f})"
    print(f"  表示 {['Re・Im', '合計', '分離'][pw]}：h1（右より前）{'一致' if ok1 else '不一致'}　h2 {'一致' if ok2 else '不一致'}　{txt} → {'一致' if okd else '不一致'}")
    return (not ok1) + (not ok2) + (not okd)

def main():
    bad = 0
    for name in ["default", "tri", "tri_obs1", "mirror", "std"]:
        print(f"== {name} =="); bad += sum(check(name, pw) for pw in (0, 1, 2))
    print("== auto6s =="); bad += check("auto6s", 0)
    print("対照 " + ("一致" if bad == 0 else f"不一致 {bad}"))
    return 0 if bad == 0 else 2
if __name__ == "__main__": sys.exit(main())
