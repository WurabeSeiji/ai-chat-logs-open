#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""HTML の対照（2026-11-07）：(1) Re・Im（pow なし）と パワー合計（pow=1）は原本 index_v4c3.html と h1・h2 が完全一致
（合計の項目名だけ「パワー表示」→「パワー合計」に変えたので、その語を揃えて比べる）。
(2) パワー分離（pow=2）は h2 が原本と一致、h1 は分離の項目より前が一致、Re² と Im² の ψ 平均の和が合計に一致し、合計は pow=1 の ⟨|Z|²⟩_ψ に一致。"""
import re, os, sys, html
base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
NUM = r"(\d+\.\d+)"
SPLIT = r"　パワー分離（max／ψ 平均）：Re² " + NUM + "／" + NUM + "　Im² " + NUM + "／" + NUM + r"　合計 ⟨\|Z\|²⟩_ψ = " + NUM + r"　\|C₀\|² \+ R² = " + NUM
SUM1 = r"⟨\|Z\|²⟩_ψ = " + NUM

def head(path):
    s = open(path, encoding="utf-8").read()
    h1 = re.search(r'id="h1">(.*?)</div>', s, re.S); h2 = re.search(r'id="h2"[^>]*>(.*?)</div>', s, re.S)
    return None if not h1 or not h2 else (html.unescape(h1.group(1)), html.unescape(h2.group(1)))

def main():
    names = sorted(set(f[9:-5] for f in os.listdir(base) if f.startswith("dom_v4c3_") and f.endswith(".html") and "pow1" not in f))
    bad = 0
    for name in names:
        print(f"== {name} ==")
        a = head(os.path.join(base, f"dom_v4c3_{name}.html")); b = head(os.path.join(base, f"dom_v5c3_{name}.html"))
        a1 = head(os.path.join(base, f"dom_v4c3_pow1_{name}.html")); b1 = head(os.path.join(base, f"dom_v5c3_pow1_{name}.html")); b2 = head(os.path.join(base, f"dom_v5c3_pow2_{name}.html"))
        if None in (a, b, a1, b1, b2): print("  読めない"); bad += 1; continue
        ok = (a == b); print(f"  Re・Im：原本と h1・h2 {'一致' if ok else '不一致'}"); bad += (not ok)
        ok = (a1[0].replace("パワー表示：", "パワー合計：") == b1[0] and a1[1] == b1[1]); print(f"  パワー合計：原本 pow=1 と h1・h2 {'一致' if ok else '不一致'}"); bad += (not ok)
        m = re.search(SPLIT, b2[0]); m1 = re.search(SUM1, b1[0])
        if not m or not m1: print("  分離の項目がない"); bad += 1; continue
        g = list(map(float, m.groups())); tot1 = float(m1.group(1))
        ok_sum = abs(g[1] + g[3] - g[4]) < 2e-4; ok_tot = abs(g[4] - tot1) < 2e-4
        ok_pre = (b2[0][:m.start()] == a[0][:a[0].index("　枠 ±")]); ok_h2 = (b2[1] == a[1])
        print(f"  パワー分離：Re² {g[0]:.4f}／{g[1]:.4f}　Im² {g[2]:.4f}／{g[3]:.4f}　合計 {g[4]:.4f}（和 {'一致' if ok_sum else '不一致'}、合計モードの ⟨|Z|²⟩ {tot1:.4f} と {'一致' if ok_tot else '不一致'}）　h1（前）{'一致' if ok_pre else '不一致'}　h2 {'一致' if ok_h2 else '不一致'}")
        bad += (not ok_sum) + (not ok_tot) + (not ok_pre) + (not ok_h2)
    print("== auto6s ==")
    a1 = head(os.path.join(base, "dom_v4c3_pow1_auto6s.html")); b1 = head(os.path.join(base, "dom_v5c3_pow1_auto6s.html")); b2 = head(os.path.join(base, "dom_v5c3_pow2_auto6s.html"))
    if None in (a1, b1, b2): print("  読めない"); bad += 1
    else:
        ok = (a1[0].replace("パワー表示：", "パワー合計：") == b1[0] and a1[1] == b1[1]); print(f"  パワー合計：原本 pow=1 と {'一致' if ok else '不一致'}"); bad += (not ok)
        m = re.search(SPLIT, b2[0]); m1 = re.search(SUM1, b1[0]); t = re.search(r"ティック (\d+)", b2[1])
        if m and m1:
            g = list(map(float, m.groups())); tot1 = float(m1.group(1)); ok_sum = abs(g[1] + g[3] - g[4]) < 2e-4; ok_tot = abs(g[4] - tot1) < 2e-4
            print(f"  ティック {t.group(1) if t else '?'}　パワー分離：Re² {g[0]:.4f}／{g[1]:.4f}　Im² {g[2]:.4f}／{g[3]:.4f}　合計 {g[4]:.4f}（和 {'一致' if ok_sum else '不一致'}、合計モードと {'一致' if ok_tot else '不一致'}）")
            bad += (not ok_sum) + (not ok_tot)
        else: print("  分離の項目がない"); bad += 1
    print("対照 " + ("一致" if bad == 0 else f"不一致 {bad}"))
    return 0 if bad == 0 else 2
if __name__ == "__main__": sys.exit(main())
