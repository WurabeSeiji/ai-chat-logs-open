#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""HTML の対照（2026-11-07）：(1) 4 成分（pow なし）は原本 index_v4q3.html と h1・h2 が完全一致。
(2) パワー分離（pow=2）は原本の pow=1（軸ごと）と h1・h2 が完全一致（項目名「パワー表示」→「パワー分離」だけ揃える）。
(3) パワー合計（pow=1）は h2 が原本と一致、h1 は合計の項目より前が一致、⟨|Q|²⟩_ψ が分離の合計に一致。"""
import re, os, sys, html
base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
NUM = r"(\d+\.\d+)"
SUM1 = r"　パワー合計：max \|Q\|² = " + NUM + r"　⟨\|Q\|²⟩_ψ = " + NUM
SPLIT_TOT = r"　合計 ⟨\|Q\|²⟩_ψ = " + NUM

def head(path):
    s = open(path, encoding="utf-8").read()
    h1 = re.search(r'id="h1">(.*?)</div>', s, re.S); h2 = re.search(r'id="h2"[^>]*>(.*?)</div>', s, re.S)
    return None if not h1 or not h2 else (html.unescape(h1.group(1)), html.unescape(h2.group(1)))

def check_sum(b1, b2, ref, tag):
    bad = 0
    m = re.search(SUM1, b1[0]); mt = re.search(SPLIT_TOT, b2[0])
    if not m or not mt: print(f"  {tag}合計か分離の項目がない"); return 1
    tmax, tmean = map(float, m.groups()); tot2 = float(mt.group(1))
    ok_tot = abs(tmean - tot2) < 2e-4
    ok_pre = (b1[0][:m.start()] == ref[0][:ref[0].index("　枠 ±")]); ok_h2 = (b1[1] == ref[1])
    print(f"  {tag}パワー合計：max |Q|² {tmax:.4f}　⟨|Q|²⟩_ψ {tmean:.4f}（分離の合計 {tot2:.4f} と {'一致' if ok_tot else '不一致'}）　h1（前）{'一致' if ok_pre else '不一致'}　h2 {'一致' if ok_h2 else '不一致'}")
    return (not ok_tot) + (not ok_pre) + (not ok_h2)

def main():
    names = sorted(set(f[9:-5] for f in os.listdir(base) if f.startswith("dom_v4q3_") and f.endswith(".html") and "pow1" not in f))
    bad = 0
    for name in names:
        print(f"== {name} ==")
        a = head(os.path.join(base, f"dom_v4q3_{name}.html")); a1 = head(os.path.join(base, f"dom_v4q3_pow1_{name}.html"))
        b = head(os.path.join(base, f"dom_v5q3_{name}.html")); b1 = head(os.path.join(base, f"dom_v5q3_pow1_{name}.html")); b2 = head(os.path.join(base, f"dom_v5q3_pow2_{name}.html"))
        if None in (a, a1, b, b1, b2): print("  読めない"); bad += 1; continue
        ok = (a == b); print(f"  4 成分：原本と h1・h2 {'一致' if ok else '不一致'}"); bad += (not ok)
        ok = (a1[0].replace("パワー表示（軸ごと", "パワー分離（軸ごと") == b2[0] and a1[1] == b2[1]); print(f"  パワー分離：原本 pow=1（軸ごと）と h1・h2 {'一致' if ok else '不一致'}"); bad += (not ok)
        bad += check_sum(b1, b2, a, "")
    print("== auto6s ==")
    a1 = head(os.path.join(base, "dom_v4q3_pow1_auto6s.html")); b1 = head(os.path.join(base, "dom_v5q3_pow1_auto6s.html")); b2 = head(os.path.join(base, "dom_v5q3_pow2_auto6s.html"))
    if None in (a1, b1, b2): print("  読めない"); bad += 1
    else:
        t = re.search(r"ティック (\d+)", b2[1]); print(f"  ティック {t.group(1) if t else '?'}")
        ok = (a1[0].replace("パワー表示（軸ごと", "パワー分離（軸ごと") == b2[0] and a1[1] == b2[1]); print(f"  パワー分離：原本 pow=1 と {'一致' if ok else '不一致'}"); bad += (not ok)
        m = re.search(SUM1, b1[0]); mt = re.search(SPLIT_TOT, b2[0])
        if m and mt:
            tmean = float(m.group(2)); tot2 = float(mt.group(1)); ok = abs(tmean - tot2) < 2e-4
            print(f"  パワー合計：⟨|Q|²⟩_ψ {tmean:.4f}（分離の合計 {tot2:.4f} と {'一致' if ok else '不一致'}）"); bad += (not ok)
        else: print("  項目がない"); bad += 1
    print("対照 " + ("一致" if bad == 0 else f"不一致 {bad}"))
    return 0 if bad == 0 else 2
if __name__ == "__main__": sys.exit(main())
