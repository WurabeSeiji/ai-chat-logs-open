#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""HTML の対照（2026-11-07）：(1) 表示 3 択のそれぞれで、index_v6c3.html の h1 は「　右枠」より前が原本 index_v5c3.html の h1 と完全一致、h2 は完全一致
（左パネルと読める量は変えていない）。(2) 右を線分にしたとき、|Z|² の山（数・ψ*・P・幅・直径）が check_lobes_v6c3.py の独立な再計算と一致し、最大の山の P は左の max |Z|²（pow=1 の h1）と一致。"""
import re, os, sys, html
base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_lobes_v6c3 as chk

def head(path):
    s = open(path, encoding="utf-8").read()
    h1 = re.search(r'id="h1">(.*?)</div>', s, re.S); h2 = re.search(r'id="h2"[^>]*>(.*?)</div>', s, re.S)
    return None if not h1 or not h2 else (html.unescape(h1.group(1)), html.unescape(h2.group(1)))

def lobes(h1):
    m = re.search(r"　線分：(.*)$", h1)
    if not m: return None
    return [(float(a), float(b), float(c), float(d)) for a, b, c, d in re.findall(r"ψ\* (-?\d+\.\d+)°, P (\d+\.\d+), 幅 (\d+\.\d+)°, 直径 (\d+\.\d+)", m.group(1))]

def main():
    bad = 0
    names = ["default", "tri", "mirror", "std"]
    for name in names:
        print(f"== {name} ==")
        for pw in (0, 1, 2):
            a = head(os.path.join(base, f"dom_v5c3_pow{pw}_{name}.html")); b = head(os.path.join(base, f"dom_v6c3_pow{pw}_{name}.html"))
            if a is None or b is None: print(f"  pow={pw} 読めない"); bad += 1; continue
            k = b[0].find("　右枠"); ok1 = (k > 0 and b[0][:k] == a[0]); ok2 = (b[1] == a[1])
            print(f"  表示 {['Re・Im', '合計', '分離'][pw]}：h1（右枠より前）{'一致' if ok1 else '不一致'}　h2 {'一致' if ok2 else '不一致'}"); bad += (not ok1) + (not ok2)
        sgm = head(os.path.join(base, f"dom_v6c3_seg_{name}.html")); a1 = head(os.path.join(base, f"dom_v5c3_pow1_{name}.html"))
        lb = lobes(sgm[0]) if sgm else None
        if not lb: print("  線分：山が読めない"); bad += 1; continue
        exp = chk.lobes(chk.CASES[name])
        okn = len(lb) == len(exp)
        okv = okn and all(abs(g[0] - e[0]) < 0.3 and abs(g[1] - e[1]) < 2e-4 and abs(g[2] - e[2]) < 0.3 and abs(g[3] - e[3]) < 2e-4 for g, e in zip(sorted(lb), sorted(exp)))
        mx = re.search(r"max \|Z\|² = (\d+\.\d+)", a1[0]); okP = mx is not None and abs(max(l[1] for l in lb) - float(mx.group(1))) < 2e-4
        print(f"  線分：山 {len(lb)} 個（Python 再計算 {len(exp)} 個）{'一致' if okn else '不一致'}　ψ*・P・幅・直径 {'一致' if okv else '不一致'}　最大の P {max(l[1] for l in lb):.4f} = 左の max |Z|² {'一致' if okP else '不一致'}")
        for l in lb: print(f"    ψ* {l[0]:.1f}°, P {l[1]:.4f}, 幅 {l[2]:.1f}°, 直径 {l[3]:.4f}")
        bad += (not okn) + (not okv) + (not okP)
    print("== auto6s ==")
    a = head(os.path.join(base, "dom_v5c3_pow0_auto6s.html")); b = head(os.path.join(base, "dom_v6c3_pow0_auto6s.html")); sgm = head(os.path.join(base, "dom_v6c3_seg_auto6s.html"))
    if None in (a, b, sgm): print("  読めない"); bad += 1
    else:
        k = b[0].find("　右枠"); ok1 = (k > 0 and b[0][:k] == a[0]); ok2 = (b[1] == a[1]); t = re.search(r"ティック (\d+)", b[1])
        print(f"  ティック {t.group(1) if t else '?'}　h1（右枠より前）{'一致' if ok1 else '不一致'}　h2 {'一致' if ok2 else '不一致'}"); bad += (not ok1) + (not ok2)
        lb = lobes(sgm[0]) or []; print("  線分：" + "；".join(f"ψ* {l[0]:.1f}°, P {l[1]:.4f}, 幅 {l[2]:.1f}°, 直径 {l[3]:.4f}" for l in lb))
    print("対照 " + ("一致" if bad == 0 else f"不一致 {bad}"))
    return 0 if bad == 0 else 2
if __name__ == "__main__": sys.exit(main())
