#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""HTML の対照（2026-11-07）：右のパワーの切替（rc=0／1）と表示（0・1・2）のすべてで、index_v8q3.html の h1・h2 が原本 index_v7q3.html と完全に一致する
（読める量・左パネル・枠は変えていない。変わるのは右パネルのパワーの曲線だけで、それは DOM の文字には出ない。描画は figs で確認）。
URL 引数なし＋rc=1 でも h1 に新しい初期値が出ていること。"""
import re, os, sys, html
base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
def head(path):
    s = open(path, encoding="utf-8").read()
    h1 = re.search(r'id="h1">(.*?)</div>', s, re.S); h2 = re.search(r'id="h2"[^>]*>(.*?)</div>', s, re.S)
    return None if not h1 or not h2 else (html.unescape(h1.group(1)), html.unescape(h2.group(1)))
def main():
    bad = 0
    names = sorted(set(re.sub(r"^dom_v7q3_pow\d_rc\d_", "", f)[:-5] for f in os.listdir(base) if f.startswith("dom_v7q3_")))
    for name in names:
        print(f"== {name} ==")
        for pw in (0, 1, 2):
            for rc in (0, 1):
                a = head(os.path.join(base, f"dom_v7q3_pow{pw}_rc{rc}_{name}.html")); b = head(os.path.join(base, f"dom_v8q3_pow{pw}_rc{rc}_{name}.html"))
                if a is None or b is None: print(f"  pow={pw} rc={rc} 読めない"); bad += 1; continue
                ok = (a == b); print(f"  表示 {pw}、右のパワー {'座標' if rc else '円周上'}：原本と h1・h2 {'一致' if ok else '不一致'}"); bad += (not ok)
    b = head(os.path.join(base, "dom_v8q3_newdefault_rc1.html"))
    ok = b is not None and "1: θ = 97°, φ = 1°, δ = 20°/t, M = 100" in b[0] and "表示 ×1 = " in b[0]
    print(f"== newdefault（URL 引数なし、rc=1）== 初期値 {'一致' if ok else '不一致'}"); bad += (not ok)
    print("対照 " + ("一致" if bad == 0 else f"不一致 {bad}"))
    return 0 if bad == 0 else 2
if __name__ == "__main__": sys.exit(main())
