#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""HTML の対照（2026-11-07）：表示 3 択のそれぞれ（4 配置＋Auto 40 ティック）で、index_v7c3.html の h1 は「　クリップ：」より前が原本 index_v6c3.html と
完全一致、h2 は完全一致（左パネル・平行移動・計算は変えていない）。クリップの範囲 θ_i ± φ_i が h1 のパルスの θ・φ と一致。"""
import re, os, sys, html
base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")

def head(path):
    s = open(path, encoding="utf-8").read()
    h1 = re.search(r'id="h1">(.*?)</div>', s, re.S); h2 = re.search(r'id="h2"[^>]*>(.*?)</div>', s, re.S)
    return None if not h1 or not h2 else (html.unescape(h1.group(1)), html.unescape(h2.group(1)))

def check(name, pw):
    a = head(os.path.join(base, f"dom_v6c3_pow{pw}_{name}.html")); b = head(os.path.join(base, f"dom_v7c3_pow{pw}_{name}.html"))
    if a is None or b is None: print(f"  pow={pw} 読めない"); return 1
    k = b[0].find("　クリップ："); ok1 = (k > 0 and b[0][:k] == a[0]); ok2 = (b[1] == a[1])
    ths = re.findall(r"(\d): θ = (-?\d+)°, φ = (\d+)°", b[0]); cl = re.findall(r"(\d): (-?\d+)° ± (\d+)°", b[0])
    okc = len(ths) == 3 and ths == cl
    print(f"  表示 {['Re・Im', '合計', '分離'][pw]}：h1（クリップより前）{'一致' if ok1 else '不一致'}　h2 {'一致' if ok2 else '不一致'}　クリップ " + "、".join(f"{i}: {t}° ± {f}°" for i, t, f in cl) + f" → {'一致' if okc else '不一致'}")
    return (not ok1) + (not ok2) + (not okc)

def main():
    bad = 0
    for name in ["default", "tri", "mirror", "std"]:
        print(f"== {name} =="); bad += sum(check(name, pw) for pw in (0, 1, 2))
    print("== auto6s =="); bad += check("auto6s", 0)
    print("対照 " + ("一致" if bad == 0 else f"不一致 {bad}"))
    return 0 if bad == 0 else 2
if __name__ == "__main__": sys.exit(main())
