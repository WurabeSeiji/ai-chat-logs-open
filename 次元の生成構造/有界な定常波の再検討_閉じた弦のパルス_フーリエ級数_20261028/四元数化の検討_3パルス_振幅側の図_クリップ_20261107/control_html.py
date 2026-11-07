#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""HTML の対照（2026-11-07）：表示 3 択のそれぞれ（4 配置＋Auto 40 ティック）で、index_v7q3.html の h1 は「　右：」より前が原本 index_v5q3.html と
完全一致（平行移動とクリップの項目はその後ろ。平行移動 = −kR(cos θ_i, sin θ_i) も検算）、h2 は完全一致（左パネル・平行移動・計算は変えていない）。クリップの範囲（各パルスの山、|Z_i|² ≥ 最大の 5%）が check_clip_v7q3.py の独立な再計算と一致。"""
import re, os, sys, html
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_clip_v7q3 as chk
base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")

def head(path):
    s = open(path, encoding="utf-8").read()
    h1 = re.search(r'id="h1">(.*?)</div>', s, re.S); h2 = re.search(r'id="h2"[^>]*>(.*?)</div>', s, re.S)
    return None if not h1 or not h2 else (html.unescape(h1.group(1)), html.unescape(h2.group(1)))

def check(name, pw):
    a = head(os.path.join(base, f"dom_v5q3_pow{pw}_{name}.html")); b = head(os.path.join(base, f"dom_v7q3_pow{pw}_{name}.html"))
    if a is None or b is None: print(f"  pow={pw} 読めない"); return 1
    k = b[0].find("　右："); ok1 = (k > 0 and b[0][:k] == a[0]); ok2 = (b[1] == a[1])   # v7q3 は v5q3 に「右：平行移動」と「クリップ」の両方を足したので「右：」より前を比べる
    import math
    mk = re.search(r"表示 ×\d+ = (\d+\.\d+)", b[0]); sh = re.findall(r"(\d): \((-?\d+\.\d+), (-?\d+\.\d+)\)", b[0]); th0 = [float(t) for t in re.findall(r"\d: θ = (-?\d+)°", b[0])]
    okd = bool(mk) and len(sh) == 3 and len(th0) == 3 and all(abs(float(x) + float(mk.group(1)) * math.cos(math.radians(th0[int(i) - 1]))) < 2e-3 and abs(float(y) + float(mk.group(1)) * math.sin(math.radians(th0[int(i) - 1]))) < 2e-3 for i, x, y in sh)
    ths = re.findall(r"(\d): θ = (-?\d+)°, φ = (\d+)°, δ = -?\d+°/t, M = (\d+)", b[0]); cl = re.findall(r"(\d): (-?\d+\.\d+)°〜\+(-?\d+\.\d+)°", b[0])
    okc = len(ths) == 3 and len(cl) == 3 and all(abs(float(x) - chk.clip(float(t), float(f), M=int(mm))[0]) < 0.3 and abs(float(y) - chk.clip(float(t), float(f), M=int(mm))[1]) < 0.3 for (i, t, f, mm), (j, x, y) in zip(ths, cl))
    print(f"  表示 {['Re・Im', '合計', '分離'][pw]}：h1（右より前）{'一致' if ok1 else '不一致'}　h2 {'一致' if ok2 else '不一致'}　平行移動 {'一致' if okd else '不一致'}　山の範囲 " + "、".join(f"{i}: {x}°〜+{y}°" for i, x, y in cl) + f" → Python 再計算と {'一致' if okc else '不一致'}")
    return (not ok1) + (not ok2) + (not okd) + (not okc)

def main():
    bad = 0
    for name in ["default", "tri", "tri_ijk", "mirror", "mirror_k_i_mi", "std"]:
        print(f"== {name} =="); bad += sum(check(name, pw) for pw in (0, 1, 2))
    print("== auto6s =="); bad += check("auto6s", 0)
    print("== newdefault（v7q3 の新しい初期値、URL 引数なし）==")
    b = head(os.path.join(base, "dom_v7q3_newdefault.html"))
    if b is None: print("  読めない"); bad += 1
    else:
        ok = ("1: θ = 97°, φ = 1°, δ = 20°/t, M = 100" in b[0] and "2: θ = 194°, φ = 1°, δ = 0°/t, M = 100" in b[0] and "3: θ = 58°, φ = 1°, δ = 10°/t, M = 100" in b[0] and "表示 ×1 = " in b[0] and "パワー合計：" in b[0] and b[0].count("n = (0.000, 0.000, 1.000)") == 3)
        cl = re.findall(r"(\d): (-?\d+\.\d+)°〜\+(-?\d+\.\d+)°", b[0])
        print(f"  初期値 θ = 97, 194, 58、φ = 1、M = 100、δ = 20, 0, 10、R× = 1、パワー合計：{'一致' if ok else '不一致'}　山の範囲 " + "、".join(f"{i}: {x}°〜+{y}°" for i, x, y in cl)); bad += (not ok)
    print("対照 " + ("一致" if bad == 0 else f"不一致 {bad}"))
    return 0 if bad == 0 else 2
if __name__ == "__main__": sys.exit(main())
