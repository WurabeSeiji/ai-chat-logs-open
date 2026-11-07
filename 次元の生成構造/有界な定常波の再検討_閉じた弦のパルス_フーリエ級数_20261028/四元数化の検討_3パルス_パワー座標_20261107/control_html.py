#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""HTML の対照（2026-11-07）：(1) 表示 0・1・2 は index_v8q3.html と原本 index_v7q3.html の h1・h2 が完全一致（6 配置）。
(2) 表示 3（パワー座標）は h2 が原本と一致、h1 は「パワー座標：」より前が原本（同じ配置の pow=1 の h1 の「パワー合計：」より前）と一致し、
max |Σ P_i n_i| が check_pn_v8q3.py の独立な再計算と一致（5e-4）。新しい初期値（軸 k、軸 i,j,k）も同じ。"""
import re, os, sys, html
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_pn_v8q3 as chk
base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")

def head(path):
    s = open(path, encoding="utf-8").read()
    h1 = re.search(r'id="h1">(.*?)</div>', s, re.S); h2 = re.search(r'id="h2"[^>]*>(.*?)</div>', s, re.S)
    return None if not h1 or not h2 else (html.unescape(h1.group(1)), html.unescape(h2.group(1)))

def check3(b, ref, name):
    m = re.search(r"　パワー座標：P_i = √\|Q_i\|² を n_i の向きに、max \|Σ P_i n_i\| = (\d+\.\d+)", b[0])
    if not m: print("  パワー座標の項目がない"); return 1
    smax = float(m.group(1)); exp = chk.smax(chk.CASES[name]); ok = abs(smax - exp) < 5e-4
    bad = (not ok)
    msg = f"  パワー座標：max |Σ P_i n_i| = {smax:.4f}（Python 再計算 {exp:.4f}）{'一致' if ok else '不一致'}"
    if ref is not None:
        pre_b = b[0][:m.start()]; k = ref[0].find("　パワー合計："); pre_r = ref[0][:k] if k > 0 else None
        ok1 = (pre_r is not None and pre_b == pre_r); ok2 = (b[1] == ref[1]); bad += (not ok1) + (not ok2)
        msg += f"　h1（前）{'一致' if ok1 else '不一致'}　h2 {'一致' if ok2 else '不一致'}"
    print(msg); return bad

def main():
    bad = 0
    for name in ["default", "tri", "tri_ijk", "mirror", "mirror_k_i_mi", "std"]:
        print(f"== {name} ==")
        for pw in (0, 1, 2):
            a = head(os.path.join(base, f"dom_v7q3_pow{pw}_{name}.html")); b = head(os.path.join(base, f"dom_v8q3_pow{pw}_{name}.html"))
            if a is None or b is None: print(f"  pow={pw} 読めない"); bad += 1; continue
            ok = (a == b); print(f"  表示 {['4 成分', '合計', '分離'][pw]}：原本と h1・h2 {'一致' if ok else '不一致'}"); bad += (not ok)
        b3 = head(os.path.join(base, f"dom_v8q3_pow3_{name}.html")); a1 = head(os.path.join(base, f"dom_v7q3_pow1_{name}.html"))
        if b3 is None or a1 is None: print("  pow=3 読めない"); bad += 1
        else: bad += check3(b3, a1, name)
    for name in ["newdefault", "newdefault_ijk"]:
        print(f"== {name}（URL 引数なしの初期値）=="); b3 = head(os.path.join(base, f"dom_v8q3_pow3_{name}.html"))
        if b3 is None: print("  読めない"); bad += 1
        else: bad += check3(b3, None, name)
    print("対照 " + ("一致" if bad == 0 else f"不一致 {bad}"))
    return 0 if bad == 0 else 2
if __name__ == "__main__": sys.exit(main())
