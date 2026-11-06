#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""HTML の対照（2026-11-07）：shoot.sh が取った DOM から右上の読める量（#h1・#h2 の文字列）を取り出し、
原本 index_v2q3.html と index_v3q3.html が同じ URL 引数で一致することを確かめる。v3q3 の h1 には「枠 ±L（いま必要な枠 ±L′…）」が
末尾に足されるので、その部分を除いて比べる。さらに Auto を仮想時間 6 秒走らせた後の v3q3 で、枠 L が auto=0 の値のまま
（表示スケールが変わっていない）で、いま必要な枠 L′ だけが変わっていることを確かめる。"""
import re, os, sys, html

base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
FRAME = r"　枠 ±(\d+\.\d+)（いま必要な枠 ±(\d+\.\d+)(、はみ出し)?）"

def head(path):
    s = open(path, encoding="utf-8").read()
    h1 = re.search(r'id="h1">(.*?)</div>', s, re.S); h2 = re.search(r'id="h2"[^>]*>(.*?)</div>', s, re.S)
    if not h1 or not h2: return None
    return html.unescape(h1.group(1)), html.unescape(h2.group(1))

def main():
    names = sorted(set(f[9:-5] for f in os.listdir(base) if f.startswith("dom_v2q3_") and f.endswith(".html")))
    if not names: print("DOM がない。先に sh shoot.sh"); return 1
    bad = 0; L0 = None
    for name in names:
        c = head(os.path.join(base, f"dom_v2q3_{name}.html")); v = head(os.path.join(base, f"dom_v3q3_{name}.html"))
        print(f"== {name} ==")
        if c is None or v is None: print("  読めない（Plotly か DOM の取得に失敗）"); bad += 1; continue
        m = re.search(FRAME, v[0])
        if not m: print("  v3q3 の h1 に枠の項目がない"); bad += 1; continue
        h1v = v[0][:m.start()]
        ok1 = (h1v == c[0]); ok2 = (v[1] == c[1])
        print(f"  h1（枠の項目を除く）: {'一致' if ok1 else '不一致'}　h2: {'一致' if ok2 else '不一致'}　枠 ±{m.group(1)}、いま必要な枠 ±{m.group(2)}{m.group(3) or ''}")
        if not ok1: print(f"    原本 h1: {c[0][:200]}\n    v3q3 h1: {h1v[:200]}")
        if not ok2: print(f"    原本 h2: {c[1][:200]}\n    v3q3 h2: {v[1][:200]}")
        if name == "default": L0 = m.group(1)
        bad += (not ok1) + (not ok2)
    a = head(os.path.join(base, "dom_v3q3_auto6s.html"))
    print("== auto6s（v3q3、Auto を仮想時間 6 秒）==")
    if a is None: print("  読めない"); bad += 1
    else:
        m = re.search(FRAME, a[0]); t = re.search(r"ティック (\d+)", a[1])
        if not m: print("  枠の項目がない"); bad += 1
        else:
            fixed = (L0 is not None and m.group(1) == L0)
            print(f"  ティック {t.group(1) if t else '?'}　枠 ±{m.group(1)}（auto=0 の枠 ±{L0}）→ {'固定されている' if fixed else '変わっている'}　いま必要な枠 ±{m.group(2)}{m.group(3) or ''}")
            bad += (not fixed)
    print("対照 " + ("一致" if bad == 0 else f"不一致 {bad}"))
    return 0 if bad == 0 else 2

if __name__ == "__main__":
    sys.exit(main())
