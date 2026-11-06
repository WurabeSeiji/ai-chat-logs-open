#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""HTML の対照（2026-11-07）：shoot.sh が取った DOM から右上の読める量（#h1・#h2）を取り出して比べる。
(1) OFF：index_v4c3.html は原本 index_v3c3.html と h1・h2 が完全に一致する（同じ URL 引数 4 例）。
(2) ON（pow=1）：h2 は原本と一致（読める量は表示モードで変わらない）。h1 は「　パワー表示：」より前が原本の「　枠」より前と一致し、
    パワー表示の項目で ⟨|Z|²⟩_ψ と |C₀|² + R² が一致する（5e-4）。Auto 40 ティック後（内部位相が回った後）も同じ。"""
import re, os, sys, html

base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
POW = r"　パワー表示：max \|Z\|² = (\d+\.\d+)　⟨\|Z\|²⟩_ψ = (\d+\.\d+)　\|C₀\|² \+ R² = (\d+\.\d+)"

def head(path):
    s = open(path, encoding="utf-8").read()
    h1 = re.search(r'id="h1">(.*?)</div>', s, re.S); h2 = re.search(r'id="h2"[^>]*>(.*?)</div>', s, re.S)
    if not h1 or not h2: return None
    return html.unescape(h1.group(1)), html.unescape(h2.group(1))

def check_pow(v, ref_h1=None, ref_h2=None, tag=""):
    """ON の DOM v を検査。ref があれば原本との一致も見る。戻り値は不一致の数"""
    bad = 0
    m = re.search(POW, v[0])
    if not m: print(f"  {tag}パワー表示の項目がない"); return 1
    pmax, pmean, chk = map(float, m.groups())
    ok = abs(pmean - chk) < 5e-4
    print(f"  {tag}max |Z|² = {pmax:.4f}　⟨|Z|²⟩_ψ = {pmean:.4f}　|C₀|² + R² = {chk:.4f} → {'一致' if ok else '不一致'}")
    bad += (not ok)
    if ref_h1 is not None:
        pre_v = v[0][:m.start()]; pre_r = ref_h1[:ref_h1.index("　枠 ±")]
        ok1 = (pre_v == pre_r); ok2 = (v[1] == ref_h2)
        print(f"  {tag}h1（パワー表示より前）: {'一致' if ok1 else '不一致'}　h2: {'一致' if ok2 else '不一致'}")
        if not ok1: print(f"    原本: {pre_r[:200]}\n    ON  : {pre_v[:200]}")
        bad += (not ok1) + (not ok2)
    return bad

def main():
    names = sorted(set(f[9:-5] for f in os.listdir(base) if f.startswith("dom_v3c3_") and f.endswith(".html")))
    if not names: print("DOM がない。先に sh shoot.sh"); return 1
    bad = 0
    for name in names:
        c = head(os.path.join(base, f"dom_v3c3_{name}.html")); v = head(os.path.join(base, f"dom_v4c3_{name}.html")); w = head(os.path.join(base, f"dom_v4c3_pow_{name}.html"))
        print(f"== {name} ==")
        if c is None or v is None or w is None: print("  読めない（Plotly か DOM の取得に失敗）"); bad += 1; continue
        ok1 = (v[0] == c[0]); ok2 = (v[1] == c[1])
        print(f"  OFF：h1 {'一致' if ok1 else '不一致'}　h2 {'一致' if ok2 else '不一致'}")
        if not ok1: print(f"    原本 h1: {c[0][:200]}\n    v4c3 h1: {v[0][:200]}")
        bad += (not ok1) + (not ok2)
        bad += check_pow(w, c[0], c[1], "ON ：")
    print("== auto6s（Auto を仮想時間 6 秒）==")
    a = head(os.path.join(base, "dom_v4c3_auto6s.html")); b = head(os.path.join(base, "dom_v4c3_pow_auto6s.html"))
    if a is None or b is None: print("  読めない"); bad += 1
    else:
        t = re.search(r"ティック (\d+)", b[1]); print(f"  ティック {t.group(1) if t else '?'}")
        bad += check_pow(b, a[0], a[1], "ON ：")
    print("対照 " + ("一致" if bad == 0 else f"不一致 {bad}"))
    return 0 if bad == 0 else 2

if __name__ == "__main__":
    sys.exit(main())
