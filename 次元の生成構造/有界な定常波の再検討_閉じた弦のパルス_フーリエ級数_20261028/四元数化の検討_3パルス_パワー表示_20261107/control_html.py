#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""HTML の対照（2026-11-07）：shoot.sh が取った DOM から右上の読める量（#h1・#h2）を取り出して比べる。
(1) OFF：index_v4q3.html は原本 index_v3q3.html と h1・h2 が完全に一致する（同じ URL 引数）。
(2) ON（pow=1）：h2 は原本と一致（読める量は表示モードで変わらない）。h1 は「　パワー表示」より前が原本の「　枠」より前と一致。
    軸ごとの ψ 平均の和が合計 ⟨|Q|²⟩_ψ に一致し（丸め 2e-4）、合計は係数だけの閉じた式（check_power_mean_v4q3.closed_mean）に一致する（5e-4）。"""
import re, os, sys, html
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_power_mean_v4q3 as chk

base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
NUM = r"(\d+\.\d+)"
POW = (r"　パワー表示（軸ごと、max／ψ 平均）：W² " + NUM + "／" + NUM + "　i² " + NUM + "／" + NUM + "　j² " + NUM + "／" + NUM + "　k² " + NUM + "／" + NUM + r"　合計 ⟨\|Q\|²⟩_ψ = " + NUM)

def head(path):
    s = open(path, encoding="utf-8").read()
    h1 = re.search(r'id="h1">(.*?)</div>', s, re.S); h2 = re.search(r'id="h2"[^>]*>(.*?)</div>', s, re.S)
    if not h1 or not h2: return None
    return html.unescape(h1.group(1)), html.unescape(h2.group(1))

def check_pow(v, ref, tag, name):
    bad = 0
    m = re.search(POW, v[0])
    if not m: print(f"  {tag}パワー表示の項目がない"); return 1
    g = list(map(float, m.groups())); means = [g[1], g[3], g[5], g[7]]; total = g[8]
    ok_sum = abs(sum(means) - total) < 2e-4
    print(f"  {tag}max／ψ平均  W² {g[0]:.4f}／{g[1]:.4f}  i² {g[2]:.4f}／{g[3]:.4f}  j² {g[4]:.4f}／{g[5]:.4f}  k² {g[6]:.4f}／{g[7]:.4f}  合計 {total:.4f}（軸の和 {sum(means):.4f} {'一致' if ok_sum else '不一致'}）")
    bad += (not ok_sum)
    if name in chk.CASES:
        cm = chk.closed_mean(name); ok_c = abs(total - cm) < 5e-4
        print(f"  {tag}合計の ψ 平均 = 係数からの閉じた式 {cm:.4f}：{'一致' if ok_c else '不一致'}"); bad += (not ok_c)
    if ref is not None:
        pre_v = v[0][:m.start()]; pre_r = ref[0][:ref[0].index("　枠 ±")]
        ok1 = (pre_v == pre_r); ok2 = (v[1] == ref[1])
        print(f"  {tag}h1（パワー表示より前）: {'一致' if ok1 else '不一致'}　h2: {'一致' if ok2 else '不一致'}")
        bad += (not ok1) + (not ok2)
    return bad

def main():
    names = sorted(set(f[9:-5] for f in os.listdir(base) if f.startswith("dom_v3q3_") and f.endswith(".html")))
    if not names: print("DOM がない。先に sh shoot.sh"); return 1
    bad = 0
    for name in names:
        c = head(os.path.join(base, f"dom_v3q3_{name}.html")); v = head(os.path.join(base, f"dom_v4q3_{name}.html")); w = head(os.path.join(base, f"dom_v4q3_pow_{name}.html"))
        print(f"== {name} ==")
        if c is None or v is None or w is None: print("  読めない（Plotly か DOM の取得に失敗）"); bad += 1; continue
        ok1 = (v[0] == c[0]); ok2 = (v[1] == c[1])
        print(f"  OFF：h1 {'一致' if ok1 else '不一致'}　h2 {'一致' if ok2 else '不一致'}"); bad += (not ok1) + (not ok2)
        bad += check_pow(w, c, "ON ：", name)
    print("== auto6s（Auto を仮想時間 6 秒）==")
    a = head(os.path.join(base, "dom_v4q3_auto6s.html")); b = head(os.path.join(base, "dom_v4q3_pow_auto6s.html"))
    if a is None or b is None: print("  読めない"); bad += 1
    else:
        t = re.search(r"ティック (\d+)", b[1]); print(f"  ティック {t.group(1) if t else '?'}"); bad += check_pow(b, a, "ON ：", "auto6s")
    print("対照 " + ("一致" if bad == 0 else f"不一致 {bad}"))
    return 0 if bad == 0 else 2

if __name__ == "__main__":
    sys.exit(main())
