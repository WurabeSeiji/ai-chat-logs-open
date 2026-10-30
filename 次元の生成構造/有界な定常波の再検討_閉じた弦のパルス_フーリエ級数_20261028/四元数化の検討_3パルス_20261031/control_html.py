#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""HTML の対照：shoot.sh が取った DOM（results/dom_c3_*.html ＝ 原本 index_v2c3.html、results/dom_q3_*.html ＝ index_v2q3.html、同じ URL 引数、n は既定の k）
から右上の読める量（#h1・#h2 の文字列）を取り出し、共通の項目の値が一致することを確かめる（2026-10-31）。
名前が変わった項目：max Im / min Im → k 成分 [min, max]、Σc²/R² の虚部 → (i, j, k) の k 成分、距離の比 Δ₁₃/Δ₂₃ → x₁₃/x₂₃、Re Z → Re Q。
四元数版にだけある項目（n、i 成分、j 成分、n·n、n×n、Δ⃗、V⃗）は対照の外（一覧だけ出す）。"""
import re, os, sys, html

base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
NUM = r"[-+−]?\d+(?:\.\d+)?"

def head(path):
    s = open(path, encoding="utf-8").read()
    h1 = re.search(r'id="h1">(.*?)</div>', s, re.S); h2 = re.search(r'id="h2"[^>]*>(.*?)</div>', s, re.S)
    if not h1 or not h2: return None
    return html.unescape(h1.group(1)), html.unescape(h2.group(1))

def nums(s):   # 文字列の中の数（全角マイナスも）
    return [float(x.replace("−", "-")) for x in re.findall(NUM, s)]

def items(h1, h2):
    """項目名 → 数のリスト。全角空白・｜で区切り、「名前 = 値」「名前：値」を拾う。対ごとの読みは 12：/23：/13：の前置きを付ける"""
    out = {}
    for seg in re.split(r"[　｜]", h1 + "　" + h2):
        seg = seg.strip()
        if not seg: continue
        m = re.match(r"^(12|23|13)：(.*)$", seg)
        if m:   # 対ごとの読み（Δ = …°、Φ = …°、cos Θ = …（…°）、tan α = …、…）
            for part in m.group(2).split("、"):
                k = re.split(r"\s*=\s*", part, 1)
                if len(k) == 2: out[f"{m.group(1)} {k[0].strip()}"] = nums(k[1])
            continue
        if seg.startswith("観測者") or seg.startswith("単位"): out["単位"] = nums(seg); continue
        m2 = re.match(r"^([ijk]) 成分 \[(.*)\]$", seg)   # 四元数版の「k 成分 [min, max]」
        if m2: out[f"{m2.group(1)} 成分"] = nums(m2.group(2)); continue
        k = re.split(r"\s*=\s*|：", seg, 1)
        if len(k) == 2 and nums(k[1]): out[k[0].strip()] = nums(k[1])
        elif nums(seg): out[seg[:12]] = nums(seg)
    return out

def main():
    names = sorted(set(f[7:-5] for f in os.listdir(base) if f.startswith("dom_c3_") and f.endswith(".html")))
    if not names: print("DOM がない。先に sh shoot.sh"); return 1
    bad = 0
    for name in names:
        c = head(os.path.join(base, f"dom_c3_{name}.html")); q = head(os.path.join(base, f"dom_q3_{name}.html"))
        print(f"== {name} ==")
        if c is None or q is None: print("  読めない（Plotly か DOM の取得に失敗）"); bad += 1; continue
        ic, iq = items(*c), items(*q)
        # 名前の対応（複素 → 四元数）
        alias = {"距離の比 Δ₁₃/Δ₂₃": "距離の比 x₁₃/x₂₃"}
        shared = 0; mism = []
        for k, v in ic.items():
            kq = alias.get(k, k)
            if re.match(r"^\d: θ", k) and kq in iq:   # 各パルスの行：θ, φ, δ, M, 奇, 偶 と最後の R_i を比べる（四元数版は間に n の 3 成分が入る）
                w = iq[kq]; ok = v[:6] == w[:6] and abs(v[-1] - w[-1]) < 5e-4; shared += 1
                if not ok: mism.append((k + "（θ, φ, δ, M, 奇, 偶, R）", v, w))
            elif kq in iq:
                w = iq[kq]
                if k.startswith("Σc²/R²"):   # 複素 a ± bi と 四元数 s + (x, y, z)
                    ok = len(w) >= 4 and abs(v[0] - w[0]) < 5e-4 and abs(abs(v[1]) - abs(w[3])) < 5e-4 and abs(w[1]) < 5e-4 and abs(w[2]) < 5e-4
                    if v[1] != 0 and len(w) >= 4:   # 符号：複素の表示は「+ b」「− b」で b は絶対値なので文字列で見る
                        sgn_c = "−" if ("− " + f"{abs(v[1]):.4f}") in c[1] else "+"
                        ok = ok and ((w[3] < 0) == (sgn_c == "−"))
                else:
                    ok = len(v) == len(w) and all(abs(a - b) < 5e-4 for a, b in zip(v, w))
                shared += 1
                if not ok: mism.append((k, v, w))
            elif k == "max Im":   # 複素「max Im = a, min Im = b」→ 四元数「k 成分 [b, a]」（n = k なので虚部は k 成分）
                kk = iq.get("k 成分", []); ok = len(kk) == 2 and len(v) == 2 and abs(v[0] - kk[1]) < 5e-4 and abs(v[1] - kk[0]) < 5e-4; shared += 1
                if not ok: mism.append(("max Im / min Im → k 成分", v, kk))
        only_q = [k for k in iq if k not in ic and k not in alias.values()]
        print(f"  共通の項目 {shared} 個：{'すべて一致' if not mism else str(len(mism)) + ' 個が不一致'}")
        for k, v, w in mism: print(f"    不一致 {k}: 複素 {v} / 四元数 {w}")
        print(f"  四元数版にだけある項目：{', '.join(only_q)}")
        print(f"  複素 h1: {c[0][:160]}")
        print(f"  四元数 h1: {q[0][:200]}")
        bad += len(mism)
    print("対照 " + ("一致" if bad == 0 else f"不一致 {bad}"))
    return 0 if bad == 0 else 2

if __name__ == "__main__":
    sys.exit(main())
