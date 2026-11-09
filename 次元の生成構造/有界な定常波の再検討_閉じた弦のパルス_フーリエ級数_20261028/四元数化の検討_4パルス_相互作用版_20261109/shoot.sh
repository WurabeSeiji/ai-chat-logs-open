#!/bin/sh
# 確認と撮影（ヘッドレス Chrome、180 秒で打ち切り）。算出根拠 §9 の K0〜K7。
#  K0：κ = 0・包絡 sinc の本版と、原本 v11q4（同じ配置：既定に dl2=2, dl4=1 を URL で）の DOM を tick 40 で比較。
#  K1〜K7：?test=1 の自己検査（<pre id="ktest">）。K1〜K6 は全量、K7 は率の通過。
#  走行：既定（κ = 1、包絡 cos²）で交換 率／全量 それぞれ tick 70, 90, 110, 130 の DOM と静止画。
set -e
cd "$(dirname "$0")"
CH="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
D="file://$(pwd)"
shot() { perl -e 'alarm shift; exec @ARGV' 180 "$CH" --headless=new --hide-scrollbars --window-size=1800,1400 --virtual-time-budget=$3 --screenshot="figs/$1.png" "$2" >/dev/null 2>&1 || echo "shot $1 timed out"; echo "shot $1"; }
dom()  { perl -e 'alarm shift; exec @ARGV' 180 "$CH" --headless=new --virtual-time-budget=$3 --dump-dom "$2" 2>/dev/null > "results/dom_$1.html" || true; echo "dom  $1 ($(wc -c < "results/dom_$1.html") bytes)"; }
rm -f results/dom_v12q4i_default_vt*.html figs/v12q4i_default_*.png
dom ktest "$D/index_v12q4i.html?test=1&auto=0" 4000
dom v11q4_k0_tick40  "$D/原本/index_v11q4.html?dl2=2&dl4=1" 6000
dom v12q4i_k0_tick40 "$D/index_v12q4i.html?kap=0&env=sinc" 6000
for XM in rate full; do for T in 10500 13500 16500 19500; do dom "v12q4i_${XM}_vt$T" "$D/index_v12q4i.html?xm=$XM" $T; done; done
shot v12q4i_rate_t0      "$D/index_v12q4i.html?auto=0" 6000
shot v12q4i_rate_tick90  "$D/index_v12q4i.html" 13500
shot v12q4i_rate_tick130 "$D/index_v12q4i.html" 19500
shot v12q4i_full_tick130 "$D/index_v12q4i.html?xm=full" 19500
python3 - <<'PY'
# -*- coding: utf-8 -*-
import re, json
res = {}
s = open("results/dom_ktest.html", encoding="utf-8").read(); m = re.search(r'KTEST (\{.*?\})</pre>', s, re.S); res["ktest"] = json.loads(m.group(1)) if m else None
a = open("results/dom_v11q4_k0_tick40.html", encoding="utf-8").read(); b = open("results/dom_v12q4i_k0_tick40.html", encoding="utf-8").read()
def strip(t, new=False):
    t = re.sub(r'<script>.*?</script>', '', t, flags=re.S)
    t = re.sub(r'(<input [^>]*?)value="[^"]*"', r'\1value=""', t); t = re.sub(r'(<input [^>]*?) checked(="")?', r'\1', t)
    m = re.search(r'plotly\.js-style-modebar-([0-9a-f]+)', t); t = t.replace(m.group(1), 'X') if m else t
    if new:
        t = re.sub(r'<div class="sl"><span>κ 相互作用</span>.*?</div>\s*', '', t, flags=re.S)
        t = re.sub(r'<span class="grp">包絡 .*?</span>\s*<span class="grp">交換 .*?</span>\s*', '', t, flags=re.S)
        t = re.sub(r'<pre id="ktest"[^>]*></pre>\s*', '', t); t = re.sub(r'φ₀ = -?[\d.]+°, ', '', t); t = re.sub(r'　　κ = [^　]*（まだ交換なし）', '', t)
        t = t.replace('<title>閉じた弦の四元数パルス 4 波 相互作用版（index_v11q4.html に重なりの核と二チャネル交換）</title>', '<title>閉じた弦の四元数パルス 4 波（index_v8q3.html の 4 パルス版）</title>')
    return t
A, B = strip(a), strip(b, True); res["K0_identical"] = (A == B)
if A != B:
    la, lb = A.splitlines(), B.splitlines(); res["K0_lines"] = [len(la), len(lb)]
    for i, (x, y) in enumerate(zip(la, lb)):
        if x != y:
            j = next((k for k in range(min(len(x), len(y))) if x[k] != y[k]), 0); res["K0_first_diff"] = {"line": i, "A": x[max(0, j-80):j+120], "B": y[max(0, j-80):j+120]}; break
runs = {}
for XM in ("rate", "full"):
    runs[XM] = []
    for T in (10500, 13500, 16500, 19500):
        s = open(f"results/dom_v12q4i_{XM}_vt{T}.html", encoding="utf-8").read()
        k = re.search(r'ティック (\d+)', s); th = re.findall(r'(\d): θ = (-?[\d.]+)°, φ = 1°, v = (-?\d+)°/t, δ = (-?\d+)°/t, φ₀ = (-?[\d.]+)°', s)
        nn = re.findall(r'A = ([\d.]+) n = \(([-\d., ]+)\) R\d = ([\d.]+)', s); F = re.search(r'核 F: ([^　]*)　Δ_F: ([^　]*)　ΣA² = ([\d.]+)', s)
        runs[XM].append({"vt": T, "tick": k.group(1) if k else None, "phi0": [x[4] for x in th], "A": [x[0] for x in nn], "n": [x[1] for x in nn], "F": F.group(1) if F else None, "dF": F.group(2) if F else None, "sumA2": F.group(3) if F else None})
res["default_runs"] = runs
json.dump(res, open("results/checks_v12q4i.json", "w"), ensure_ascii=False, indent=2)
kt = res["ktest"]; print("K0 identical:", res["K0_identical"], res.get("K0_first_diff", ""))
for k in sorted(kt): print(" ", k, "=", kt[k])
for XM in ("rate", "full"):
    print("==", XM)
    for r in runs[XM]: print(" tick", r["tick"], "| φ₀:", r["phi0"], "| A:", r["A"], "| n:", r["n"], "| F:", r["F"], "| ΣA²:", r["sumA2"])
PY
echo done
