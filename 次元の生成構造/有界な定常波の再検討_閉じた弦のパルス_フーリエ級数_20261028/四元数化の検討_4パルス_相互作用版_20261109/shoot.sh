#!/bin/sh
# 確認と撮影（ヘッドレス Chrome、180 秒で打ち切り）。算出根拠 §9 の K0〜K6。
#  K0：κ = 0 の本版と、原本 v11q4（同じ配置：既定に dl2=2, dl4=1 を URL で）の DOM を tick 40 で比較。
#  K1〜K6：?test=1 の自己検査（<pre id="ktest">）。
#  走行：既定（κ = 1）で tick 70, 85, 90, 95, 110 の DOM と、t0 / tick 90 / tick 130 の静止画。
set -e
cd "$(dirname "$0")"
CH="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
D="file://$(pwd)"
shot() { perl -e 'alarm shift; exec @ARGV' 180 "$CH" --headless=new --hide-scrollbars --window-size=1800,1400 --virtual-time-budget=$3 --screenshot="figs/$1.png" "$2" >/dev/null 2>&1 || echo "shot $1 timed out"; echo "shot $1"; }
dom()  { perl -e 'alarm shift; exec @ARGV' 180 "$CH" --headless=new --virtual-time-budget=$3 --dump-dom "$2" 2>/dev/null > "results/dom_$1.html" || true; echo "dom  $1 ($(wc -c < "results/dom_$1.html") bytes)"; }
dom ktest "$D/index_v12q4i.html?test=1&auto=0" 3000
dom v11q4_k0_tick40  "$D/原本/index_v11q4.html?dl2=2&dl4=1" 6000
dom v12q4i_k0_tick40 "$D/index_v12q4i.html?kap=0" 6000
for T in 10500 12750 13500 14250 16500; do dom "v12q4i_default_vt$T" "$D/index_v12q4i.html" $T; done
shot v12q4i_default_t0     "$D/index_v12q4i.html?auto=0" 6000
shot v12q4i_default_tick90 "$D/index_v12q4i.html" 13500
shot v12q4i_default_tick130 "$D/index_v12q4i.html" 19500
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
        t = re.sub(r'<pre id="ktest"[^>]*></pre>\s*', '', t); t = re.sub(r'φ₀ = -?[\d.]+°, ', '', t); t = re.sub(r'　　κ = [^　]*（まだ交換なし）', '', t)   # desc　　κ…　　和 → desc　　和
        t = t.replace('<title>閉じた弦の四元数パルス 4 波 相互作用版（index_v11q4.html に重なりの核と二チャネル交換）</title>', '<title>閉じた弦の四元数パルス 4 波（index_v8q3.html の 4 パルス版）</title>')
    return t
A, B = strip(a), strip(b, True); res["K0_identical"] = (A == B)
if A != B:
    import difflib; d = [l for l in difflib.unified_diff(A.splitlines(), B.splitlines(), lineterm="", n=0) if l[:1] in "+-" and l[:3] not in ("+++", "---")]; res["K0_diff_head"] = [x[:240] for x in d[:6]]
runs = []
for T in (10500, 12750, 13500, 14250, 16500):
    s = open(f"results/dom_v12q4i_default_vt{T}.html", encoding="utf-8").read()
    k = re.search(r'ティック (\d+)', s); th = re.findall(r'(\d): θ = (-?[\d.]+)°, φ = 1°, v = (-?\d+)°/t, δ = (-?\d+)°/t, φ₀ = (-?[\d.]+)°', s)
    nn = re.findall(r'A = ([\d.]+) n = \(([-\d., ]+)\) R\d = ([\d.]+)', s); F = re.search(r'核 F: ([^　]*)　Δ_F: ([^　]*)　ΣA² = ([\d.]+)', s)
    runs.append({"vt": T, "tick": k.group(1) if k else None, "theta_v_delta_phi0": th, "A_n_R": nn, "F": F.group(1) if F else None, "dF": F.group(2) if F else None, "sumA2": F.group(3) if F else None})
res["default_runs"] = runs
json.dump(res, open("results/checks_v12q4i.json", "w"), ensure_ascii=False, indent=2)
print(json.dumps(res["ktest"], ensure_ascii=False, indent=1)); print("K0_identical:", res["K0_identical"], res.get("K0_diff_head", ""))
for r in runs: print("tick", r["tick"], "| φ₀:", [x[4] for x in r["theta_v_delta_phi0"]], "| A:", [x[0] for x in r["A_n_R"]], "| n:", [x[1] for x in r["A_n_R"]], "| F:", r["F"], "| ΔF:", r["dF"], "| ΣA²:", r["sumA2"])
PY
echo done
