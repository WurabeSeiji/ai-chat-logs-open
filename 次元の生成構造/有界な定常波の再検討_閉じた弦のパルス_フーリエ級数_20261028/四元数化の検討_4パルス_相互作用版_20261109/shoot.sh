#!/bin/sh
# 確認と撮影（ヘッドレス Chrome、180 秒で打ち切り）。算出根拠 v0.5 §9 の K0〜K7。
#  K0：κ = 0 の本版と原本 v11q4 を、どちらも URL なしの既定で tick 40 の DOM 比較（既定値は同一）。
#  K1〜K7：?test=1 の自己検査（<pre id="ktest">）。K1 運動量・運動エネルギーの保存、K3 正面 Δ_F=180 の反転、K4 恒等、K5 質量 A² の違う正面＝古典の弾性衝突、K6 置換、K7 率／全量の通過を相対速度 2 と 1 で（tick ごとの α と散乱面の法線の記録つき）。
#  走行：既定（X 軸の対・Y 軸の対）を率／全量で tick 50, 60, 70, 90, 100, 130 の DOM。静止画は率と全量。
set -e
cd "$(dirname "$0")"
CH="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
D="file://$(pwd)"
shot() { perl -e 'alarm shift; exec @ARGV' 180 "$CH" --headless=new --hide-scrollbars --window-size=1800,1400 --virtual-time-budget=$3 --screenshot="figs/$1.png" "$2" >/dev/null 2>&1 || echo "shot $1 timed out"; echo "shot $1"; }
dom()  { perl -e 'alarm shift; exec @ARGV' 180 "$CH" --headless=new --virtual-time-budget=$3 --dump-dom "$2" 2>/dev/null > "results/dom_$1.html" || true; echo "dom  $1 ($(wc -c < "results/dom_$1.html") bytes)"; }
dom ktest "$D/index_v12q4i.html?test=1&auto=0" 4000
# K0：本版の既定（斜めの衝突・振幅ばらつき）を原本には URL で与え、κ = 0 の本版の既定と比べる
DEF='th1=20&phi1=20&v1=1&dl1=0&A1=0.5&ni1=1&nj1=0&nk1=0&th2=200&phi2=20&v2=-1&dl2=2&A2=1&ni2=1&nj2=0&nk2=0.05&th3=220&phi3=20&v3=1&dl3=0&A3=1.5&ni3=0&nj3=1&nk3=0&th4=340&phi4=20&v4=-1&dl4=2&A4=2&ni4=0&nj4=1&nk4=0.05'
dom v11q4_default_tick40  "$D/原本/index_v11q4.html?$DEF" 6000
dom v12q4i_k0_tick40      "$D/index_v12q4i.html?kap=0" 6000
# 既定（URL なし）の走行：対 3–4 が tick 60、対 1–2 が tick 90 に出会う。tick 50, 60, 70, 90, 100, 130
for XM in rate full; do for T in 7500 9000 10500 13500 15000 19500; do dom "v12q4i_${XM}_vt$T" "$D/index_v12q4i.html?xm=$XM" $T; done; done
shot v12q4i_rate_t0      "$D/index_v12q4i.html?auto=0" 3000
shot v12q4i_rate_tick60  "$D/index_v12q4i.html" 9000
shot v12q4i_rate_tick90  "$D/index_v12q4i.html" 13500
shot v12q4i_rate_tick130 "$D/index_v12q4i.html" 19500
shot v12q4i_full_tick130 "$D/index_v12q4i.html?xm=full" 19500
python3 - <<'PY'
# -*- coding: utf-8 -*-
import re, json
res = {}
s = open("results/dom_ktest.html", encoding="utf-8").read(); m = re.search(r'KTEST (\{.*?\})</pre>', s, re.S); res["ktest"] = json.loads(m.group(1)) if m else None
a = open("results/dom_v11q4_default_tick40.html", encoding="utf-8").read(); b = open("results/dom_v12q4i_k0_tick40.html", encoding="utf-8").read()
def strip(t, new=False):
    t = re.sub(r'<script>.*?</script>', '', t, flags=re.S)
    t = re.sub(r'(<input [^>]*?)value="[^"]*"', r'\1value=""', t); t = re.sub(r'(<input [^>]*?) checked(="")?', r'\1', t)
    m = re.search(r'plotly\.js-style-modebar-([0-9a-f]+)', t); t = t.replace(m.group(1), 'X') if m else t
    if new:   # 本版だけにある要素と文字：κ スライダー、交換のラジオ、検査欄、φ₀ の表示、κ の行、題
        t = re.sub(r'<div class="sl"><span>κ 相互作用</span>.*?</div>\s*', '', t, flags=re.S)
        t = re.sub(r'<span class="grp">交換 .*?</span>\s*', '', t, flags=re.S)
        t = re.sub(r'<pre id="ktest"[^>]*></pre>\s*', '', t); t = re.sub(r'φ₀ = -?[\d.]+°, ', '', t); t = re.sub(r'　　κ = [^　]*（まだ交換なし）', '', t)
        t = t.replace('<title>閉じた弦の四元数パルス 4 波 相互作用版（index_v11q4.html に重なりの核とゲージ回転）</title>', '<title>閉じた弦の四元数パルス 4 波（index_v8q3.html の 4 パルス版）</title>')
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
    for T in (7500, 9000, 10500, 13500, 15000, 19500):
        s = open(f"results/dom_v12q4i_{XM}_vt{T}.html", encoding="utf-8").read()
        k = re.search(r'ティック (\d+)', s); th = re.findall(r'(\d): θ = (-?[\d.]+)°, φ = 20°, v = (-?[\d.]+)°/t, δ = (-?\d+)°/t, φ₀ = (-?[\d.]+)°', s)
        nn = re.findall(r'A = ([\d.]+) n = \(([-\d., ]+)\) R\d = ([\d.]+)', s); F = re.search(r'核 F: ([^　]*)　Δ_F: ([^　]*)　α: ([^　]*)　Σ m v n = \(([^)]*)\)　Σ ½ m v² = ([\d.]+)', s)
        runs[XM].append({"vt": T, "tick": k.group(1) if k else None, "theta": [x[1] for x in th], "v": [x[2] for x in th], "phi0": [x[4] for x in th], "A": [x[0] for x in nn], "n": [x[1] for x in nn], "F": F.group(1) if F else None, "dF": F.group(2) if F else None, "alpha": F.group(3) if F else None, "P": F.group(4) if F else None, "K": F.group(5) if F else None})
res["default_runs_xy"] = runs
json.dump(res, open("results/checks_v12q4i.json", "w"), ensure_ascii=False, indent=2)
kt = res["ktest"]; print("K0 identical:", res["K0_identical"], res.get("K0_first_diff", ""))
for k in sorted(kt): print(" ", k, "=", kt[k])
for XM in ("rate", "full"):
    print("==", XM)
    for r in runs[XM]: print(" tick", r["tick"], "| θ:", r["theta"], "| v:", r["v"], "| n:", r["n"], "\n      F:", r["F"], "| ΔF:", r["dF"], "| α:", r["alpha"], "| P:", r["P"], "| K:", r["K"])
PY
echo done
