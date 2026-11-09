#!/bin/sh
# 確認と撮影（ヘッドレス Chrome、180 秒で打ち切り）。算出根拠 v0.6 §9 の K0〜K7。
#  K0：相互作用を外した本版（?kap=0）と原本 v11q4（同じ初期値を URL で）を tick 40 の DOM で比較。
#  K1〜K7：?test=1 の自己検査（<pre id="ktest">）。K1 保存則、K3 一歩の角と門、K4 δ が同じなら恒等、K5 通過と閉包（n = 20、δ₂ = 18 / 9 / 2、同質量と 1:4）、K6 置換、K7 δ₂ の走査。
#  走行：既定（δ = 0, 2, 0, 2）と ?dl2=9&dl4=9、?dl2=18&dl4=18 を tick 70, 90, 100, 130 の DOM。静止画は tick 130。
set -e
cd "$(dirname "$0")"
CH="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
D="file://$(pwd)"
shot() { perl -e 'alarm shift; exec @ARGV' 180 "$CH" --headless=new --hide-scrollbars --window-size=1800,1400 --virtual-time-budget=$3 --screenshot="figs/$1.png" "$2" >/dev/null 2>&1 || echo "shot $1 timed out"; echo "shot $1"; }
dom()  { perl -e 'alarm shift; exec @ARGV' 180 "$CH" --headless=new --virtual-time-budget=$3 --dump-dom "$2" 2>/dev/null > "results/dom_$1.html" || true; echo "dom  $1 ($(wc -c < "results/dom_$1.html") bytes)"; }
dom ktest "$D/index_v12q4i.html?test=1&auto=0" 6000
DEF='th1=20&phi1=20&v1=1&dl1=0&A1=0.5&ni1=1&nj1=0&nk1=0&th2=200&phi2=20&v2=-1&dl2=2&A2=1&ni2=1&nj2=0&nk2=0.05&th3=220&phi3=20&v3=1&dl3=0&A3=1.5&ni3=0&nj3=1&nk3=0&th4=340&phi4=20&v4=-1&dl4=2&A4=2&ni4=0&nj4=1&nk4=0.05'
dom v11q4_default_tick40  "$D/原本/index_v11q4.html?$DEF" 6000
dom v12q4i_k0_tick40      "$D/index_v12q4i.html?kap=0" 6000
for V in "default:" "dl9:dl2=9&dl4=9" "dl18:dl2=18&dl4=18"; do NAME="${V%%:*}"; Q="${V#*:}"; for T in 10500 13500 15000 19500; do dom "v12q4i_${NAME}_vt$T" "$D/index_v12q4i.html?$Q" $T; done; done
shot v12q4i_t0            "$D/index_v12q4i.html?auto=0" 3000
shot v12q4i_default_tick130 "$D/index_v12q4i.html" 19500
shot v12q4i_dl9_tick130   "$D/index_v12q4i.html?dl2=9&dl4=9" 19500
shot v12q4i_dl18_tick130  "$D/index_v12q4i.html?dl2=18&dl4=18" 19500
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
    if new:   # 本版だけにある要素と文字：相互作用のチェック、検査欄、φ₀ の表示、相互作用の行、題
        t = re.sub(r'<label class="ck"><input type="checkbox" id="kap"[^>]*>相互作用（[^<]*）</label>\s*', '', t)
        t = re.sub(r'<pre id="ktest"[^>]*></pre>\s*', '', t); t = re.sub(r'φ₀ = -?[\d.]+°, ', '', t); t = re.sub(r'　　相互作用 なし（まだ交換なし）', '', t)
        t = t.replace('<title>閉じた弦の四元数パルス 4 波 相互作用版（index_v11q4.html に、重なっている tick の有限位数回転）</title>', '<title>閉じた弦の四元数パルス 4 波（index_v8q3.html の 4 パルス版）</title>')
    return t
A, B = strip(a), strip(b, True); res["K0_identical"] = (A == B)
if A != B:
    la, lb = A.splitlines(), B.splitlines(); res["K0_lines"] = [len(la), len(lb)]
    for i, (x, y) in enumerate(zip(la, lb)):
        if x != y:
            j = next((k for k in range(min(len(x), len(y))) if x[k] != y[k]), 0); res["K0_first_diff"] = {"line": i, "A": x[max(0, j-80):j+120], "B": y[max(0, j-80):j+120]}; break
runs = {}
for NAME in ("default", "dl9", "dl18"):
    runs[NAME] = []
    for T in (10500, 13500, 15000, 19500):
        s = open(f"results/dom_v12q4i_{NAME}_vt{T}.html", encoding="utf-8").read()
        k = re.search(r'ティック (\d+)', s); th = re.findall(r'(\d): θ = (-?[\d.]+)°, φ = 20°, v = (-?[\d.]+)°/t, δ = (-?\d+)°/t, φ₀ = (-?[\d.]+)°', s)
        nn = re.findall(r'A = ([\d.]+) n = \(([-\d., ]+)\) R\d = ([\d.]+)', s); F = re.search(r'相互作用 あり　d: ([^　]*)　α: ([^　]*)　核 F（表示）: ([^　]*)　Σ m v n = \(([^)]*)\)　Σ ½ m v² = ([\d.]+)', s)
        runs[NAME].append({"vt": T, "tick": k.group(1) if k else None, "theta": [x[1] for x in th], "v": [x[2] for x in th], "delta": [x[3] for x in th], "A": [x[0] for x in nn], "n": [x[1] for x in nn], "d": F.group(1) if F else None, "alpha": F.group(2) if F else None, "F": F.group(3) if F else None, "P": F.group(4) if F else None, "K": F.group(5) if F else None})
res["runs"] = runs
json.dump(res, open("results/checks_v12q4i.json", "w"), ensure_ascii=False, indent=2)
kt = res["ktest"]; print("K0 identical:", res["K0_identical"], res.get("K0_first_diff", ""))
for k in sorted(kt):
    if k.startswith("K7_scan") or k.startswith("K5"): continue
    print(" ", k, "=", kt[k])
for k in [k for k in kt if k.startswith("K5")]: print(" ", k, {kk: vv for kk, vv in kt[k].items() if kk != "trace"})
for NAME in ("default", "dl9", "dl18"):
    print("==", NAME)
    for r in runs[NAME]: print(" tick", r["tick"], "| v:", r["v"], "| n:", r["n"], "\n      d:", r["d"], "| α:", r["alpha"], "| P:", r["P"], "| K:", r["K"])
PY
echo done
