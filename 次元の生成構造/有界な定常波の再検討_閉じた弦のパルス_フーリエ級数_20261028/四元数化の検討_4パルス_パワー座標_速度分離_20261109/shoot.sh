#!/bin/sh
# 対照と撮影（ヘッドレス Chrome、180 秒で打ち切り）。
# 対照：既定（v = φ = 1）で原本 v10q4 と本版の DOM を比べる。本版だけにある v スライダーの要素・「v = 1°/t」の文字・<script> 本文・plotly の乱数 id を除けば一致するはず。
# 撮影：正面衝突の例（対 1–2 は軸 x で θ₁=20,v₁=+1 / θ₂=200,v₂=−1 → tick 90、対 3–4 は軸 y で θ₃=40,v₃=+2 / θ₄=280,v₄=−1 → tick 80）
set -e
cd "$(dirname "$0")"
CH="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
D="file://$(pwd)"
shot() { perl -e 'alarm shift; exec @ARGV' 180 "$CH" --headless=new --hide-scrollbars --window-size=1800,1400 --virtual-time-budget=$3 --screenshot="figs/$1.png" "$2" >/dev/null 2>&1 || echo "shot $1 timed out"; echo "shot $1"; }
dom()  { perl -e 'alarm shift; exec @ARGV' 180 "$CH" --headless=new --virtual-time-budget=$3 --dump-dom "$2" 2>/dev/null > "results/dom_$1.html" || true; echo "dom  $1 ($(wc -c < "results/dom_$1.html") bytes)"; }
# 対照の配置は v10q4 の既定（θ=97,194,58,300、φ=1、M=100、δ=20,0,10,−10、n=x,y,z,(1,1,1)）。v11q4 の既定は正面衝突に変えたので URL で明示し、v は 1 を渡す（v10q4 は v を無視）
OLD='th1=97&phi1=1&v1=1&M1=100&dl1=20&ni1=1&nj1=0&nk1=0&th2=194&phi2=1&v2=1&M2=100&dl2=0&ni2=0&nj2=1&nk2=0&th3=58&phi3=1&v3=1&M3=100&dl3=10&ni3=0&nj3=0&nk3=1&th4=300&phi4=1&v4=1&M4=100&dl4=-10&ni4=1&nj4=1&nk4=1&rc=1&auto=0'
dom v10q4_default_rc1 "$D/原本/index_v10q4.html?$OLD" 6000
dom v11q4_default_rc1 "$D/index_v11q4.html?$OLD" 6000
python3 - <<'PY'
import re
a = open("results/dom_v10q4_default_rc1.html", encoding="utf-8").read(); b = open("results/dom_v11q4_default_rc1.html", encoding="utf-8").read()
b = re.sub(r'<div class="sl"><span>v[₁₂₃₄]</span><input id="vv\d"[^>]*><span class="v" id="vv\dv">1°/t</span></div>\s*', '', b)
b = b.replace("v = 1°/t, ", "")
# 比べるのは描画結果の文字と要素。<script> の本文（意図した差）と plotly の乱数 id は除く
def strip(t):
    t = re.sub(r'<script>.*?</script>', '', t, flags=re.S)
    t = re.sub(r'(<input [^>]*?)value="[^"]*"', r'\1value=""', t)   # input の value 属性は HTML の初期値で、URL で与えた生の状態ではない（表示の数値は別の span で比べる）
    t = re.sub(r'(<input [^>]*?) checked(="")?', r'\1', t)   # radio の checked 属性も HTML の既定（本版は座標が既定）。描画は rcoord で決まり、それは URL の rc で揃えてある
    m = re.search(r'plotly\.js-style-modebar-([0-9a-f]+)', t)   # plotly が描画ごとに振る乱数 id（modebar と defs に同じ値）
    return t.replace(m.group(1), 'X') if m else t
a, b = strip(a), strip(b)
same = (a == b); print("対照（既定、v スライダーを除く）:", "一致" if same else "不一致")
if not same:
    import difflib; d = [l for l in difflib.unified_diff(a.splitlines(), b.splitlines(), lineterm="", n=0) if l[:1] in "+-" and l[:3] not in ("+++", "---")]; print("\n".join(x[:160] for x in d[:12]))
open("results/control_default.txt", "w", encoding="utf-8").write("一致\n" if same else "不一致\n")
PY
HEAD='th1=20&phi1=1&v1=1&M1=100&dl1=0&ni1=1&nj1=0&nk1=0&th2=200&phi2=1&v2=-1&M2=100&dl2=0&ni2=1&nj2=0&nk2=0&th3=40&phi3=1&v3=2&M3=100&dl3=0&ni3=0&nj3=1&nk3=0&th4=280&phi4=1&v4=-1&M4=100&dl4=0&ni4=0&nj4=1&nk4=0&rc=1'
shot v11q4_old_default_rc1_auto0 "$D/index_v11q4.html?$OLD" 6000
# 既定（URL なし）＝正面衝突の配置、右は座標表示
shot v11q4_default_t0      "$D/index_v11q4.html?auto=0" 6000
shot v11q4_default_tick90  "$D/index_v11q4.html" 13500
dom  v11q4_default_tick90  "$D/index_v11q4.html" 13500
shot v11q4_headon_t0        "$D/index_v11q4.html?$HEAD&auto=0" 6000
shot v11q4_headon_auto6s    "$D/index_v11q4.html?$HEAD" 6000
shot v11q4_headon_tick90    "$D/index_v11q4.html?$HEAD" 13500   # tick 90：対 1–2 が Θ = 110° で重なる
echo done
