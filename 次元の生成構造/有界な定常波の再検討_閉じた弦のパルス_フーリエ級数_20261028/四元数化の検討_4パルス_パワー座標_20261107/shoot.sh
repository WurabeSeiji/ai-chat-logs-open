#!/bin/sh
# 描画記録（ヘッドレス Chrome、180 秒で打ち切り）と、原本 index_v8q3.html との DOM 対照。
# 対照は A₄ = 0・n₄ = k（4 が和に寄与せず、巻き数の軸判定も変えない）で、1〜3 の配置は M・k・n を含めて全部 URL で明示（n は全部 k）
set -e
cd "$(dirname "$0")"
CH="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
D="file://$(pwd)"
shot() { perl -e 'alarm shift; exec @ARGV' 180 "$CH" --headless=new --hide-scrollbars --window-size=1800,1400 --virtual-time-budget=$3 --screenshot="figs/$1.png" "$2" >/dev/null 2>&1 || echo "shot $1 timed out"; echo "shot $1"; }
dom()  { perl -e 'alarm shift; exec @ARGV' 180 "$CH" --headless=new --virtual-time-budget=$3 --dump-dom "$2" 2>/dev/null > "results/dom_$1.html" || true; echo "dom  $1 ($(wc -c < "results/dom_$1.html") bytes)"; }
Q_DEF='th1=180&phi1=20&M1=10&dl1=20&th2=180&phi2=10&M2=10&dl2=0&th3=180&phi3=10&M3=10&dl3=10&ni1=0&nj1=0&nk1=1&ni2=0&nj2=0&nk2=1&ni3=0&nj3=0&nk3=1&k=10&auto=0'
Q_TRI='th1=180&phi1=10&M1=10&dl1=0&th2=150&phi2=10&M2=10&dl2=0&th3=100&phi3=10&M3=10&dl3=0&ni1=0&nj1=0&nk1=1&ni2=0&nj2=0&nk2=1&ni3=0&nj3=0&nk3=1&k=10&auto=0'
Q4='th4=300&phi4=1&M4=100&dl4=-10&A4=0&ni4=0&nj4=0&nk4=1'
for pair in "default:$Q_DEF" "tri:$Q_TRI"; do
  name=${pair%%:*}; q=${pair#*:}
  for pw in 0 1 2; do for rc in 0 1; do
    dom "v8q3_pow${pw}_rc${rc}_$name" "$D/原本/index_v8q3.html?$q&pow=$pw&rc=$rc" 6000
    dom "v9q4_pow${pw}_rc${rc}_$name" "$D/index_v9q4.html?$q&$Q4&pow=$pw&rc=$rc" 6000
  done; done
done
dom "v9q4_newdefault_rc1" "$D/index_v9q4.html?auto=0&rc=1" 6000
dom "v9q4_newdefault_rc0" "$D/index_v9q4.html?auto=0" 6000
dom "v9q4_auto6s_rc1"     "$D/index_v9q4.html?rc=1" 6000
dom "v9q4_auto6s_rc0"     "$D/index_v9q4.html" 6000
shot v9q4_default_rc0 "$D/index_v9q4.html?auto=0" 6000
shot v9q4_default_rc1 "$D/index_v9q4.html?auto=0&rc=1" 6000
shot v9q4_auto6s_rc1  "$D/index_v9q4.html?rc=1" 6000
echo done
