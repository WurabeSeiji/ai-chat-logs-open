#!/bin/sh
# 描画記録（ヘッドレス Chrome、Auto のものは 180 秒で打ち切り）と、原本 index_v7q3.html との DOM 対照。配置は M・k を含めて全部 URL で明示
set -e
cd "$(dirname "$0")"
CH="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
D="file://$(pwd)"
shot() { perl -e 'alarm shift; exec @ARGV' 180 "$CH" --headless=new --hide-scrollbars --window-size=1800,1100 --virtual-time-budget=$3 --screenshot="figs/$1.png" "$2" >/dev/null 2>&1 || echo "shot $1 timed out"; echo "shot $1"; }
dom()  { perl -e 'alarm shift; exec @ARGV' 180 "$CH" --headless=new --virtual-time-budget=$3 --dump-dom "$2" 2>/dev/null > "results/dom_$1.html" || true; echo "dom  $1 ($(wc -c < "results/dom_$1.html") bytes)"; }
Q_DEF='th1=180&phi1=20&M1=10&dl1=20&th2=180&phi2=10&M2=10&dl2=0&th3=180&phi3=10&M3=10&dl3=10&k=10&auto=0'
Q_TRI='th1=180&phi1=10&M1=10&dl1=0&th2=150&phi2=10&M2=10&dl2=0&th3=100&phi3=10&M3=10&dl3=0&k=10&auto=0'
for pair in "default:$Q_DEF" "tri:$Q_TRI"; do
  name=${pair%%:*}; q=${pair#*:}
  for pw in 0 1 2; do for rc in 0 1; do
    dom "v7q3_pow${pw}_rc${rc}_$name" "$D/原本/index_v7q3.html?$q&pow=$pw&rc=$rc" 6000
    dom "v8q3_pow${pw}_rc${rc}_$name" "$D/index_v8q3.html?$q&pow=$pw&rc=$rc" 6000
  done; done
done
dom "v8q3_newdefault_rc1" "$D/index_v8q3.html?auto=0&rc=1" 6000
dom "v8q3_auto6s_rc1"     "$D/index_v8q3.html?rc=1" 6000
dom "v8q3_auto6s_ijk_rc1" "$D/index_v8q3.html?rc=1&ni1=1&nk1=0&nj2=1&nk2=0" 6000
shot v8q3_default_rc1 "$D/index_v8q3.html?auto=0&rc=1" 6000
shot v8q3_auto6s_rc1  "$D/index_v8q3.html?rc=1" 6000
shot v8q3_auto6s_ijk_rc1 "$D/index_v8q3.html?rc=1&ni1=1&nk1=0&nj2=1&nk2=0" 6000
echo done
