#!/bin/sh
# 描画記録（ヘッドレス Chrome）と、原本 index_v6c3.html との DOM 対照（右上の読める量の文字列）を取る。
set -e
cd "$(dirname "$0")"
CH="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
D="file://$(pwd)"
shot() { "$CH" --headless=new --hide-scrollbars --window-size=1800,1100 --virtual-time-budget=$3 --screenshot="figs/$1.png" "$2" >/dev/null 2>&1; echo "shot $1"; }
dom()  { "$CH" --headless=new --virtual-time-budget=$3 --dump-dom "$2" 2>/dev/null > "results/dom_$1.html"; echo "dom  $1 ($(wc -c < "results/dom_$1.html") bytes)"; }
Q_DEF='auto=0'
Q_TRI='th1=180&phi1=10&dl1=0&th2=150&phi2=10&dl2=0&th3=100&phi3=10&dl3=0&auto=0'
Q_MIR='th1=180&th2=200&th3=160&phi1=10&phi2=10&phi3=10&dl1=0&dl2=0&dl3=0&p02=45&p03=-45&auto=0'
Q_STD='th3=90&unit=1&auto=0'
for pair in "default:$Q_DEF" "tri:$Q_TRI" "mirror:$Q_MIR" "std:$Q_STD"; do
  name=${pair%%:*}; q=${pair#*:}
  for pw in 0 1 2; do
    dom "v6c3_pow${pw}_$name" "$D/原本/index_v6c3.html?$q&pow=$pw" 6000
    dom "v7c3_pow${pw}_$name" "$D/index_v7c3.html?$q&pow=$pw" 6000
  done
done
dom "v6c3_pow0_auto6s" "$D/原本/index_v6c3.html" 6000
dom "v7c3_pow0_auto6s" "$D/index_v7c3.html" 6000
shot v7c3_default   "$D/index_v7c3.html?auto=0" 6000
shot v7c3_tri       "$D/index_v7c3.html?$Q_TRI" 6000
shot v7c3_tri_pow1  "$D/index_v7c3.html?$Q_TRI&pow=1" 6000
shot v7c3_tri_pow2  "$D/index_v7c3.html?$Q_TRI&pow=2" 6000
shot v7c3_mirror    "$D/index_v7c3.html?$Q_MIR" 6000
echo done
