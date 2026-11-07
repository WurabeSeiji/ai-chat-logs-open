#!/bin/sh
# 描画記録（ヘッドレス Chrome）と、原本 index_v4q3.html との DOM 対照（右上の読める量の文字列）を取る。
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
for pair in "default:$Q_DEF" "tri:$Q_TRI" "mirror:$Q_MIR" "std:$Q_STD" "tri_ijk:$Q_TRI&ni1=1&nk1=0&nj2=1&nk2=0" "mirror_k_i_mi:$Q_MIR&ni2=1&nk2=0&ni3=-1&nk3=0"; do
  name=${pair%%:*}; q=${pair#*:}
  dom "v4q3_$name"      "$D/原本/index_v4q3.html?$q" 6000
  dom "v4q3_pow1_$name" "$D/原本/index_v4q3.html?$q&pow=1" 6000
  dom "v5q3_$name"      "$D/index_v5q3.html?$q" 6000
  dom "v5q3_pow1_$name" "$D/index_v5q3.html?$q&pow=1" 6000
  dom "v5q3_pow2_$name" "$D/index_v5q3.html?$q&pow=2" 6000
done
dom "v4q3_pow1_auto6s" "$D/原本/index_v4q3.html?pow=1" 6000
dom "v5q3_pow1_auto6s" "$D/index_v5q3.html?pow=1" 6000
dom "v5q3_pow2_auto6s" "$D/index_v5q3.html?pow=2" 6000
shot v5q3_pow1_tri_ijk "$D/index_v5q3.html?$Q_TRI&ni1=1&nk1=0&nj2=1&nk2=0&pow=1" 6000
shot v5q3_pow2_tri_ijk "$D/index_v5q3.html?$Q_TRI&ni1=1&nk1=0&nj2=1&nk2=0&pow=2" 6000
shot v5q3_pow1_auto6s  "$D/index_v5q3.html?pow=1" 6000
shot v5q3_pow2_auto6s  "$D/index_v5q3.html?pow=2" 6000
echo done
