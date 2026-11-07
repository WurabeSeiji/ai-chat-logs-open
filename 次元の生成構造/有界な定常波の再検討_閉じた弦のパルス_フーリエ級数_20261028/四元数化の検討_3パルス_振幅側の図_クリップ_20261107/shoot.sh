#!/bin/sh
# 描画記録（ヘッドレス Chrome）と、原本 index_v5q3.html との DOM 対照（右上の読める量の文字列）を取る。配置は M・k を含めて全部 URL で明示（v7q3 は初期値を変えたので）
set -e
cd "$(dirname "$0")"
CH="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
D="file://$(pwd)"
shot() { "$CH" --headless=new --hide-scrollbars --window-size=1800,1100 --virtual-time-budget=$3 --screenshot="figs/$1.png" "$2" >/dev/null 2>&1; echo "shot $1"; }
dom()  { "$CH" --headless=new --virtual-time-budget=$3 --dump-dom "$2" 2>/dev/null > "results/dom_$1.html"; echo "dom  $1 ($(wc -c < "results/dom_$1.html") bytes)"; }
Q_DEF='th1=180&phi1=20&M1=10&dl1=20&th2=180&phi2=10&M2=10&dl2=0&th3=180&phi3=10&M3=10&dl3=10&k=10&auto=0'
Q_TRI='th1=180&phi1=10&M1=10&dl1=0&th2=150&phi2=10&M2=10&dl2=0&th3=100&phi3=10&M3=10&dl3=0&k=10&auto=0'
Q_MIR='th1=180&th2=200&th3=160&phi1=10&phi2=10&phi3=10&M1=10&M2=10&M3=10&dl1=0&dl2=0&dl3=0&p02=45&p03=-45&k=10&auto=0'
Q_STD='th1=180&phi1=20&M1=10&dl1=20&th2=180&phi2=10&M2=10&dl2=0&th3=90&phi3=10&M3=10&dl3=10&k=10&unit=1&auto=0'
for pair in "default:$Q_DEF" "tri:$Q_TRI" "tri_ijk:$Q_TRI&ni1=1&nk1=0&nj2=1&nk2=0" "mirror:$Q_MIR" "mirror_k_i_mi:$Q_MIR&ni2=1&nk2=0&ni3=-1&nk3=0" "std:$Q_STD"; do
  name=${pair%%:*}; q=${pair#*:}
  for pw in 0 1 2; do
    dom "v5q3_pow${pw}_$name" "$D/原本/index_v5q3.html?$q&pow=$pw" 6000
    dom "v7q3_pow${pw}_$name" "$D/index_v7q3.html?$q&pow=$pw" 6000
  done
done
Q_OLD='th1=180&phi1=20&M1=10&dl1=20&th2=180&phi2=10&M2=10&dl2=0&th3=180&phi3=10&M3=10&dl3=10&k=10&pow=0'
dom "v5q3_pow0_auto6s" "$D/原本/index_v5q3.html?$Q_OLD" 6000
dom "v7q3_pow0_auto6s" "$D/index_v7q3.html?$Q_OLD" 6000
dom "v7q3_newdefault" "$D/index_v7q3.html?auto=0" 6000
shot v7q3_default      "$D/index_v7q3.html?auto=0" 6000
shot v7q3_default_4c   "$D/index_v7q3.html?auto=0&pow=0" 6000
shot v7q3_default_ijk  "$D/index_v7q3.html?auto=0&ni1=1&nk1=0&nj2=1&nk2=0" 6000
shot v7q3_tri_ijk      "$D/index_v7q3.html?$Q_TRI&ni1=1&nk1=0&nj2=1&nk2=0" 6000
shot v7q3_tri_ijk_pow2 "$D/index_v7q3.html?$Q_TRI&ni1=1&nk1=0&nj2=1&nk2=0&pow=2" 6000
shot v7q3_mirror       "$D/index_v7q3.html?$Q_MIR" 6000
echo done
