#!/bin/sh
# 描画記録（ヘッドレス Chrome）と、原本 index_v2q3.html との DOM 対照（右上の読める量の文字列）を取る。
# 使い方: sh shoot.sh        （Chrome は /Applications/Google Chrome.app、Plotly は cdnjs から読むので通信が要る）
set -e
cd "$(dirname "$0")"
CH="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
D="file://$(pwd)"
shot() { "$CH" --headless=new --hide-scrollbars --window-size=1800,1100 --virtual-time-budget=$3 --screenshot="figs/$1.png" "$2" >/dev/null 2>&1; echo "shot $1"; }
dom()  { "$CH" --headless=new --virtual-time-budget=$3 --dump-dom "$2" 2>/dev/null > "results/dom_$1.html"; echo "dom  $1 ($(wc -c < "results/dom_$1.html") bytes)"; }

# ---- 対照：同じ URL 引数で原本と v3q3。読める量は同じで、v3q3 の h1 に「枠 ±L（いま必要な枠 ±L′）」だけが足される
Q_DEF='auto=0'
Q_TRI='th1=180&phi1=10&dl1=0&th2=150&phi2=10&dl2=0&th3=100&phi3=10&dl3=0&auto=0'
Q_MIR='th1=180&th2=200&th3=160&phi1=10&phi2=10&phi3=10&dl1=0&dl2=0&dl3=0&p02=45&p03=-45&auto=0'
Q_STD='th3=90&unit=1&auto=0'
for pair in "default:$Q_DEF" "tri:$Q_TRI" "mirror:$Q_MIR" "std:$Q_STD"; do
  name=${pair%%:*}; q=${pair#*:}
  dom "v2q3_$name" "$D/原本/index_v2q3.html?$q" 6000
  dom "v3q3_$name" "$D/index_v3q3.html?$q" 6000
done
# ---- Auto を仮想時間 6 秒（約 40 ティック）走らせた後の DOM：v3q3 の「枠」が auto=0 の値のまま、「いま必要な枠」だけ変わる
dom "v3q3_auto6s" "$D/index_v3q3.html" 6000

# ---- 描画記録：既定（auto=0）、Auto 6 秒後（原本と v3q3、同じ仮想時間）
shot v3q3_default "$D/index_v3q3.html?auto=0" 6000
shot v2q3_auto6s  "$D/原本/index_v2q3.html" 6000
shot v3q3_auto6s  "$D/index_v3q3.html" 6000
echo done
