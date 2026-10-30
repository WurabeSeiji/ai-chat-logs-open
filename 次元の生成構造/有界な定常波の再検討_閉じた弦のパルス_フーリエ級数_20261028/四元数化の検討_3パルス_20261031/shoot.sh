#!/bin/sh
# 描画記録（ヘッドレス Chrome）と、原本 index_v2c3.html との DOM 対照（読める量の文字列）を取る。
# 使い方: sh shoot.sh        （Chrome は /Applications/Google Chrome.app、Plotly は cdnjs から読むので通信が要る）
set -e
cd "$(dirname "$0")"
CH="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
D="file://$(pwd)"
shot() { "$CH" --headless=new --hide-scrollbars --window-size=1800,1100 --virtual-time-budget=6000 --screenshot="figs/$1.png" "$2" >/dev/null 2>&1; echo "shot $1"; }
dom()  { "$CH" --headless=new --virtual-time-budget=6000 --dump-dom "$2" 2>/dev/null > "results/dom_$1.html"; echo "dom  $1 ($(wc -c < "results/dom_$1.html") bytes)"; }

# ---- 対照：同じ URL 引数で原本（複素）と四元数版（n は既定の k）。読める量の文字列を control_html.py で突き合わせる
Q_DEF='auto=0'
Q_TRI='th1=180&phi1=10&dl1=0&th2=150&phi2=10&dl2=0&th3=100&phi3=10&dl3=0&auto=0'
Q_MIR='th1=180&th2=200&th3=160&phi1=10&phi2=10&phi3=10&dl1=0&dl2=0&dl3=0&p02=45&p03=-45&auto=0'
Q_STD='th3=90&unit=1&auto=0'
for pair in "default:$Q_DEF" "tri:$Q_TRI" "mirror:$Q_MIR" "std:$Q_STD"; do
  name=${pair%%:*}; q=${pair#*:}
  dom "c3_$name" "$D/原本/index_v2c3.html?$q"
  dom "q3_$name" "$D/index_v2q3.html?$q"
done

# ---- 描画記録
shot v2q3_D0        "$D/index_v2q3.html?phi1=10&dl1=0&dl2=0&dl3=0&auto=0"                                   # 原本 v2c3_D0 と同じ配置（n は全部 k）
shot v2q3_tri_kkk   "$D/index_v2q3.html?$Q_TRI"                                                             # 三角形、軸 k, k, k（原本 v2c3_tri と同じ）
shot v2q3_tri_ijk   "$D/index_v2q3.html?$Q_TRI&ni1=1&nk1=0&nj2=1&nk2=0"                                     # 三角形、軸 i, j, k
shot v2q3_same_ijk  "$D/index_v2q3.html?th1=150&th2=150&th3=150&phi1=10&dl1=0&dl2=0&dl3=0&auto=0&ni1=1&nk1=0&nj2=1&nk2=0"   # 同じ位置、軸 i, j, k
shot v2q3_mirror_k_i_i  "$D/index_v2q3.html?$Q_MIR&ni2=1&nk2=0&ni3=1&nk3=0"                                 # 鏡映配置、対が i を共有（Σc² 実）
shot v2q3_mirror_k_i_mi "$D/index_v2q3.html?$Q_MIR&ni2=1&nk2=0&ni3=-1&nk3=0"                                # 鏡映配置、対の軸が逆 i, −i（Σc² 閉じない）。ファイル名は大文字小文字を使わない（macOS は区別しない）
shot v2q3_cone      "$D/index_v2q3.html"                                                                    # 既定（Auto）、軸 k, k, k
shot v2q3_cone_ijk  "$D/index_v2q3.html?ni1=1&nk1=0&nj2=1&nk2=0"                                            # 既定（Auto）、軸 i, j, k
echo done
