#!/bin/sh
# 撮影（ヘッドレス Chrome、180 秒で打ち切り）：座標表示 rc=1、auto=0 の静止画を原本 v9q4 と v10q4 で撮る。変わるのは球の半径だけ
set -e
cd "$(dirname "$0")"
CH="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
D="file://$(pwd)"
shot() { perl -e 'alarm shift; exec @ARGV' 180 "$CH" --headless=new --hide-scrollbars --window-size=1800,1400 --virtual-time-budget=$3 --screenshot="figs/$1.png" "$2" >/dev/null 2>&1 || echo "shot $1 timed out"; echo "shot $1"; }
shot v9q4_rc1_auto0  "$D/原本/index_v9q4.html?auto=0&rc=1" 6000
shot v10q4_rc1_auto0 "$D/index_v10q4.html?auto=0&rc=1" 6000
shot v10q4_rc1_auto6s "$D/index_v10q4.html?rc=1" 6000
echo done
