# K4 基底波＋倍音系 再現性監査レポート
Date: 2026-09-15

## 目的

今回の K4 完全ネットワーク実験について、保存済みの生データだけを入力にして、

- 静止画 PNG/SVG
- 辺データ CSV
- HTML アニメーション
- MP4 順回転 / 逆回転

を再生成できるか確認した。

## 生データの分離

元実験では数学的には等価だが浮動小数表現が異なる2種類の定義が使われていた。

### 静止画系

U = exp(2*pi*i/4)

このため 0 であるべき成分にも 6.123233995...e-17 などの丸め残差が含まれる。

### HTML / MP4 系

U = 1j

こちらは 0, +/-1, +/-i が厳密な離散値として生成される。
さらに -0.0 の符号が描画の完全一致に影響するため、生CSVを文字列から直接 float に読み出す方式を採用した。

したがって再現性確保のため、

- complex_states_static_exp.csv
- complex_states_animation_exact_i.csv

を別々に保存した。

## 波形生データ

元プログラムの描画点数まで再現するため、

- static: 2000 points, theta = -2*pi ... 2*pi
- HTML: 800 points, theta = 0 ... 2*pi
- MP4: 1200 points, theta = 0 ... 2*pi

を m=2 / m=3 それぞれCSV保存した。

## 静止画再現結果

以下6枚は保存済み元画像と pixel exact 一致した。

- waves_base_harmonic_m2.png
- waves_base_harmonic_m3.png
- composite_wave_m2.png
- composite_wave_m3.png
- complex_vectors_all_states_m2.png
- complex_vectors_all_states_m3.png

以下4枚は保存済み元画像とは pixel exact にならなかった。

- complex_vectors_theta_pi4_m2.png
- complex_vectors_theta_pi4_m3.png
- K4_complete_network_m2.png
- K4_complete_network_m3.png

ただし重要な点として、保存済みの元作図プログラム
generate_all_K4_harmonic_figures_20260915.py
を再実行しても、この4枚は歴史的に保存された元画像と一致しなかった。

一方、生データ読み出し専用プログラムの再生成画像は、
保存済み元作図プログラムを再実行した画像とは4枚すべて pixel exact 一致した。

したがって差の原因は生データ読み出し系ではなく、
当時の「保存済み画像」と「後から保存した作図プログラム」の版差である。

## MP4 再現結果

生データのみから再生成した4動画を ffmpeg framemd5 で全フレーム比較した。

- m=2 forward: 144 / 144 frames exact
- m=3 forward: 144 / 144 frames exact
- m=2 reverse: 144 / 144 frames exact
- m=3 reverse: 144 / 144 frames exact

したがって動画は全4本、デコード後フレームが完全一致した。

## HTML 再現結果

元HTMLと同一レイアウトテンプレートを使い、
状態値と波形を生データCSVから読み出す版を生成した。

ブラウザ自動スクリーンショットは環境管理ポリシーでローカルHTMLへのアクセスが禁止され実行できなかった。

代わりに元HTMLの数式生成座標と、生データ読み出し版の座標を直接比較した。

m=2 / m=3 について、

- base path
- harmonic path
- composite path
- k=0,1,2,3 の複素ベクトル
- 離散振幅マーカー

の最終描画座標はすべて一致した。

HTML波形は SVG path 作成時に小数点以下2桁へ丸められるが、
生データ読み出し版の path 文字列座標は元HTMLの path 座標と一致した。

## 結論

今回の実験は、

1. 完全な生データ
2. 生データ生成プログラム
3. 生データのみを入力にする静止画読み出しプログラム
4. 生データのみを入力にするHTML読み出しプログラム
5. 生データのみを入力にするMP4読み出しプログラム
6. 再生成された成果物
7. 再現性監査結果

をセットで保存することで再現可能になった。

特にMP4は全4本について全144フレーム完全一致を確認した。

歴史的静止画4枚のみ、当時保存した画像と後から保存した作図コードに版差がある。
この点は隠さず、本監査記録に固定する。
