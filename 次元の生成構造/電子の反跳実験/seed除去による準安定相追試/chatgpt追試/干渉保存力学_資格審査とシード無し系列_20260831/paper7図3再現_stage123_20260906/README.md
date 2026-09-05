# 論文7 図3（5色占有）の段1+2+3 データによる再現（読み出しのみ）

作成日: 2026-09-06

## 目的

種なし系論文（第8論文系、note「タネがなくても、インフレーションは起きた」）の結論
「第3の方向まで含めた三方向構造が急拡大の停止とともに定着」が、旧エンジン＋5色読出しの
アーティファクトかどうかを、同じ読出しを今回のデータ（段1+2+3、10,000歩、N=5, 40）に
適用して判定する。

## 方法（系列規約準拠）

- 基底構成（make_parent / parent_plane_split_exact / gram_reduce / dominant_plane /
  zero_closure_kernel_seed）は第5論文原本エンジンからそのまま import。
- 5色読出し（occ / s4_new_dirs / align_2d、E_d3・E_d4・残余・核の式）は
  run_paper7_5color_timeseries.py の逐語コピー。
- 時間発展は走らせず、今回の 10,000歩状態 npz（SHA台帳照合）を読むだけ。
- **対照ゲート全合格**: 入力 npz SHA 一致、make_parent の v・初期 Z が今回の静的親 npz
  および states[0] と bit 一致（旧フレームワークと今回のフレームワークの親が同一である
  ことの直接証明）。

## 結果: 第3の方向はアーティファクトではない——今回のエンジンでも再現

fig3_compare_stage123_N5_N40.png（正本図と同一様式・床1e-6）:

- **N=40**: crossing（step 358）とともに direction 3（赤）と direction 4（橙）が床から
  立ち上がり 0.022 / 0.025 で定着。kernel（緑）0.062、P1（青）0.89。
  **remaining other-rotation（灰）は床（10⁻⁶）へ落ちたまま**——成長は新2方向＋核に
  完全に閉じ込められ、other 空間の残りには漏れない。終端 f = 0.109。
- **N=5**: d3 = 0.11、d4 = 0.075、kernel 0.21、P1 0.55、f ≈ 0.40。周辺系らしい
  ゆらぎを伴うが方向構造は同型。灰は即座に床。
- 旧論文の図（paper7_f_projection_v1/figures_control）と定性的に同一の構造
  （P1高・d3≈d4 の立ち上がりと定着・灰の消滅・核の定常帯）が、**全く別の力学
  （σ時計+Cayley → 固定時計スペクトル写像）で再現**された。

判定: 「三方向構造の定着」は読出し・エンジンに依存しない実在の構造である。
第3・第4方向（=新しい1枚の平面）は確かに立ち上がり、それ以外の方向はほぼ空のまま
——今回の低ランク成長（支配平面ランク2）・位相ガラスの知見とも整合する。

## ファイル

- make_fig3_5color_stage123_v1.py（ゲート内蔵）
- fig3_compare_stage123_N5_N40.png / .svg
- bands5_stage123_N00005.csv / bands5_stage123_N00040.csv
- make_fig3_5color_stage123_v1_meta.json
- run_all.sh / SHA256SUMS.txt
