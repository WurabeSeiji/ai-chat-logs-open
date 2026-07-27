# 系統A 局在性交換R近傍斉一スイープ 予備実験結果 v1

## 1. 実験条件

```text
coarse R range = 0.6 .. 0.9
coarse Delta R = 0.01
fine R range = 0.68 .. 0.71
fine Delta R = 0.001
R_137 = 0.6971778791282474
R_128 = 0.686671465671125
max_collision = 256
pairs = ((1, 1), (1, 2), (1, 3), (1, 5), (1, 10), (1, 15), (1, 31), (1, 63), (31, 31), (63, 63))
cases = ['custom_packet|A=1|B=1,2;p=0,0.628318530718']
```

## 2. 系統A 判定サマリー

| mode | A packet | B packet | A weights | B weights | A phases | B phases | A wavelengths | B wavelengths | A shift | B shift | N_A | N_B | kind | R*_L | col_L | R*_N | col_N | R*_transfer | col_transfer | R*_joint | col_joint | band5 | band10 | d137 | d128 | coarse R | fine R | control R |
|---|---|---|---|---|---|---|---|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| custom_packet | 1 | 1,2 | 1 | 1,1 | 0 | 0,0.628318530718 | 1 | 1,1 | 0 | 0 | 1 | 2 | one_side_high_harmonic | 0.697177879128 | 31 | 0.697177879128 | 31 | 0.697177879128 | 62 | 0.697177879128 | 31 | 0 | 0 | 2.47357689886e-13 | 0.0105064134569 | 0.74 | 0.697177879128 | 0.5 |

## 3. 読み

本予備実験では、片側高次倍音条件と同次数対照を同じ `R` 掃引で比較した。

`R_star_joint` は、広域粗スイープと局所精密スイープを合わせた候補範囲の中で、局在性差、実効次数差、B側初期倍音分布のA側への類似度を合わせて読んだ値である。

`R=0.0`, `R=0.5`, `R=1.0` は対照であり、`R_star` 判定からは除外した。

同次数対照では `R_band_width` が主掃引幅全体に広がり、`R` に対して鋭い判定点を持たない。

一方、片側高次倍音条件では、`R_star_L`, `R_star_N`, `R_star_transfer`, `R_star_joint` が同じ位置に集まった。

ただし、広域粗スイープでは複数の浅い谷候補も見える。

特に片側高次倍音条件では、粗スイープ上で `R=0.74` 付近が最良の粗候補となり、`R=0.63`, `R=0.67`, `R=0.78`, `R=0.86`, `R=0.89` 付近にも浅い谷が現れた。

したがって、`R_137` 近傍だけを細分化するのではなく、これらの谷候補も次段階の局所精密スイープ対象に含める。

`R_137` と `R_128` は判定値ではなく、全系統で同じ位置を読むための固定プローブ点である。

本予備実験では `R_137` を明示的にプローブ点として含めたため、次段階では `R_137` 近傍を特別扱いせず、より細かい一様掃引で谷幅を確認する必要がある。

また、点としての `R_star` だけでは、反射係数の最小点と衝突回数方向の干渉縞を分離できない。

そのため、本実装では `R` ごとの全衝突回を縦方向へ並べる分布図を追加した。

分布図では、`L_gap` と `N_eff_gap` から作った `gap_terrain_score` を `depth=-log10(gap_terrain_score)` として縦軸に置く。

完全同一対照では数値的に深さが過大になるため、図の表示では `depth` を 6 で上限表示する。CSV には上限前の `gap_depth` も保存する。

灰色点は各 `R` における全衝突回の分布、黒線は各 `R` で最も深い衝突回の包絡線である。

最深部の拡大分布図は、次段階で局所スイープ幅を決めるための候補図である。

## 4. 出力

| 種類 | ファイル |
|---|---|
| rows | `system_A_phase10_B12_custom_packet_A1_B1-2-p-0-0-628318530718_Rdefault_C256_rows_v1.csv` |
| summary | `system_A_phase10_B12_custom_packet_A1_B1-2-p-0-0-628318530718_Rdefault_C256_summary_v1.csv` |
| best | `system_A_phase10_B12_custom_packet_A1_B1-2-p-0-0-628318530718_Rdefault_C256_best_v1.csv` |
| terrain | `system_A_phase10_B12_custom_packet_A1_B1-2-p-0-0-628318530718_Rdefault_C256_collision_terrain_v1.csv` |
| json | `system_A_phase10_B12_custom_packet_A1_B1-2-p-0-0-628318530718_Rdefault_C256_result_v1.json` |
| report | `system_A_phase10_B12_custom_packet_A1_B1-2-p-0-0-628318530718_Rdefault_C256_report_v1.md` |
