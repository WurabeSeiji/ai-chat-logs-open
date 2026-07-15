# 系統B 灰色猫準安定界面R近傍スイープ 予備実験結果 v1

## 1. 実験条件

```text
coarse R range = 0.6 .. 0.9
coarse Delta R = 0.01
fine R range = 0.68 .. 0.71
fine Delta R = 0.001
R_137 = 0.6971778791282474
R_alpha128_nominal = 0.686671465671125
steps = 4096
C = off
D = off
```

## 2. 系統B 判定サマリー

| mode | A packet | B packet | N_A | N_B | kind | R*_metastable | metastable count | R*_eigen | eigen count | R*_selection | selection count | R*_joint | phase@joint | depth | S_amp | band5 | band10 | d137 | d128 | coarse R | fine R | control R |
|---|---|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| alternating_packet_2 | 1 | 1,3,4,5,6 | 1 | 6 | one_side_high_harmonic | 0.683 | 11 | 0.683 | 6 | 0.6 | 0 | 0.683 | gray_metastable | 5.90237109527 | 0.0199999824606 | 0.017 | 0.08 | 0.0141778791282 | 0.00367146567112 | 0.64 | 0.683 | 0 |

## 3. R地形ピーク候補

`R*_joint` は単一の代表点である。ただし系統BではR地形が多峰的になるため、局所ピークを別表として保存する。

ここでは、各ケースで `joint_gray_score` の局所最大点を順位付きで記録する。

`R_137` と `alpha~128 nominal` は、局所最大でない場合でも固定プローブ点として表に残す。

| case | rank | kind | R_peak | score | phase | depth | S_amp | d137 | d128 |
|---|---:|---|---:|---:|---|---:|---:|---:|---:|
| alternating_packet_2\|A=1\|B=1,3,4,5,6 | 1 | local_max | 0.683 | 10.1523710953 | gray_metastable | 5.90237109527 | 0.0199999824606 | 0.0141778791282 | 0.00367146567112 |
| alternating_packet_2\|A=1\|B=1,3,4,5,6 | 2 | local_max | 0.7 | 10.0091671605 | gray_metastable | 5.7591671605 | 0.019999996261 | 0.00282212087175 | 0.0133285343289 |
| alternating_packet_2\|A=1\|B=1,3,4,5,6 | 3 | local_max+probe_R137 | 0.697177879128 | 9.78658100097 | gray_metastable | 5.53658100097 | 0.0199999994916 | 2.47357689886e-13 | 0.0105064134569 |
| alternating_packet_2\|A=1\|B=1,3,4,5,6 | 4 | local_max | 0.64 | 9.43239500505 | gray_metastable | 5.18239500505 | 0.0199999879521 | 0.0571778791282 | 0.0466714656711 |
| alternating_packet_2\|A=1\|B=1,3,4,5,6 | 5 | local_max | 0.694 | 9.42691879981 | gray_metastable | 5.17691879981 | 0.0199999362998 | 0.00317787912825 | 0.00732853432887 |
| alternating_packet_2\|A=1\|B=1,3,4,5,6 | 6 | local_max | 0.66 | 9.40048851059 | gray_metastable | 5.15048851059 | 0.0199999990794 | 0.0371778791282 | 0.0266714656711 |
| alternating_packet_2\|A=1\|B=1,3,4,5,6 | 7 | local_max | 0.696 | 9.24201198105 | gray_metastable | 4.99201198105 | 0.0199999830725 | 0.00117787912825 | 0.00932853432887 |
| alternating_packet_2\|A=1\|B=1,3,4,5,6 | 8 | local_max | 0.708 | 9.16910092528 | gray_metastable | 4.91910092528 | 0.0199999781426 | 0.0108221208718 | 0.0213285343289 |
| alternating_packet_2\|A=1\|B=1,3,4,5,6 | 9 | local_max | 0.72 | 9.14607894579 | gray_metastable | 4.89607894579 | 0.019999206968 | 0.0228221208718 | 0.0333285343289 |
| alternating_packet_2\|A=1\|B=1,3,4,5,6 | 10 | local_max | 0.692 | 9.01729950252 | gray_metastable | 4.76729950252 | 0.0199999587608 | 0.00517787912825 | 0.00532853432887 |
| alternating_packet_2\|A=1\|B=1,3,4,5,6 | 11 | local_max | 0.686 | 8.98354031964 | gray_metastable | 4.73354031964 | 0.0199999663371 | 0.0111778791282 | 0.000671465671125 |
| alternating_packet_2\|A=1\|B=1,3,4,5,6 | 12 | local_max | 0.76 | 8.96019159712 | gray_metastable | 4.71019159712 | 0.019999989844 | 0.0628221208718 | 0.0733285343289 |
| alternating_packet_2\|A=1\|B=1,3,4,5,6 | 13 | probe_alpha128_nominal | 0.686671465671 | 8.58203865501 | gray_metastable | 4.33203865501 | 0.0199999769547 | 0.0105064134572 | 1.25011112573e-13 |

## 4. 読み

本予備実験では、20260714 の `S=p_A-p_B` 準安定力学を母体とし、AB交換を `R` 指定の散乱写像に置き換えた。

Nと倍音パケットは、系統Aと同じ比較軸として記録した。

この二成分準安定モデルでは、Nと奇偶は `S` 写像に直接入らないため、N依存性が出る場合は、R写像または条件選別を通じた間接効果として読む。

`R_137` は低エネルギー側の精密プローブ、`alpha~128 nominal` は高エネルギー側の名目プローブである。

`alpha~128 nominal` は精密値ではないため、固定点そのものだけでなく、その近傍局所ピークをあわせて読む。

したがって、最深ピークが別のRにある場合でも、alpha近傍に独立したピークが立つなら、それ自体をR地形上の重要な候補として扱う。

分布図では、`gray_error=|S_mean|+|S_amp-S_amp_target|+S_drift+penalty` を `gray_depth=-log10(gray_error)` として縦軸に置く。

灰色点は各 `R` における各初期条件、黒線は各 `R` で最も深い条件の包絡線である。

## 5. 出力

| 種類 | ファイル |
|---|---|
| rows | `system_B_alternating_packet_2_A1_B1_3_4_5_6_rows_v1.csv` |
| samples | `system_B_alternating_packet_2_A1_B1_3_4_5_6_samples_v1.csv` |
| summary | `system_B_alternating_packet_2_A1_B1_3_4_5_6_summary_v1.csv` |
| best | `system_B_alternating_packet_2_A1_B1_3_4_5_6_best_v1.csv` |
| peaks | `system_B_alternating_packet_2_A1_B1_3_4_5_6_peaks_v1.csv` |
| result | `system_B_alternating_packet_2_A1_B1_3_4_5_6_result_v1.json` |
| report | `system_B_alternating_packet_2_A1_B1_3_4_5_6_report_v1.md` |
| scores | `system_B_alternating_packet_2_A1_B1_3_4_5_6_scores_v1.png` |
| gray_depth_overview | `system_B_alternating_packet_2_A1_B1_3_4_5_6_gray_depth_distribution_overview_v1.png` |
| gray_depth_deep | `system_B_alternating_packet_2_A1_B1_3_4_5_6_gray_depth_distribution_deep_v1.png` |
| peak_zoom | `system_B_alternating_packet_2_A1_B1_3_4_5_6_peak_zoom_v1.png` |
