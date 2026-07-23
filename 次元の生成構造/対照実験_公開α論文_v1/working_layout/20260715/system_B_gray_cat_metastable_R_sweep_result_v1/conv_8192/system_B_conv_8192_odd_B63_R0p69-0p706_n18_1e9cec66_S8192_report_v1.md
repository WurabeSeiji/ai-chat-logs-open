# 系統B 灰色猫準安定界面R近傍スイープ 予備実験結果 v1

## 1. 実験条件

```text
coarse R range = 0.6 .. 0.9
coarse Delta R = 0.01
fine R range = 0.68 .. 0.71
fine Delta R = 0.001
R_137 = 0.6971778791282474
R_alpha128_nominal = 0.686671465671125
steps = 8192
C = off
D = off
```

## 2. 系統B 判定サマリー

| mode | A packet | B packet | N_A | N_B | kind | R*_metastable | metastable count | R*_eigen | eigen count | R*_selection | selection count | R*_joint | phase@joint | depth | S_amp | band5 | band10 | d137 | d128 | coarse R | fine R | control R |
|---|---|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| odd_kernel | 1 | 63 | 1 | 63 | one_side_high_harmonic | 0.69 | 9 | 0.69 | 8 | 0.69 | 0 | 0.69 | gray_metastable | 5.36660722429 | 0.0199999565391 | 0.013 | 0.016 | 0.00717787912825 | 0.00332853432887 | nan | 0.69 | nan |

## 3. R地形ピーク候補

`R*_joint` は単一の代表点である。ただし系統BではR地形が多峰的になるため、局所ピークを別表として保存する。

ここでは、各ケースで `joint_gray_score` の局所最大点を順位付きで記録する。

`R_137` と `alpha~128 nominal` は、局所最大でない場合でも固定プローブ点として表に残す。

| case | rank | kind | R_peak | score | phase | depth | S_amp | d137 | d128 |
|---|---:|---|---:|---:|---|---:|---:|---:|---:|
| odd_kernel\|A=1\|B=63 | 1 | local_max | 0.69 | 9.61660722429 | gray_metastable | 5.36660722429 | 0.0199999565391 | 0.00717787912825 | 0.00332853432887 |
| odd_kernel\|A=1\|B=63 | 2 | local_max | 0.703 | 9.46476696429 | gray_metastable | 5.21476696429 | 0.0199999967665 | 0.00582212087175 | 0.0163285343289 |
| odd_kernel\|A=1\|B=63 | 3 | local_max | 0.696 | 9.40271042234 | gray_metastable | 5.15271042234 | 0.0199999990155 | 0.00117787912825 | 0.00932853432887 |
| odd_kernel\|A=1\|B=63 | 4 | local_max+probe_R137 | 0.697177879128 | 9.22734215743 | gray_metastable | 4.97734215743 | 0.0199999980542 | 2.47357689886e-13 | 0.0105064134569 |
| odd_kernel\|A=1\|B=63 | 5 | local_max | 0.7 | 9.20291452678 | gray_metastable | 4.95291452678 | 0.019999997138 | 0.00282212087175 | 0.0133285343289 |
| odd_kernel\|A=1\|B=63 | 6 | local_max | 0.693 | 9.11042333521 | gray_metastable | 4.86042333521 | 0.0199999789862 | 0.00417787912825 | 0.00632853432887 |
| odd_kernel\|A=1\|B=63 | 7 | local_max | 0.706 | 8.8695201982 | gray_metastable | 4.6195201982 | 0.0199999965335 | 0.00882212087175 | 0.0193285343289 |

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
| rows | `system_B_conv_8192_odd_B63_R0p69-0p706_n18_1e9cec66_S8192_rows_v1.csv` |
| samples | `system_B_conv_8192_odd_B63_R0p69-0p706_n18_1e9cec66_S8192_samples_v1.csv` |
| summary | `system_B_conv_8192_odd_B63_R0p69-0p706_n18_1e9cec66_S8192_summary_v1.csv` |
| best | `system_B_conv_8192_odd_B63_R0p69-0p706_n18_1e9cec66_S8192_best_v1.csv` |
| peaks | `system_B_conv_8192_odd_B63_R0p69-0p706_n18_1e9cec66_S8192_peaks_v1.csv` |
| result | `system_B_conv_8192_odd_B63_R0p69-0p706_n18_1e9cec66_S8192_result_v1.json` |
| report | `system_B_conv_8192_odd_B63_R0p69-0p706_n18_1e9cec66_S8192_report_v1.md` |
