# 白猫黒猫灰色猫準安定界面におけるC弱読出しとD強観測選択予備実験総括 v1

**日付:** 2026-07-14  
**著者:** 木原 範昭  
**位置づけ:** 波の情報読出しシリーズ・閉鎖系内A/B配分準安定状態と観測選択の予備実験  
**Version DOI:** 10.5281/zenodo.21353209  
**Concept DOI:** 10.5281/zenodo.21353208

---

## 0. 結論

本稿は、白猫状態 `A` と黒猫状態 `B` の配分を、閉鎖系内の複素振幅

```math
a,\quad b
```

として表し、

```math
p_A=|a|^2,\qquad p_B=|b|^2
```

```math
Q=p_A+p_B,\qquad S=p_A-p_B
```

を読む予備実験をまとめる。

ここでいう白猫・黒猫・灰色猫は、標準的な猫実験そのものではなく、A/B配分状態を読むための比喩である。

本実験で確認したことは次である。

```text
1. AB二体だけで、灰色猫固有相、灰色猫準安定相、自然選択相を分離できた。
2. C弱読出しでは、灰色猫準安定相を一方選択せずに読む窓が存在した。
3. D強観測では、灰色猫準安定相が白猫または黒猫へ選択された。
4. 灰色猫固有相では、D強観測を入れても灰色状態が保持された。
5. Dなし対照では白猫・黒猫への同等の落下は起きず、D起因選択として分離できた。
```

したがって、本実験の結果は次のようにまとめられる。

```text
灰色猫準安定相:
  弱く読むと灰色として読める。
  強く観測すると白猫または黒猫へ落ちる。

灰色猫固有相:
  弱く読んでも灰色として読める。
  強く観測しても灰色のまま保持される。
```

白猫に落ちた場合は、

```math
S \simeq +1,\qquad p_A\simeq 1,\qquad p_B\simeq 0
```

黒猫に落ちた場合は、

```math
S \simeq -1,\qquad p_A\simeq 0,\qquad p_B\simeq 1
```

である。

これは、白猫と黒猫の二匹へ分岐したというより、A/B配分が片側へほぼ全移動した状態である。

ただし、閉鎖量

```math
Q=p_A+p_B
```

は保存される。

---

## 1. 背景と目的

波の情報読出しシリーズでは、背景空間を先に置かず、複素位相波の閉鎖関係と干渉読出しから、位置様、運動量様、エネルギー様、加速度様に見える量を構成してきた。

前段の加速度基底・局在性交換実験では、交換干渉散乱行列によるAB相互作用のもとで、局在性と有効倍音次数が相互に再配分される条件を調べた。

本稿では、その流れを、A/B配分状態の選択問題へ移す。

問いは次である。

```text
白猫Aと黒猫Bが閉鎖系内で混在する状態は、
単なる未選択状態なのか。
それとも、灰色猫固有状態として保持される場合があるのか。

また、弱い読出しCと強い観測Dでは、
その状態はどう変わるのか。
```

本実験では、この問いを四段階に分ける。

```text
Stage 1:
  AB二体だけで、灰色猫固有相、灰色猫準安定相、自然選択相を探す。

Stage 2:
  C弱読出しで、灰色猫準安定相を壊さずに読めるか調べる。

Stage 3:
  D強観測で、灰色猫準安定相が白猫または黒猫へ落ちるか調べる。

Stage 4:
  D観測開始ステップとD利得を掃引し、選択境界を読む。
```

## 1.1 状態遷移図

本実験の三つの代表的な読みを、AB、ABC、ABCD の段階に沿って図示すると次のようになる。

![白猫・黒猫・灰色猫 AB-ABC-ABCD 状態遷移図](gray_cat_state_transition_figures_v1/gray_cat_ab_abc_abcd_three_scenarios_v1.png)

上段は、ABで作った白+黒の準安定混在をCで弱く読み、D強観測で白または黒へ選択する場合である。

中段は、ABの段階で灰色猫固有相になり、C弱読出しでもD強観測でも灰色として保持される場合である。

下段は、Cが強すぎるため、ABCの段階ですでに白または黒への選択が起こる場合である。

## 1.2 実測値による状態遷移図

同じ三つの読みを、実際に計算された `p_A`, `p_B`, `S` の時系列で図示すると次のようになる。

![白猫・黒猫・灰色猫 実測値 AB-ABC-ABCD 遷移](gray_cat_observed_value_transition_figures_v1/gray_cat_ab_abc_abcd_observed_values_three_scenarios_v1.png)

実線は内部状態の `p_A`, `p_B` である。

破線は、比較のために

```math
\frac{S+1}{2}
```

として同じ縦軸へ写した選択秩序変数である。

点線は、ABC区間ではC読出し、ABCD区間ではD読出しを示す。

この図では、準安定混在、灰色固有相、Cが強すぎる場合の三者が、同じ `p_A,p_B,S` の読出し上で分離している。

---

## 2. 変数と相分類

## 2.1 A/B配分

A/B配分は、

```math
p_A=|a|^2,\qquad p_B=|b|^2
```

で読む。

閉鎖量は、

```math
Q=p_A+p_B
```

である。

A/B選択秩序変数を、

```math
S=p_A-p_B
```

とする。

このとき、

```math
p_A=\frac{1+S}{2},
\qquad
p_B=\frac{1-S}{2}
```

である。

## 2.2 白猫・黒猫・灰色猫

本稿では、次の読替えを用いる。

| 状態 | 数値条件の目安 | 読み |
|---|---|---|
| 白猫 | `S≈+1` | `p_A≈1`, `p_B≈0` |
| 黒猫 | `S≈-1` | `p_A≈0`, `p_B≈1` |
| 灰色猫 | `S≈0` | `p_A≈0.5`, `p_B≈0.5` |

## 2.3 相分類

AB二体だけの発展では、次の相を分類する。

```text
gray_eigen:
  S_mean ≈ 0
  S_amp ≈ 0
  S_drift ≈ 0

gray_metastable:
  S_mean ≈ 0
  0 < S_amp < S_gray_limit
  S_drift ≈ 0

natural_selection:
  CもDもないのに S -> +1 または S -> -1 へ落ちる。

large_oscillation:
  Sが大きく振動し、A/B偏りが明確に現れる。
```

本稿の中心対象は `gray_metastable` と `gray_eigen` である。

---

## 3. モデル

## 3.1 AB交換相互作用

AB二体の最小交換写像は、

```math
\begin{pmatrix}
a_{k+1}\\
b_{k+1}
\end{pmatrix}
=
U_\epsilon
\begin{pmatrix}
a_k\\
b_k
\end{pmatrix}
```

```math
U_\epsilon
=
\begin{pmatrix}
\cos\epsilon & i\sin\epsilon\\
i\sin\epsilon & \cos\epsilon
\end{pmatrix}
```

である。

必要に応じて、閉鎖系内の安定性を調べるため、弱い復元項または弱い非線形項を加える。

ただし、Cなし・Dなしで自然選択する条件は、観測による選択の候補から除外する。

## 3.2 C弱読出し

Cは、A/B配分を弱く読む内部読出し波である。

可視度を、

```math
v_C=\frac{g_C}{g_C+\kappa_C}
```

とし、

```math
S_C=v_CS
```

```math
C_A=\frac{1+S_C}{2},
\qquad
C_B=\frac{1-S_C}{2}
```

として読む。

Cは、A/B配分を読むが、一方選択を起こさない窓を探すために用いる。

## 3.3 D強観測

Dは、現在の `S` を読んだ方向へ増幅する強観測写像である。

可視度を、

```math
v_D=\frac{g_D}{g_D+\kappa_D}
```

とし、

```math
S_D=v_DS
```

を読む。

Dバックアクションは、

```math
S_{k+1}
=
S_k
+
G_D S_D(1-S_k^2)
```

として実装する。

ここで `G_D` はD利得である。

同じD開始状態からDなし対照を並走させ、Dなし対照でも同じ白猫または黒猫選択が起こる場合は、D起因選択とは数えない。

---

## 4. Stage 1: AB準安定界面探索

## 4.1 条件

```text
C = off
D = off
steps = 4096
S_gray_limit = 0.05
selection_limit = 0.95
```

掃引した主なパラメータは次である。

```text
epsilon_values = (0.0, 1e-05, 3e-05, 0.0001, 0.0003, 0.001, 0.003, 0.01)
stability_gain_values = (-0.01, -0.002, 0.0, 0.002, 0.01)
```

## 4.2 結果

| phase | count |
|---|---:|
| `gray_eigen` | `1248` |
| `gray_metastable` | `733` |
| `large_oscillation` | `470` |
| `natural_selection` | `1070` |
| `unstable_or_drifting` | `2079` |

AB二体だけで、灰色猫固有相、灰色猫準安定相、大振幅振動領域、自然選択相が分離して観測された。

代表的な灰色猫準安定候補は次である。

| epsilon | phi/pi | s0 | gain | S_mean | S_amp | S_drift |
|---:|---:|---:|---:|---:|---:|---:|
| `0.01` | `0` | `0.01` | `0` | `0.000173315` | `0.02` | `0.00205793` |
| `0.01` | `1` | `0.01` | `0` | `0.000173315` | `0.02` | `0.00205793` |
| `0.003` | `0.0833333` | `0` | `-0.002` | `-0.00213507` | `0.020193` | `0.00227521` |

---

## 5. Stage 2: C弱読出し窓

## 5.1 条件

```text
C = on
D = off
steps = 4096
readout_kappa = 0.02
g_C_values = (0.0, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0)
backaction_scale_values = (0.0, 1e-05, 5e-05, 0.0001, 0.0005, 0.001, 0.002, 0.005, 0.01)
```

## 5.2 結果

```text
total_cases = 1080
C_window_count = 247
C_informative_window_count = 144
C_nonzero_backaction_window_count = 114
```

C後の相分類は次であった。

| phase_after_C | count |
|---|---:|
| `gray_eigen` | `90` |
| `gray_metastable` | `786` |
| `large_oscillation` | `141` |
| `natural_selection` | `2` |
| `unstable_or_drifting` | `61` |

灰色猫準安定候補にCを加え、A/B配分を読めるが一方選択を起こさない結合窓を確認した。

非ゼロCバックアクションを持つ代表的読出し窓は次である。

| epsilon | phi/pi | s0 | base_gain | g_C | c_gain | C_rel_err | C_bias_delta | phase_after_C |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `0.003` | `0.0833333` | `0` | `-0.002` | `1` | `5e-05` | `0.0196078` | `0.00122196` | `gray_metastable` |
| `0.003` | `0.0833333` | `0` | `-0.002` | `1` | `1e-05` | `0.0196078` | `0.000220996` | `gray_metastable` |
| `0.003` | `0.0833333` | `0.001` | `-0.002` | `1` | `1e-05` | `0.0196078` | `0.000222243` | `gray_metastable` |

---

## 6. Stage 3: D強観測応答

## 6.1 条件

```text
Stage 1: AB準安定界面探索済み
Stage 2: C読出し窓確認済み
Stage 3: D強観測応答
pre_steps_values = (0, 1, 2, 5, 10, 20, 50, 100, 250, 500, 1000, 2000)
d_steps = 2048
c_modes = (record_only, weak_C_window)
g_D_values = (0.0, 0.05, 0.1, 0.2, 0.5, 1.0)
d_backaction_scale_values = (0.0, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0)
```

## 6.2 全体結果

```text
total_cases = 7056
D_induced_selection_count = 2016
```

D結果分類は次である。

| D_outcome | count |
|---|---:|
| `white_selected` | `1244` |
| `black_selected` | `772` |
| `gray_kept_eigen` | `1062` |
| `unresolved` | `3978` |

Dなし対照では、次であった。

| baseline_outcome | count |
|---|---:|
| `gray_kept_eigen` | `1176` |
| `unresolved` | `5880` |

Dなし対照では、白猫または黒猫への同等の落下は検出されなかった。

## 6.3 灰色猫固有相

灰色猫固有相では、強D条件でも次が確認された。

| pre | C_mode | S_start | outcome | S_mean_after_D | S_amp_after_D | Q_err |
|---:|---|---:|---|---:|---:|---:|
| `0` | `record_only` | `0` | `gray_kept_eigen` | `0` | `0` | `2.22e-16` |
| `0` | `weak_C_window` | `0` | `gray_kept_eigen` | `0` | `0` | `2.22e-16` |
| `20` | `record_only` | `0` | `gray_kept_eigen` | `0` | `0` | `2.22e-16` |
| `20` | `weak_C_window` | `0` | `gray_kept_eigen` | `0` | `0` | `2.22e-16` |

灰色猫固有相は、D強観測によって白猫または黒猫へ落ちなかった。

## 6.4 灰色猫準安定相

灰色猫準安定相では、D強観測により白猫または黒猫への選択が起きた。

代表例は次である。

| epsilon | phi/pi | s0 | gain | pre | S_start | g_D | D_gain | outcome | baseline | S_mean_after_D |
|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|
| `0.003` | `0.0833333` | `0` | `-0.002` | `0` | `0` | `1` | `1` | `black_selected` | `unresolved` | `-1` |
| `0.003` | `0.0833333` | `0` | `-0.002` | `0` | `0` | `0.5` | `0.25` | `black_selected` | `unresolved` | `-0.999878` |
| `0.01` | `0` | `0.01` | `0` | `0` | `0.02` | `1` | `0.2` | `white_selected` | `unresolved` | `0.997467` |

## 6.5 大振幅分離領域

大振幅分離領域では、D結果がC読出し符号と一致する条件を確認した。

| pre | S_start | C_sign | outcome | S_mean_after_D | agreement |
|---:|---:|---|---|---:|---|
| `0` | `0.06` | `A` | `white_selected` | `1` | `same_sign` |
| `1` | `0.0556528` | `A` | `white_selected` | `1` | `same_sign` |
| `2` | `0.0513121` | `A` | `white_selected` | `1` | `same_sign` |

この領域では、Dが新たに選択を作ったというより、すでに系内に現れたA/B偏りを強く読んだものとして扱う。

---

## 7. Stage 4: D選択境界

## 7.1 条件

対象は、Stage 3でD起因選択が確認された灰色猫準安定候補である。

```text
target = gray_metastable candidates only
d_steps = 2048
pre_steps_values_count = 73
C_modes = (record_only, weak_C_window)
```

D利得は、境界近傍を細かく刻んだ。

```text
D_gain_values =
(0.0, 0.005, 0.01, 0.015, 0.02,
 0.0225, 0.025, 0.0275, 0.03,
 0.0325, 0.035, 0.0375, 0.04,
 0.045, 0.05, 0.055, 0.06,
 0.065, 0.07, 0.0725, 0.075,
 0.0775, 0.08, 0.0825, 0.085,
 0.0875, 0.09, 0.095, 0.1,
 0.12, 0.15, 0.2, 0.3,
 0.5, 0.75, 1.0)
```

## 7.2 結果

```text
total_rows = 15768
boundary_points = 438
selection_possible_boundary_points = 438
no_selection_boundary_points = 0
min_D_gain_overall = 0.0225
max_min_D_gain_overall = 0.065
```

候補別の境界は次である。

| case_id | boundary_points | selection_possible | no_selection | min_gain | max_min_gain | sign_counts |
|---|---:|---:|---:|---:|---:|---|
| `gray_metastable_0_eps0.01_phi0_s0.01_g0` | `146` | `146` | `0` | `0.065` | `0.065` | `A:90, B:56` |
| `gray_metastable_1_eps0.01_phi1_s0.01_g0` | `146` | `146` | `0` | `0.065` | `0.065` | `A:90, B:56` |
| `gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002` | `146` | `146` | `0` | `0.0225` | `0.0225` | `B:114, A:32` |

今回の掃引範囲では、対象とした全ての観測開始ステップでD起因選択が可能であった。

候補ごとの最小D利得は二段に分かれた。

```text
weak threshold candidate:
  min_D_gain = 0.0225

strong threshold candidates:
  min_D_gain = 0.065
```

---

## 8. 考察

本実験で重要なのは、灰色猫が一種類ではなかったことである。

灰色猫準安定相は、弱く読めば灰色として読める。

しかし、強いD観測を入れると、白猫または黒猫へ落ちる。

このとき、白猫になれば

```math
p_A\simeq 1,\qquad p_B\simeq 0
```

黒猫になれば

```math
p_A\simeq 0,\qquad p_B\simeq 1
```

である。

したがって、白猫と黒猫が二匹同時に出るのではない。

A/B配分が片側へ移り、反対側の成分は読出し上ほぼ消える。

一方、灰色猫固有相では、強いD観測を入れても

```math
p_A\simeq 0.5,\qquad p_B\simeq 0.5
```

が保持された。

これは、灰色猫準安定相と灰色猫固有相を、観測応答で区別できることを示す。

また、C弱読出しでは、選択を起こさずにA/B配分を読む窓があった。

したがって、本実験では、次の三つが分離された。

```text
AB二体だけで自然に落ちる自然選択。
C弱読出しで壊さず読める準安定灰色状態。
D強観測で白猫または黒猫へ落ちる選択状態。
```

---

## 9. 結論

本実験により、閉鎖系内のA/B配分状態について、次を確認した。

```text
1. AB交換相互作用だけで、灰色猫固有相、灰色猫準安定相、自然選択相を分類できる。
2. C弱読出しには、灰色猫準安定相を選択せずに読む窓がある。
3. D強観測では、灰色猫準安定相が白猫または黒猫へ落ちる。
4. 灰色猫固有相は、D強観測でも灰色状態を保持する。
5. D選択境界では、対象とした全ての準安定候補でD起因選択が可能だった。
6. 最小D利得は、候補によって 0.0225 と 0.065 の二段に分かれた。
```

したがって、白猫・黒猫・灰色猫という三状態は、単一の未選択状態ではなく、閉鎖系内の相分類と観測強度によって分かれる。

特に、

```text
灰色猫準安定相:
  弱読出しでは灰色。
  強観測では白または黒。

灰色猫固有相:
  強観測でも灰色。
```

という分岐を得た。

---

## 10. 出力

| 種別 | ファイル |
|---|---|
| 実験仕様 | `白猫黒猫_灰色猫準安定界面_AB-C-D段階実験仕様_v1.md` |
| Stage 1 スクリプト | `run_gray_cat_ab_metastable_interface_preliminary_v1.py` |
| Stage 1 JSON | `gray_cat_ab_metastable_interface_preliminary_result_v1/gray_cat_ab_metastable_interface_preliminary_result_v1.json` |
| Stage 1 CSV | `gray_cat_ab_metastable_interface_preliminary_result_v1/gray_cat_ab_metastable_interface_rows_v1.csv` |
| Stage 1 レポート | `gray_cat_ab_metastable_interface_preliminary_result_v1/gray_cat_ab_metastable_interface_report_v1.md` |
| Stage 2 スクリプト | `run_gray_cat_c_readout_window_preliminary_v1.py` |
| Stage 2 JSON | `gray_cat_c_readout_window_preliminary_result_v1/gray_cat_c_readout_window_preliminary_result_v1.json` |
| Stage 2 CSV | `gray_cat_c_readout_window_preliminary_result_v1/gray_cat_c_readout_window_rows_v1.csv` |
| Stage 2 レポート | `gray_cat_c_readout_window_preliminary_result_v1/gray_cat_c_readout_window_report_v1.md` |
| Stage 3 スクリプト | `run_gray_cat_d_observation_response_preliminary_v1.py` |
| Stage 3 JSON | `gray_cat_d_observation_response_preliminary_result_v1/gray_cat_d_observation_response_preliminary_result_v1.json` |
| Stage 3 CSV | `gray_cat_d_observation_response_preliminary_result_v1/gray_cat_d_observation_response_rows_v1.csv` |
| Stage 3 レポート | `gray_cat_d_observation_response_preliminary_result_v1/gray_cat_d_observation_response_report_v1.md` |
| Stage 4 スクリプト | `run_gray_cat_d_selection_boundary_preliminary_v1.py` |
| Stage 4 JSON | `gray_cat_d_selection_boundary_preliminary_result_v1/gray_cat_d_selection_boundary_preliminary_result_v1.json` |
| Stage 4 CSV | `gray_cat_d_selection_boundary_preliminary_result_v1/gray_cat_d_selection_boundary_rows_v1.csv` |
| Stage 4 境界表 | `gray_cat_d_selection_boundary_preliminary_result_v1/gray_cat_d_selection_boundary_table_v1.csv` |
| Stage 4 レポート | `gray_cat_d_selection_boundary_preliminary_result_v1/gray_cat_d_selection_boundary_report_v1.md` |
| 図化スクリプト | `draw_gray_cat_ab_abc_abcd_state_diagrams_v1.py` |
| AB-ABC-ABCD 状態遷移図 | `gray_cat_state_transition_figures_v1/gray_cat_ab_abc_abcd_three_scenarios_v1.png` |
| 準安定混在図 | `gray_cat_state_transition_figures_v1/gray_cat_metastable_mix_to_white_or_black_v1.png` |
| 灰色固有相図 | `gray_cat_state_transition_figures_v1/gray_cat_eigen_gray_to_gray_v1.png` |
| C強読出し選択図 | `gray_cat_state_transition_figures_v1/gray_cat_strong_c_selects_before_d_v1.png` |
| 実測値図化スクリプト | `draw_gray_cat_observed_value_transition_figures_v1.py` |
| 実測値 AB-ABC-ABCD 遷移図 | `gray_cat_observed_value_transition_figures_v1/gray_cat_ab_abc_abcd_observed_values_three_scenarios_v1.png` |
| 実測値時系列CSV | `gray_cat_observed_value_transition_figures_v1/gray_cat_ab_abc_abcd_observed_values_timeseries_v1.csv` |

---

# 参考文献

## 自己引用

1. 木原範昭「無名等振幅複合波モデル基本公理系 v4」Version DOI: `10.5281/zenodo.21316620`, Concept DOI: `10.5281/zenodo.21315735`, 2026.
2. 木原範昭「背景空間を仮定しない閉じた位相系におけるフェルミオン的二局所波の完全弾性反射の構成実験」Version DOI: `10.5281/zenodo.21332866`, Concept DOI: `10.5281/zenodo.21291018`, 2026.
3. 木原範昭「フェルミオン的逆相核による完全反射写像の干渉構成」Version DOI: `10.5281/zenodo.21332867`, Concept DOI: `10.5281/zenodo.21295479`, 2026.
4. 木原範昭「曲率付き閉鎖定常波による曲率繰り込みと完全反射安定性」Version DOI: `10.5281/zenodo.21332874`, Concept DOI: `10.5281/zenodo.21304039`, 2026.
5. 木原範昭「ABC閉鎖位相系における多ゲージ干渉読出し保存量の構成実験」Version DOI: `10.5281/zenodo.21332875`, Concept DOI: `10.5281/zenodo.21308049`, 2026.
6. 木原範昭「AB二体閉鎖位相系における調和読出しとc=1面積スイープ予備実験総括」Version DOI: `10.5281/zenodo.21332876`, Concept DOI: `10.5281/zenodo.21318696`, 2026.
7. 木原範昭「交換干渉散乱行列フェルミオン的衝突における加速度基底と局在性交換予備実験総括」Version DOI: `10.5281/zenodo.21333768`, Concept DOI: `10.5281/zenodo.21333766`, 2026.

## 外部参考文献

外部参考文献は、本稿の導出根拠ではなく、観測確率と猫型思考実験に関する標準的背景を示すために最小限に用いる。

8. Max Born, "Zur Quantenmechanik der Stossvorgaenge", *Zeitschrift fuer Physik* 37, 863-867, 1926. DOI: `10.1007/BF01397477`.
9. Erwin Schroedinger, "Die gegenwaertige Situation in der Quantenmechanik", *Naturwissenschaften* 23, 807-812, 823-828, 844-849, 1935.
