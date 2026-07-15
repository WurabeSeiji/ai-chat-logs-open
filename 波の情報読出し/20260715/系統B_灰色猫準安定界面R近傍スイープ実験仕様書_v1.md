# 系統B 灰色猫準安定界面R近傍スイープ実験仕様書 v1

**日付:** 2026-07-15  
**著者:** 木原 範昭  
**位置づけ:** 交換散乱係数R集中実験群・系統B仕様書  

---

## 1. 目的

本仕様書の目的は、20260714 の白猫黒猫・灰色猫準安定界面実験を `R` 依存性の観点から読み直し、灰色猫準安定相がどの反射率範囲で現れやすいかを調べることである。

系統Aが局在性と倍音移乗を見るのに対し、本系統では A/B 配分状態の相分類を見る。

本系統では、`R=0.70` 近傍が準安定界面で特別かどうかを確認する。

---

## 2. 直接の前提

直接の前提は次である。

```text
../20260714/白猫黒猫_灰色猫準安定界面_AB-C-D段階実験仕様_v1.md
../20260714/白猫黒猫灰色猫準安定界面におけるC弱読出しとD強観測選択予備実験総括 v1.md
../20260714/run_gray_cat_ab_metastable_interface_preliminary_v1.py
../20260714/gray_cat_ab_metastable_interface_preliminary_result_v1/gray_cat_ab_metastable_interface_report_v1.md
```

20260713 の局在性交換実験は、AB相互作用の前段参照として用いる。

---

## 3. 状態量

A/B配分を、複素振幅

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

を読む。

相分類は既存実験に合わせる。

```text
gray_eigen
gray_metastable
large_oscillation
natural_selection
unstable_or_drifting
```

---

## 4. R版AB交換写像

既存のAB交換相互作用を、`R` で制御される二チャネル交換写像として掃引する。

```math
\begin{pmatrix}
a_{k+1}\\
b_{k+1}
\end{pmatrix}
=
S_R
\begin{pmatrix}
a_k\\
b_k
\end{pmatrix}
```

```math
S_R
=
\begin{pmatrix}
r & t\\
t & r
\end{pmatrix}
```

```math
R=|r|^2,\qquad T=|t|^2,\qquad R+T=1
```

既存実験で使った弱い復元項または弱い非線形項は、同じ条件で保持する。

ここで見るのは、AB配分相が `R` によってどう変わるかである。

---

## 5. R掃引

主掃引は次である。

```text
0.680 <= R <= 0.710
Delta R = 0.001
```

対照として次を置く。

```text
R = 0.0
R = 0.5
R = 1.0
```

alpha 接続の共通プローブ点として、次を必ず含める。

```text
R_137 = 1 - sqrt(4 pi / 137.035999084)
R_128 = 1 - sqrt(4 pi / 128)
```

```text
R_137 ≈ 0.6971778791
R_128 ≈ 0.6866714657
```

これらは灰色猫準安定界面の判定値ではなく、系統A/Cと同じ位置を読むための固定評価点である。

有効境界が見えた場合は、近傍を細分化する。

```text
Delta R = 0.0005
Delta R = 0.00025
```

---

## 6. 掃引条件

既存のAB準安定界面探索条件を基本にする。

```text
C = off
D = off
steps = 4096
S_gray_limit = 0.05
selection_limit = 0.95
```

既存実験で用いた `s0`, `phi`, `stability_gain`, `noise` 条件を保持し、そこに `R` を加える。

---

## 7. 記録量

各 `R` と初期条件について次を記録する。

```text
R
T
condition_id
steps
S_mean
S_amp
S_drift
p_A_final
p_B_final
Q_error
phase
```

相分類ごとの件数を集約する。

```text
gray_eigen_count
gray_metastable_count
large_oscillation_count
natural_selection_count
unstable_or_drifting_count
```

---

## 8. R_star

本系統では、次を読む。

```text
R_star_gray_metastable = argmax_R gray_metastable_score
R_star_gray_eigen      = argmax_R gray_eigen_score
R_star_selection       = argmax_R selection_score
condition_at_R_star_gray_metastable
condition_at_R_star_gray_eigen
condition_at_R_star_selection
```

主判定は `R_star_gray_metastable` である。

`gray_metastable_score` は、灰色準安定相の件数、`S_amp` の小ささ、`S_drift` の小ささを合わせて読む。

自然選択相は、C/Dなしで白猫または黒猫へ落ちる相であり、観測選択の候補からは外す。

各 `R_star` について、次も記録する。

```text
distance_to_R_137
distance_to_R_128
```

---

## 9. C/D確認は後段に回す

本仕様書の第一実験では、C弱読出しとD強観測を入れない。

まず AB二体だけで、`R` によって灰色準安定界面がどう動くかを読む。

有効な `R_star_gray_metastable` が得られた場合に限り、第二実験として次を確認する。

```text
C弱読出しで灰色準安定相を壊さず読めるか。
D強観測で白猫または黒猫へ選択されるか。
```

---

## 10. R_band

各相分類スコアについて、最良値から次の範囲を記録する。

```text
within_5_percent
within_10_percent
```

`R_band_width` を記録し、系統Aの `R_band_width` と比較する。

---

## 11. 出力

出力ディレクトリは次とする。

```text
system_B_gray_cat_metastable_R_sweep_result_v1
```

出力ファイルは次を基本とする。

```text
system_B_gray_cat_metastable_R_sweep_rows_v1.csv
system_B_gray_cat_metastable_R_sweep_summary_v1.csv
system_B_gray_cat_metastable_R_sweep_result_v1.json
system_B_gray_cat_metastable_R_sweep_report_v1.md
system_B_gray_cat_metastable_R_phase_counts_v1.png
system_B_gray_cat_metastable_R_score_v1.png
```

---

## 12. 読み方

本系統で `R_star_gray_metastable` が系統Aの `R_star` と近い場合、局在性交換とA/B配分準安定界面が同じ交換散乱係数に敏感である可能性がある。

一方、系統Bの `R_star` が大きく異なる場合、灰色猫準安定界面は局在性交換とは異なる相図条件で決まる。

どちらの場合も、結果は総括論文で、`R` がどの幾何量に効いているかを分析する材料になる。
