# 系統B 灰色猫準安定界面R近傍スイープ実験仕様書 v1

**日付:** 2026-07-15  
**著者:** 木原 範昭  
**位置づけ:** 交換散乱係数R集中実験群・系統B仕様書  

---

## 1. 目的

本仕様書の目的は、20260714 の白猫黒猫・灰色猫準安定界面実験を `R` 依存性の観点から読み直し、灰色猫準安定相がどの反射率範囲で現れやすいかを調べることである。

系統Aでは、片側倍音条件において、`N=1` は対照となり、`N>=2` では奇数・偶数・混在パケットに依らず同じ `R` 近傍へ評価点が集中した。

系統Bでは、この知見を受けて、灰色猫準安定界面でも次を確認する。

```text
N=1 対照では R 地形が鈍感か。
N>=2 条件では灰色準安定相の R 地形が現れるか。
奇数・偶数・混在パケットで R 地形が変わるか。
```

本系統の主判定は、局在性ではなく A/B 配分状態の相分類である。

ただし、20260714 の灰色猫実験は `S=p_A-p_B` を読む二成分準安定モデルであり、倍音波形を直接展開するモデルではない。

したがって本仕様書では、N条件と奇偶パケットを、系統Aと比較するための外部比較軸として記録する。Nまたは奇偶による差が出ない場合、それは本系統が倍音構成ではなくA/B配分力学に支配されていることを示す。

---

## 2. 直接の前提

直接の前提は次である。

```text
../20260714/白猫黒猫_灰色猫準安定界面_AB-C-D段階実験仕様_v1.md
../20260714/白猫黒猫灰色猫準安定界面におけるC弱読出しとD強観測選択予備実験総括 v1.md
../20260714/run_gray_cat_ab_metastable_interface_preliminary_v1.py
../20260714/gray_cat_ab_metastable_interface_preliminary_result_v1/gray_cat_ab_metastable_interface_report_v1.md
../20260715/run_system_A_localization_exchange_R_sweep_preliminary_v1.py
```

20260714 のAB二体準安定界面実験は、灰色猫相分類の母体である。

20260715 系統Aは、`R` 掃引、N系列、奇偶パケット、出力図構成の母体である。

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

既存の `epsilon` 型AB交換を、`R` で制御される二チャネル散乱写像へ置き換える。

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

実装では、系統Aと同じく、各交換後にABノルムを再正規化し、次の交換回の入力に戻す。

20260714 で用いた弱い復元項または弱い非線形項は、同じ条件で保持する。

---

## 5. R掃引

最初に広域粗スイープを行い、複数ピークや平坦領域を確認する。

```text
0.600 <= R <= 0.900
Delta R = 0.010
```

その後、局所精密スイープを行う。

```text
0.680 <= R <= 0.710
Delta R = 0.001
```

この局所範囲は、`R_137` 近傍を含む共通確認窓である。

ただし、最終的な精密化はこの窓だけに固定しない。

広域粗スイープで複数ピークが立つ場合は、各ピーク近傍を別々に精密化候補とする。

対照として、次を別に記録する。

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

`R_137` は低エネルギー側の精密な `1/alpha` を用いた固定評価点である。

一方、`R_128` は高エネルギー領域で `1/alpha` が 128 近傍へ走ることを読むための名目プローブであり、精密な測定値そのものではない。

したがって、`R_128` については、点そのものだけでなく、その近傍に局所ピークが立つかも重視する。

これらは判定値ではなく、系統A/Cと同じ位置を読むための固定評価点である。

系統Bでは、単一の最良Rだけを採用しない。

灰色準安定界面は多峰的なR地形を持つ可能性があるため、広域スイープで立った局所ピークを順位付きで保存する。

最深ピークが `R_137` 近傍でなくても、`R_137` 近傍に独立したピークが立つなら、それ自体を重要な候補として扱う。

したがって、第一段階では次を同時に記録する。

```text
global best peak
local peak candidates
R_137 probe point
alpha~128 nominal probe point
```

---

## 6. N条件

主条件は、A側を基底波相当 `N_A=1` に固定し、B側の次数を変える片側次数条件である。

```text
(N_A, N_B) = (1, 1)
(N_A, N_B) = (1, 2)
(N_A, N_B) = (1, 3)
(N_A, N_B) = (1, 5)
(N_A, N_B) = (1, 15)
(N_A, N_B) = (1, 63)
```

`(1,1)` は基底対照である。

系統Bでは、`N` はA/B配分の名前ではない。

相分類は常に `S=p_A-p_B` から読む。

この二成分モデルでは、N条件は `S` 写像へ直接入らない。N条件は、系統Aと同じ実験表を作るための比較軸であり、N依存性が出ないことも有効な結果として扱う。

---

## 7. 倍音組立モード

単一次数スイープに加えて、B側の波形を複数倍音パケットとして組み立てる対照実験を行う。

### 7.1 single-N mode

```text
A: [1]
B: [1], [2], [3], [5], [15], [63]
```

### 7.2 even packet mode

```text
A: [1]
B: [1, 2, 4, 6]
```

### 7.3 alternating packet mode 1

```text
A: [1]
B: [1, 2, 3, 4, 5]
```

### 7.4 alternating packet mode 2

```text
A: [1]
B: [1, 3, 4, 5, 6]
```

各パケットは総ノルムを正規化して比較する。

ただし、本系統の第一実装ではパケット構成は相分類表の比較軸であり、`S` 更新式そのものは同じである。

---

## 8. 掃引条件

既存のAB準安定界面探索条件を基本にする。

```text
C = off
D = off
steps = 4096
S_gray_limit = 0.05
selection_limit = 0.95
```

既存実験で用いた `s0`, `phi`, `stability_gain`, `noise` 条件を保持し、そこに `R` と倍音構成を加える。

---

## 9. 記録量

各 `R`, `N`, 倍音構成、初期条件について次を記録する。

```text
case_id
mode
harmonic_packet_A
harmonic_packet_B
N_A
N_B
R
T
condition_id
step
S_mean
S_amp
S_drift
S_max_abs
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

## 10. R_star

本系統では、次を読む。

```text
R_star_gray_metastable = argmax_R gray_metastable_score
R_star_gray_eigen      = argmax_R gray_eigen_score
R_star_selection       = argmax_R selection_score
R_star_joint           = argmax_R joint_gray_score
```

主判定は `R_star_joint` とする。

ただし、`R_star_joint` は地形の代表点であり、唯一の結論点ではない。

系統Bでは、R地形そのものを読むため、次も同時に記録する。

```text
peak_rank
peak_kind
R_peak
joint_gray_score
best_phase
best_gray_depth
best_S_amp
distance_to_R_137
distance_to_R_128
```

`peak_kind` は次を含む。

```text
local_max
probe_R137
probe_alpha128_nominal
local_max+probe_R137
local_max+probe_alpha128_nominal
```

これにより、最深ピーク、低エネルギー alpha 近傍ピーク、高エネルギー alpha 近傍ピーク、その他の局所ピークを分けて読む。

`joint_gray_score` は、灰色準安定相の件数、`S_amp` が準安定候補振幅に近いこと、`S_drift` の小ささ、自然選択相の少なさを合わせて読む。

完全な `S=0` 固有灰色相は、灰色猫としては重要だが、準安定界面のR地形を読む主判定からは分ける。

第一実装では、20260714 の準安定候補に合わせて、目標振幅を次に置く。

```text
S_amp_target = 0.02
```

自然選択相は、C/Dなしで白猫または黒猫へ落ちる相であり、観測選択の候補からは外す。

各 `R_star` について、次も記録する。

```text
distance_to_R_137
distance_to_R_128
```

---

## 11. R地形図

系統Aと比較しやすくするため、色図ではなく、分布図を用いる。

各 `R` に対し、全条件の相分類スコアを縦方向に散布する。

```text
x axis: R
y axis: gray depth
```

`gray depth` は、灰色準安定相としての良さを対数的に読める量とする。

```text
gray_error = |S_mean| + |S_amp - S_amp_target| + S_drift + phase_penalty
gray_depth = -log10(gray_error)
```

灰色点は全条件、黒線は各 `R` における最良条件の包絡線とする。

スコア図では、上位局所ピークも点で示す。

これにより、最大値を一つだけ読むのではなく、R地形上に複数ピークが立つかどうかを確認する。

さらに、ピーク近傍だけを拡大した図を出す。

この図では、上位3ピーク、`R_137` 固定プローブ、`alpha~128` 名目プローブ、`alpha~128` 近傍の局所ピークを、`R` と `joint_gray_score` の数値付きで表示する。

系統Aと同じく、次の二種類の図を出す。

```text
overview distribution
deep distribution
```

---

## 12. R_band

各相分類スコアについて、最良値から次の範囲を記録する。

```text
within_5_percent
within_10_percent
```

`R_band_width` を記録し、系統Aの `R_band_width` と比較する。

---

## 13. C/D確認は後段に回す

本仕様書の第一実験では、C弱読出しとD強観測を入れない。

まず AB二体だけで、`R` によって灰色準安定界面がどう動くかを読む。

有効な `R_star_gray_metastable` が得られた場合に限り、第二実験として次を確認する。

```text
C弱読出しで灰色準安定相を壊さず読めるか。
D強観測で白猫または黒猫へ選択されるか。
```

---

## 14. 出力

出力ディレクトリは次とする。

```text
system_B_gray_cat_metastable_R_sweep_result_v1
```

出力ファイルは、系統Aと同じく `file_stem` を付けて上書きを避ける。

```text
<file_stem>_rows_v1.csv
<file_stem>_samples_v1.csv
<file_stem>_summary_v1.csv
<file_stem>_best_v1.csv
<file_stem>_peaks_v1.csv
<file_stem>_result_v1.json
<file_stem>_report_v1.md
<file_stem>_scores_v1.png
<file_stem>_gray_depth_distribution_overview_v1.png
<file_stem>_gray_depth_distribution_deep_v1.png
<file_stem>_peak_zoom_v1.png
```
