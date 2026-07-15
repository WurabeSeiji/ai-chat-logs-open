# 系統A 局在性交換R近傍斉一スイープ実験仕様書 v1

**日付:** 2026-07-15  
**著者:** 木原 範昭  
**位置づけ:** 交換散乱係数R集中実験群・系統A仕様書  

---

## 1. 目的

本仕様書の目的は、20260713 の交換干渉散乱行列フェルミオン的衝突実験を精密化し、片側高次倍音条件で局在性と倍音分布が相手側チャネルへどの `R` で移乗しやすいかを読むことである。

本系統は、三系統の中で `R` に最も敏感であることが期待される。

ただし、`R=0.70` 近傍に有効点があるとは仮定しない。

---

## 2. 直接の前提

直接の前提は次である。

```text
../20260713/run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1.py
../20260713/交換干渉散乱行列フェルミオン的衝突における低局在性・倍音移乗読出し実験仕様書 v1.md
../20260713/交換干渉散乱行列フェルミオン的衝突における低局在性・倍音移乗予備実験総括 v1.md
../20260713/交換干渉散乱行列フェルミオン的衝突における加速度基底・低奇数倍音底・片側高次倍音予備実験検証メモ_v1.md
```

20260710 のフェルミオン的逆相核反射論文は、`R,T` 読出しの基礎参照とする。

20260711 の加速度様読出しは、本系統では直接使わない。

---

## 3. 基本写像

20260713 の二チャネル散乱行列をそのまま用いる。

```math
\begin{pmatrix}
A_{k+1}\\
B_{k+1}
\end{pmatrix}
=
S_R
\begin{pmatrix}
A_k\\
B_k
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

実装では、各衝突回で出射チャネルを再正規化し、次の衝突回の入力に戻す。

---

## 4. R掃引

本系統では、最初に広域粗スイープを行い、複数谷が存在しないかを確認する。

広域粗スイープは次である。

```text
0.600 <= R <= 0.900
Delta R = 0.010
```

この段階では、`R=0.70` 近傍が特別であるとは仮定しない。

広域粗スイープで谷候補を抽出した後、局所精密スイープを行う。

```text
0.680 <= R <= 0.710
Delta R = 0.001
```

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

これらは最小値を先に仮定するための点ではない。`R_star` がどちらに近いかを後で読むための固定評価点である。

必要に応じて、最良点近傍を細分化する。

```text
Delta R = 0.0005
Delta R = 0.00025
```

したがって、判定順は次である。

```text
Stage A0:
  0.600 から 0.900 の粗スイープで地形と複数谷を確認する。

Stage A1:
  粗スイープで得た谷候補の近傍を細かく掃引する。

Stage A2:
  R_137, R_128 との距離を診断量として読む。
```

---

## 5. N条件

主条件は、A側を基底波 `N_A=1` に固定し、B側の次数 `N_B` を掃引する片側次数条件である。

```text
(N_A, N_B) = (1, 1)
(N_A, N_B) = (1, 2)
(N_A, N_B) = (1, 3)
(N_A, N_B) = (1, 5)
(N_A, N_B) = (1, 10)
(N_A, N_B) = (1, 15)
(N_A, N_B) = (1, 31)
(N_A, N_B) = (1, 63)
```

このうち `(1,1)` は基底対照である。

追加対照として、高次側の同次数条件を置く。

```text
(N_A, N_B) = (31, 31)
(N_A, N_B) = (63, 63)
```

ここで `A` は低局在性の基底波、`B` は次数を変える対象波である。

低次 `N=2,3,5,10` を含めることで、`R` の谷が高次倍音で初めて現れるのか、低次からすでに現れるのかを確認する。

---

## 6. 倍音組立モード

単一次数スイープに加えて、B側の波形を複数倍音パケットとして組み立てる対照実験を行う。

目的は、`R` の谷が単一の次数差だけで決まるのか、偶数・奇数を含む倍音列の構造に依存するのかを確認することである。

### 6.1 single-N mode

A側は基底波、B側は単一次数波とする。

```text
A: [1]
B: [1], [2], [3], [5], [10], [15], [31], [63]
```

### 6.2 even packet mode

B側を基底波と偶数倍音のみで構成する。

```text
A: [1]
B: [1, 2, 4, 6]
```

### 6.3 alternating packet mode 1

B側を、基底、偶数、奇数、偶数、奇数の順で構成する。

```text
A: [1]
B: [1, 2, 3, 4, 5]
```

### 6.4 alternating packet mode 2

B側を、基底、奇数、偶数、奇数、偶数の順で構成する。

```text
A: [1]
B: [1, 3, 4, 5, 6]
```

各パケットは総ノルムを正規化して比較する。

初期位相と振幅規則は全パケットで共通にする。

---

## 7. 衝突回数

最大衝突回数は次とする。

```text
max_collision = 256
```

判定は全衝突回で行う。

図化する場合は、次を代表点として使う。

```text
collision = 0, 1, 2, 3, 5, 10, 20, 42, 64, 128, 256
```

---

## 8. 記録量

各 `R`, `N` 条件、衝突回ごとに次を記録する。

```text
R
T
Delta_F
N_A
N_B
collision
channel
mode
harmonic_packet_A
harmonic_packet_B
L
N_eff
H(n)
origin_A
origin_B
P_m_A
P_m_B
chi_center_cell_mass
```

集約量として次を記録する。

```text
L_gap_min
L_gap_min_collision
N_eff_gap_min
N_eff_gap_min_collision
joint_score_min
joint_score_min_collision
max_B_harmonic_transfer_to_A
max_A_harmonic_transfer_to_B
harmonic_transfer_collision
```

---

## 9. 判定量

基本量は次である。

```text
L_gap = abs(L_A - L_B)
N_eff_gap = abs(N_eff_A - N_eff_B)
```

倍音移乗は、初期倍音分布との類似度で読む。

```text
sim_A_to_A(k) = similarity(H_A(k), H_A(0))
sim_B_to_A(k) = similarity(H_A(k), H_B(0))
sim_A_to_B(k) = similarity(H_B(k), H_A(0))
sim_B_to_B(k) = similarity(H_B(k), H_B(0))
```

片側高次倍音条件での主読出しは次である。

```text
B_to_A_transfer(k) = sim_B_to_A(k)
A_to_B_transfer(k) = sim_A_to_B(k)
```

---

## 10. R_star

本系統では、次を読む。

```text
R_star_L        = argmin_R L_gap_min
R_star_N        = argmin_R N_eff_gap_min
R_star_transfer = argmax_R max_k B_to_A_transfer(k)
R_star_joint    = argmin_R joint_score_min
collision_at_R_star_L
collision_at_R_star_N
collision_at_R_star_transfer
collision_at_R_star_joint
```

`R_star_joint` は、正規化した `L_gap`, `N_eff_gap`, `1 - B_to_A_transfer` の合成スコアから読む。

ただし、各 `R_star` が大きく分離する場合は、単一の有効 `R` を定義しない。

各 `R_star` について、次も記録する。

```text
distance_to_R_137
distance_to_R_128
```

---

## 11. R_band

各 `N` 条件について、`joint_score_min` の最小値から次の範囲を記録する。

```text
within_5_percent
within_10_percent
```

それぞれの幅を、

```text
R_band_width_5
R_band_width_10
```

として記録する。

これにより、`N` が高くなると許容される `R` 範囲が狭くなるかを確認する。

---

## 12. 出力

出力ディレクトリは次とする。

```text
system_A_localization_exchange_R_sweep_result_v1
```

出力ファイルは、実行条件から作る `<stem>` を付けて保存する。

```text
<stem>_rows_v1.csv
<stem>_summary_v1.csv
<stem>_best_v1.csv
<stem>_collision_terrain_v1.csv
<stem>_result_v1.json
<stem>_report_v1.md
<stem>_scores_v1.png
<stem>_gap_depth_distribution_overview_v1.png
<stem>_gap_depth_distribution_deep_v1.png
```

`<stem>` には、倍音条件、`R` 条件、衝突回数を含める。

深さの主判定には、横軸を `R`、縦軸を `depth=-log10(gap_terrain_score)` とする分布図を用いる。

分布図では、各 `R` に対する全衝突回の深さ分布と、各 `R` の最深包絡線を読む。

完全同一対照では数値的に深さが過大になるため、図では表示用の上限を置く。生の深さはCSVに保存する。

これにより、単一の `R_star` だけでなく、反射係数方向に現れる谷の幅と、衝突回数方向の分布を確認する。

---

## 13. 分割実行

全量実行に時間がかかる場合は、同じスクリプトを条件別に分割して実行する。

対象スクリプトは次である。

```text
run_system_A_localization_exchange_R_sweep_preliminary_v1.py
```

主な分割引数は次である。

```text
--pairs        実行する (N_A,N_B) 条件。例: 1:2,1:3,1:63
--pair-start   既定 pairs の開始 index
--pair-stop    既定 pairs の終了 index
--r-min        実行する R の下限
--r-max        実行する R の上限
--r-values     明示的に実行する R 値
--max-collision 衝突回数の上限
--run-id       既定結果ディレクトリ下の分割出力名
--file-stem    出力ファイル名の明示stem
--output-dir   出力先ディレクトリ
--no-plots     CSV/JSON/レポートのみ作成し、図を省略する
```

実行例は次である。

```text
python3 run_system_A_localization_exchange_R_sweep_preliminary_v1.py --run-id lowN --pairs 1:1,1:2,1:3 --r-min 0.60 --r-max 0.90
python3 run_system_A_localization_exchange_R_sweep_preliminary_v1.py --run-id highN --pairs 1:15,1:31,1:63 --r-min 0.68 --r-max 0.71
python3 run_system_A_localization_exchange_R_sweep_preliminary_v1.py --run-id probe --pairs 1:2 --r-values 0.697177879128,0.686671465671 --max-collision 64 --no-plots
python3 run_system_A_localization_exchange_R_sweep_preliminary_v1.py --file-stem odd_kernel_N_1_2_3_5_15_63 --pairs 1:1,1:2,1:3,1:5,1:15,1:63
```

分割実行では、出力ファイル名に、倍音条件、R範囲、衝突回数から作ったstemを付ける。

`--file-stem` を指定した場合は、そのstemを用いる。

`--run-id` を指定した場合は、出力ディレクトリとファイルstemの両方に反映する。

総括時には、分割出力の `summary` と `best` を統合して比較する。

---

## 14. 読み方

本系統で `R_star` が鋭く出る場合、局在性と倍音分布の交換に対して `R` が敏感であると読む。

同次数対照でも同じ構造が出る場合は、片側高次倍音由来の移乗ではなく、単なるチャネル混合として扱う。

`R=0.5` で良い結果が出る場合は、一回平均化に近い単純混合と区別する。

本系統の結果は、三系統比較の第一基準になる。
