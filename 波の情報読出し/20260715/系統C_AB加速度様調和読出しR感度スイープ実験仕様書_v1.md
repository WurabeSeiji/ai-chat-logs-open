# 系統C AB加速度様調和読出しR感度スイープ実験仕様書 v1

**日付:** 2026-07-15  
**著者:** 木原 範昭  
**位置づけ:** 交換散乱係数R集中実験群・系統C仕様書  

---

## 1. 目的

本仕様書の目的は、20260711 のAB二体加速度様調和読出しに部分反射 `R` を導入し、調和読出しが `R` に敏感か鈍感かを調べることである。

本系統では、`R=0.70` 近傍に有効点が出ることを期待しすぎない。

重要なのは、加速度様調和読出しが `R` に鋭く反応するのか、それとも広い `R` 範囲で安定するのかを分類することである。

系統Aでは、片側倍音条件において `N=1` が対照となり、`N>=2` では奇数・偶数・混在パケットに依らず同じ `R` 近傍へ評価点が集中した。

系統Bでは、灰色猫準安定界面において、最良点が一つに限られず、`R_137` 近傍、`alpha~128` 名目点近傍、その他の局所ピークを分けて読む必要が出た。

系統Cでは、この知見を受けて、加速度様調和読出しが同じ `R` 地形を持つのか、それとも `R` に鈍感な保存読出しとして振る舞うのかを調べる。

---

## 2. 直接の前提

直接の前提は次である。

```text
../20260711/AB二体閉鎖位相系における調和読出しとc=1面積スイープ予備実験総括 v4.md
../20260711/run_ab_two_body_fermionic_reflection_harmonic_readout_v4.py
../20260711/ab_two_body_fermionic_reflection_harmonic_readout_preliminary_result_v4/ab_two_body_fermionic_reflection_harmonic_readout_report_v4.md
run_system_C_ab_acceleration_harmonic_R_sensitivity_preliminary_v1.py
```

20260713 の散乱行列写像は、部分反射 `R` の導入に用いる。

---

## 3. 基本方針

既存のAB二体加速度様調和読出しでは、完全反射極限として、

```text
R = 1
T ≈ 0
q_out_factor_diagnostic = -1
```

が診断量として読まれる。

本実験では、これを固定せず、`R` を掃引する。

部分反射の場合、透過率と反射率から読まれる診断量は、

```text
q_out_factor_diagnostic = T - R = 1 - 2R
```

として読む。

ただし、実装ではこの量を `chi_read` に掛けない。

A/B二チャネル散乱行列を入射チャネルへ作用させ、出射チャネル差から `chi_read` と `eta_read` を読む。

ただし、主目的は、局在性交換や倍音移乗を読むことではない。

主目的は、次である。

```text
AB二体の調和読出しが、部分反射Rにより壊れるか。
chi-tau 面積読出しが、部分反射Rにより壊れるか。
R=0.70近傍だけで特別に安定するか。
または、広いR範囲で安定するか。
```

---

## 4. R掃引

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

対照として次を置く。

```text
R = 0.0
R = 0.5
R = 1.0
```

alpha 接続の共通プローブ点として、次を必ず含める。

```text
R_137 = 1 - sqrt(4 pi / 137.035999084)
R_128_nominal = 1 - sqrt(4 pi / 128)
```

```text
R_137 ≈ 0.6971778791
R_128_nominal ≈ 0.6866714657
```

`R_137` は低エネルギー側の精密な `1/alpha` を用いた固定評価点である。

一方、`R_128_nominal` は高エネルギー領域で `1/alpha` が 128 近傍へ走ることを読むための名目プローブであり、精密な測定値そのものではない。

したがって、`R_128_nominal` については、点そのものだけでなく、その近傍に局所ピークが立つかも重視する。

これらは加速度様読出しの安定点を仮定するものではなく、系統A/Bと同じ位置で感度を確認するための固定評価点である。

本系統では、単一の最良Rだけを採用しない。

第一段階では次を同時に記録する。

```text
global best peak
local peak candidates
R_137 probe point
alpha~128 nominal probe point
```

さらに、第一段階の広域粗スイープで得られた上位ローカルピークを、そのまま結論にしない。

実験Bと同様に複数ピークが立つ可能性があるため、次の二段階で実行する。

```text
Stage 1:
  R = 0.60 から 0.90 を広域粗スイープする。
  R_137 と R_128_nominal も固定プローブとして含める。

Stage 2:
  Stage 1 で得られた上位ローカルピークそれぞれの近傍を局所精密化する。
  α近傍だけを精密化するのではなく、広域で立った複数ピークを保持する。
```

この処理により、`R≈0.70` 近傍だけでなく、別の谷が出た場合も見落とさない。

---

## 5. 実験条件

既存のV4条件をできるだけ保持する。

```text
harmonic cases: 既存48条件
c1 cases: 既存144条件
external_c_used_any = false
f_A_or_f_B_used_any = false
```

変更するのは、フェルミオン型反跳写像に入る `R` のみである。

系統A/Bで用いた `N` 条件および奇数・偶数・混在パケットは、第一実装の主変数にはしない。

理由は、20260711 のAB二体加速度様調和読出しが、局在波形や倍音パケットの移乗ではなく、調和読出しと `chi-tau` 面積読出しの成立性を読む実験だからである。

したがって、系統Cの第一実験では、既存V4の条件群を保持したまま `R` だけを掃引する。

後続実装で、調和読出し側に明示的な倍音カーネル軸を導入する場合のみ、N系列と奇偶パケットを追加対照として扱う。

---

## 6. 記録量

各 `R` について、既存レポートと同等の量を記録する。

```text
sweep_region
R
T
Delta_F
q_out_factor_diagnostic
q_out_factor_applied
full_two_channel_scattering_used_rate
max_scattering_unitarity_error
max_Q_closed_abs
harmonic_case_count
c1_case_count
fermionic_regular_cell_harmonic_consistent_all_cases
fermionic_regular_cell_harmonic_consistent_nonstrong_modes
fermionic_strong_readout_perturbs_harmonic_projection
fermionic_max_f_AB_projection_error_regular_nonstrong
fermionic_max_f_AB_projection_error_regular_strong
readout_off_decay_max_abs
readout_strong_decay_min_abs
fermionic_c1_area_sweep_detected_all_cases
score_C
classification_C
peak_rank
peak_kind
distance_to_R_137
distance_to_R_128_nominal
```

---

## 7. 判定量

本系統では、次を主判定とする。

```text
harmonic_valid_rate
c1_area_valid_rate
projection_error_nonstrong
projection_error_strong
```

総合スコアは、調和読出しと `c=1` 面積読出しが壊れないことを確認するための管理用集約値である。

```text
score_C =
  harmonic_valid_penalty
  + c1_area_penalty
  + normalized_projection_error
```

ただし、`score_C` は異なる種類の判定量を足した管理用指標であり、物理的なR地形そのものとしては読まない。

R地形は、少なくとも次を分けて読む。

```text
projection_penalty:
  加速度様調和読出しの射影誤差地形。
```

---

## 8. R_starとR感度

本系統では、次を読む。

```text
R_star_C = argmin_R score_C
condition_at_R_star_C
distance_to_R_137
distance_to_R_128_nominal
```

ただし、`R_star_C` は管理用集約値の代表点であり、唯一の結論点ではない。

系統Cでは、管理用の局所点として次も同時に記録する。

```text
peak_rank
peak_kind
R_peak
score_C
condition_at_peak
distance_to_R_137
distance_to_R_128_nominal
```

`peak_kind` は次を含む。

```text
local_min
probe_R137
probe_alpha128_nominal
local_min+probe_R137
local_min+probe_alpha128_nominal
```

主判定は `R_star_C` そのものではなく、次の分類である。

```text
sharp:
  R=0.70 近傍でのみ調和読出しが安定する。

broad:
  広い R 範囲で調和読出しが安定する。

flat:
  掃引範囲内でほぼ差がない。

multi_peak:
  複数の局所安定点が立つ。

outside:
  共通局所窓の外側に主要安定点が立つ。

broken:
  部分反射 R を入れると調和読出しが成立しない。
```

`broad` または `flat` は失敗ではない。

その場合、加速度様調和読出しは局在性交換ほど `R` に敏感ではないと読む。

---

## 9. R地形図

系統A/Bと比較しやすくするため、色図ではなく、Rに沿った地形図を用いる。

ただし、系統Cでは、条件別の最良点を結ぶ包絡線は主図に使わない。

理由は、`c=1` 面積読出しに、広い `R` 範囲でほぼ完全一致する条件が含まれるためである。

その条件を包絡線に入れると、図全体が上限に張り付き、加速度様読出しのR感度が読めなくなる。

したがって、系統Cの主図では、合成スコアをR地形として描かない。

主図は次の単独図とする。

```text
x axis: R
y axis: projection_penalty
```

`projection_penalty` は小さいほどよい。

この図では、加速度様調和読出しの射影誤差が `R` に対してどのような地形を持つかだけを読む。

`score_C` は出力表とピーク管理には残すが、この図の縦軸には入れない。

```text
projection_penalty
```

粗い `R` 刻みだけでは、狭い局所地形を見落とす可能性がある。

そのため、主要候補の近傍では、必要に応じて追加の高分解能掃引を行う。

例:

```text
0.748 <= R <= 0.752
Delta R = 0.0001
```

この追加掃引では、`projection_penalty` の局所地形を確認する。

さらに、ピーク近傍だけを拡大した図を出す。

この図では、上位局所安定点、`R_137` 固定プローブ、`alpha~128` 名目プローブ、`alpha~128` 近傍の局所安定点を、`R` と管理用 `score_C` の数値付きで表示する。

---

## 10. 出力

出力ディレクトリは次とする。

```text
system_C_ab_acceleration_harmonic_R_sensitivity_result_v1
```

出力ファイルは次を基本とする。

```text
<file_stem>_rows_v1.csv
<file_stem>_condition_rows_v1.csv
<file_stem>_summary_v1.csv
<file_stem>_best_v1.csv
<file_stem>_peaks_v1.csv
<file_stem>_result_v1.json
<file_stem>_report_v1.md
<file_stem>_scores_v1.png
<file_stem>_stability_depth_distribution_overview_v1.png
<file_stem>_stability_depth_distribution_deep_v1.png
<file_stem>_peak_zoom_v1.png
```

`file_stem` には、掃引条件を含め、再実行時の上書きを避ける。

図の構成は、スコア概観図、射影誤差R地形図、ピーク拡大図を基本とする。

`stability_depth_distribution` というファイル名は既存ワークフローとの互換のため保持するが、系統Cでは射影誤差R地形図を意味する。

---

## 11. 読み方

本系統で `R=0.70` 近傍に鋭い安定点または独立した局所安定点が出る場合、交換散乱係数は加速度様調和読出しにも直接関係している可能性がある。

一方、広い `R` 範囲で調和読出しが保たれる場合、`R=0.70` 近傍は加速度様読出し一般ではなく、局在性、倍音移乗、状態選択に関係する係数として読む。

この区別が、三系統比較の重要な対照になる。
