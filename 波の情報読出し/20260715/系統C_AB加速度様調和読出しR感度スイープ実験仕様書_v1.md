# 系統C AB加速度様調和読出しR感度スイープ実験仕様書 v1

**日付:** 2026-07-15  
**著者:** 木原 範昭  
**位置づけ:** 交換散乱係数R集中実験群・系統C仕様書  

---

## 1. 目的

本仕様書の目的は、20260711 のAB二体加速度様調和読出しに部分反射 `R` を導入し、調和読出しが `R` に敏感か鈍感かを調べることである。

本系統では、`R=0.70` 近傍に有効点が出ることを期待しすぎない。

重要なのは、加速度様調和読出しが `R` に鋭く反応するのか、それとも広い `R` 範囲で安定するのかを分類することである。

---

## 2. 直接の前提

直接の前提は次である。

```text
../20260711/AB二体閉鎖位相系における調和読出しとc=1面積スイープ予備実験総括 v2.md
../20260711/run_ab_two_body_fermionic_reflection_harmonic_readout_v2.py
../20260711/ab_two_body_fermionic_reflection_harmonic_readout_preliminary_result_v2/ab_two_body_fermionic_reflection_harmonic_readout_report_v2.md
```

20260713 の散乱行列写像は、部分反射 `R` の導入に用いる。

---

## 3. 基本方針

既存のAB二体加速度様調和読出しでは、完全反射極限として、

```text
R = 1
T ≈ 0
q_out_factor = -1
```

が使われていた。

本実験では、これを固定せず、`R` を掃引する。

部分反射の場合、方向読出しの圧縮表示は、

```text
q_out_factor = T - R = 1 - 2R
```

として読む。

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

これらは加速度様読出しの安定点を仮定するものではなく、系統A/Bと同じ位置で感度を確認するための固定評価点である。

必要に応じて、安定性境界が見えた箇所を細分化する。

---

## 5. 実験条件

既存のV2条件をできるだけ保持する。

```text
harmonic cases: 既存48条件
c1 cases: 既存144条件
external_c_used_any = false
f_A_or_f_B_used_any = false
```

変更するのは、フェルミオン型反跳写像に入る `R` のみである。

---

## 6. 記録量

各 `R` について、既存レポートと同等の量を記録する。

```text
R
T
Delta_F
q_out_factor
max_Q_closed_abs
harmonic_case_count
c1_case_count
fermionic_regular_cell_harmonic_consistent_all_cases
fermionic_regular_cell_harmonic_consistent_nonstrong_modes
fermionic_strong_readout_perturbs_harmonic_projection
fermionic_max_f_AB_projection_error_regular_nonstrong
fermionic_max_f_AB_projection_error_regular_strong
label_free_pass_vs_fermionic_match_all_cases
label_free_display_vs_fermionic_match_all_cases
readout_off_decay_max_abs
readout_strong_decay_min_abs
fermionic_c1_area_sweep_detected_all_cases
```

---

## 7. 判定量

本系統では、次を主判定とする。

```text
harmonic_valid_rate
c1_area_valid_rate
projection_error_nonstrong
projection_error_strong
label_free_match_rate
```

総合スコアは、調和読出しと `c=1` 面積読出しが壊れないことを優先する。

```text
score_C =
  harmonic_valid_penalty
  + c1_area_penalty
  + normalized_projection_error
  + label_free_mismatch_penalty
```

ただし、`score_C` が平坦な場合は、無理に `R_star` を定義しない。

---

## 8. R_starとR感度

本系統では、次を読む。

```text
R_star_C = argmin_R score_C
condition_at_R_star_C
distance_to_R_137
distance_to_R_128
```

ただし、主判定は `R_star_C` そのものではなく、次の分類である。

```text
sharp:
  R=0.70 近傍でのみ調和読出しが安定する。

broad:
  広い R 範囲で調和読出しが安定する。

flat:
  掃引範囲内でほぼ差がない。

broken:
  部分反射 R を入れると調和読出しが成立しない。
```

`broad` または `flat` は失敗ではない。

その場合、加速度様調和読出しは局在性交換ほど `R` に敏感ではないと読む。

---

## 9. 出力

出力ディレクトリは次とする。

```text
system_C_ab_acceleration_harmonic_R_sensitivity_result_v1
```

出力ファイルは次を基本とする。

```text
system_C_ab_acceleration_harmonic_R_sensitivity_rows_v1.csv
system_C_ab_acceleration_harmonic_R_sensitivity_summary_v1.csv
system_C_ab_acceleration_harmonic_R_sensitivity_result_v1.json
system_C_ab_acceleration_harmonic_R_sensitivity_report_v1.md
system_C_ab_acceleration_harmonic_R_sensitivity_score_v1.png
```

---

## 10. 読み方

本系統で `R=0.70` 近傍に鋭い安定点が出る場合、交換散乱係数は加速度様調和読出しにも直接関係している可能性がある。

一方、広い `R` 範囲で調和読出しが保たれる場合、`R=0.70` 近傍は加速度様読出し一般ではなく、局在性、倍音移乗、状態選択に関係する係数として読む。

この区別が、三系統比較の重要な対照になる。
