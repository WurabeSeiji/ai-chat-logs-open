# AB二体閉鎖位相系におけるc=1内部較正と空間位相・時間位相面積スイープ予備実験検証メモ v1

**日付:** 2026-07-12  
**著者:** 木原 範昭  
**位置づけ:** 波の情報読出しシリーズ・予備実験検証メモ  

---

## 1. 対象

本メモは、次の実験仕様書に対応する予備実験の実行結果を記録する。

[AB二体閉鎖位相系におけるc=1内部較正と空間位相・時間位相面積スイープ読出し実験仕様書 v1.md](AB二体閉鎖位相系におけるc=1内部較正と空間位相・時間位相面積スイープ読出し実験仕様書%20v1.md)

本予備実験の目的は、一角度 AB 調和読出しで距離減衰が出なかった結果を受けて、実験ステップ `s` と時間位相読出し `tau_read` を分離し、内部 `c=1` 較正のもとで `chi-tau` 面積スイープが成立するかを検査することである。

---

## 2. 実行コマンド

```text
python3 run_ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_v1.py
```

出力ディレクトリ:

```text
ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1/
```

---

## 3. 出力ファイル

| 種類 | ファイル |
|---|---|
| レポート | [ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_report_v1.md](ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1/ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_report_v1.md) |
| JSON | [ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1.json](ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1/ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1.json) |
| 時系列 CSV | [ab_two_body_c1_internal_calibration_chi_tau_area_sweep_series_v1.csv](ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1/ab_two_body_c1_internal_calibration_chi_tau_area_sweep_series_v1.csv) |
| ケース要約 CSV | [ab_two_body_c1_internal_calibration_chi_tau_area_sweep_case_summary_v1.csv](ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1/ab_two_body_c1_internal_calibration_chi_tau_area_sweep_case_summary_v1.csv) |
| 逆冪候補 CSV | [ab_two_body_c1_internal_calibration_chi_tau_area_sweep_power_candidates_v1.csv](ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1/ab_two_body_c1_internal_calibration_chi_tau_area_sweep_power_candidates_v1.csv) |
| `chi-tau` 軌跡図 | [ab_two_body_c1_internal_calibration_chi_tau_surface_v1.png](ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1/ab_two_body_c1_internal_calibration_chi_tau_surface_v1.png) |
| 面積スイープ図 | [ab_two_body_c1_internal_calibration_area_sweep_v1.png](ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1/ab_two_body_c1_internal_calibration_area_sweep_v1.png) |
| `c=1` 較正誤差図 | [ab_two_body_c1_internal_calibration_error_v1.png](ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1/ab_two_body_c1_internal_calibration_error_v1.png) |
| 冪候補図 | [ab_two_body_c1_internal_calibration_power_candidate_v1.png](ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1/ab_two_body_c1_internal_calibration_power_candidate_v1.png) |
| 読出しリーク図 | [ab_two_body_c1_internal_calibration_readout_leak_v1.png](ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1/ab_two_body_c1_internal_calibration_readout_leak_v1.png) |

---

## 4. 結果サマリー

主要結果は次である。

```text
c1_internal_calibration_chi_tau_area_sweep_preliminary_valid: true
case_summary_count: 288
power_candidate_count: 48
max_Q_closed_abs: 0.0
disabled_max_area: 0.0
locked_max_area: 7.105427357601002e-15
independent_min_area: 0.0024367633602631385
c1_max_epsilon_c_abs: 0.003062925647902448
c1_readout_off_max_epsilon_c_abs: 2.3314683517128287e-15
c1_area_sweep_detected_all_cases: true
tau_disabled_rank1_all_cases: true
tau_locked_rank1_all_cases: true
tau_independent_rank2_all_cases: true
readout_off_decay_max_abs: 4.4544900995234305e-18
readout_strong_decay_min_abs: 0.00040004000533408923
c1_readout_off_power_candidate_alpha_values: [-2.0000000000000013, -2.0000000000000013]
tau_is_step_used_any: false
external_c_used_any: false
f_A_or_f_B_used_any: false
```

---

## 5. 読み

### 5.1 `s` と `tau_read` の分離

本実験では、実験ステップ `s` を `tau_read` として使っていない。

`tau_is_step_used_any` は `false` である。

したがって、本予備実験は、外部時間をそのまま持ち込む実験ではなく、AB 閉鎖位相系内で時間位相候補を別に読む実験である。

### 5.2 tau disabled / tau locked は面積を作らない

`tau_disabled_control` では、

```text
disabled_max_area = 0.0
```

となった。

また、`tau_locked_0deg`, `tau_locked_90deg` では、

```text
locked_max_area = 7.105427357601002e-15
```

であり、数値丸め範囲で面積ゼロと読める。

さらに、

```text
tau_disabled_rank1_all_cases = true
tau_locked_rank1_all_cases = true
```

であった。

これは、`tau` という名前を付けただけでは `chi-tau` 面が立たないことを示す。

### 5.3 tau independent は面積を作る

独立時間位相候補を置いた場合、

```text
tau_independent_rank2_all_cases = true
independent_min_area = 0.0024367633602631385
```

となった。

したがって、`chi_read` と `tau_read` が独立な読出しとして構成された場合、`chi-tau` 面積スイープは数値的に検出された。

### 5.4 c=1 内部較正

`tau_independent_c1` の `readout_off` 条件では、

```text
c1_readout_off_max_epsilon_c_abs = 2.3314683517128287e-15
```

であった。

これは、一周期 RMS 交換比としての内部 `c=1` 較正が、読出し波によるリークを止めた条件で成立することを示す。

一方、全読出し条件を含めると、

```text
c1_max_epsilon_c_abs = 0.003062925647902448
```

まで歪む。

これは、読出し波が `c=1` 内部較正をわずかに乱す候補として記録する。

### 5.5 読出し波による減衰

`readout_off` では、

```text
readout_off_decay_max_abs = 4.4544900995234305e-18
```

であり、実質的に減衰なしである。

一方、`readout_strong` では、

```text
readout_strong_decay_min_abs = 0.00040004000533408923
```

となった。

したがって、今回も読出し波による包絡減衰候補が確認された。

### 5.6 逆冪候補

`chi-tau` 面積スイープが成立した後にのみ、冪候補をフィットした。

`tau_independent_c1`, `readout_off` では、

```text
alpha = -2
```

であった。

これは、読出し強度候補が距離で減衰するのではなく、初期位相偏差に対して二乗的に増加することを意味する。

したがって、本予備実験は `chi-tau` 面の成立を確認したが、逆比例型または逆二乗型の距離減衰を確認したものではない。

---

## 6. 判定

本予備実験は `valid` である。

成立したことは次である。

1. `s` と `tau_read` は分離された。
2. `tau_disabled` と `tau_locked` は面積を作らない。
3. `tau_independent` は `chi-tau` 面積を作る。
4. `tau_independent_c1` は `readout_off` 条件で内部 `c=1` 較正を満たす。
5. 読出し波が強いほど、包絡減衰と較正歪みが現れる。

成立していないことは次である。

1. 逆比例型の距離減衰。
2. 逆二乗型の距離減衰。
3. 標準重力対応。

---

## 7. 次の課題

次に検査すべき候補は、次である。

1. `tau_independent_c1` の位相差、振幅比、周期比をさらに掃引する。
2. `chi-tau` 面積が読出し波リークでどのように縮むかを詳細化する。
3. 面積スイープが成立しても逆冪減衰が出ない理由を整理する。
4. 位置位相自由度を二つに増やした `S^2` 型読出しに進むか、または `chi-tau` 面に追加の伝搬写像を入れる必要があるかを判定する。

現時点の最重要結論は、

```text
時間位相を独立に読めば chi-tau 面は立つ。
しかし、それだけでは距離減衰は出ない。
```

である。
