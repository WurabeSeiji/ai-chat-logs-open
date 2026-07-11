# ABCベースライン定常波C 準安定G1窓スイープ予備実験 v1

## 目的

`G1`: `R*aχ ≈ -∂χE_read` 候補が、定常的な恒等式ではなく準安定遷移窓として現れるかを調べる。

前段の G1 予備実験では、C単独では整合し、R2/R3反力対照では棄却され、合成状態から R2/R3 を差し引いた residual は early window で整合したが full window では破れた。

本予備実験では、C記憶減衰、反力記憶減衰、C源強度、C戻り強度、反力候補強度を掃引し、G1 residual の準安定窓長を分類する。

## 統合判定

- case_count: `17`
- all_cases_classified: `True`
- single_gauge_only_used: `False`
- C_only_full_window_valid_all_cases: `True`
- residual_cases_with_finite_break: `14`
- residual_cases_without_break: `3`
- base_residual_window_length: `27`
- base_residual_first_validity_break_step: `28`
- min_residual_window_length: `7`
- median_residual_window_length: `27.0`
- max_residual_window_length: `46`
- min_C_only_window_length: `46`
- max_C_only_window_length: `46`
- max_C_only_error: `9.736820538544116e-23`
- max_residual_inside_window_error: `6.920064035294847e-11`
- max_reaction_dominance_ratio_abs: `8243.181592047982`
- memory_decay_C_window_correlation: `0.6914965350282571`
- memory_decay_R_window_correlation: `-0.3251567909658586`
- reaction_scale_window_correlation: `-0.5165645763062486`
- metastable_window_sweep_preliminary_valid: `True`

## ケース別サマリー

| case | C decay | R decay | C source | C return | reaction scale | residual window | break | reaction/C | residual after R2/R3 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| base | 0.920 | 0.900 | 1.00e-06 | 1.00e-03 | 1.000 | 27 | 28 | 1.814796e+03 | -1.1666047494696130e-07 |
| C_decay_0_78 | 0.780 | 0.900 | 1.00e-06 | 1.00e-03 | 1.000 | 11 | 12 | 4.121261e+03 | -4.1567180919521718e-08 |
| C_decay_0_86 | 0.860 | 0.900 | 1.00e-06 | 1.00e-03 | 1.000 | 16 | 17 | 2.785073e+03 | -6.9915269462761387e-08 |
| C_decay_0_97 | 0.970 | 0.900 | 1.00e-06 | 1.00e-03 | 1.000 | 46 |  | 1.078700e+03 | -2.0821756419042714e-07 |
| R_decay_0_78 | 0.920 | 0.780 | 1.00e-06 | 1.00e-03 | 1.000 | 43 | 44 | 9.388735e+02 | -1.2837818064159734e-07 |
| R_decay_0_96 | 0.920 | 0.960 | 1.00e-06 | 1.00e-03 | 1.000 | 20 | 21 | 3.180885e+03 | -9.4791389182447006e-08 |
| C_source_half | 0.920 | 0.900 | 5.00e-07 | 1.00e-03 | 1.000 | 21 | 22 | 3.629591e+03 | -4.9568505300356946e-08 |
| C_source_double | 0.920 | 0.900 | 2.00e-06 | 1.00e-03 | 1.000 | 34 | 35 | 9.073979e+02 | -2.5084439847500306e-07 |
| C_return_half | 0.920 | 0.900 | 1.00e-06 | 5.00e-04 | 1.000 | 21 | 22 | 3.629591e+03 | -4.9568505300356946e-08 |
| C_return_double | 0.920 | 0.900 | 1.00e-06 | 2.00e-03 | 1.000 | 34 | 35 | 9.073979e+02 | -2.5084439847500306e-07 |
| reaction_quarter | 0.920 | 0.900 | 1.00e-06 | 1.00e-03 | 0.250 | 46 |  | 4.536717e+02 | -1.3305930457452320e-07 |
| reaction_half | 0.920 | 0.900 | 1.00e-06 | 1.00e-03 | 0.500 | 42 | 43 | 9.073615e+02 | -1.2978355456993995e-07 |
| reaction_double | 0.920 | 0.900 | 1.00e-06 | 1.00e-03 | 2.000 | 15 | 16 | 3.629882e+03 | -6.4125397103431681e-08 |
| slow_C_weak_reaction | 0.970 | 0.900 | 1.00e-06 | 1.00e-03 | 0.250 | 46 |  | 2.696587e+02 | -2.2460302784299557e-07 |
| slow_C_strong_reaction | 0.970 | 0.900 | 1.00e-06 | 1.00e-03 | 2.000 | 30 | 31 | 2.157572e+03 | -1.5570030909017873e-07 |
| fast_C_weak_reaction | 0.780 | 0.900 | 1.00e-06 | 1.00e-03 | 0.250 | 20 | 21 | 1.030253e+03 | -5.7980404921664785e-08 |
| fast_C_strong_reaction | 0.780 | 0.900 | 1.00e-06 | 1.00e-03 | 2.000 | 7 | 8 | 8.243182e+03 | 1.0987090570679925e-08 |

## 解釈

- `C_only_full_window_valid` が保たれる一方で、`combined_minus_R2_R3_residual` が有限窓で破れる場合、G1 は単純な常時成立式ではなく、相互作用後の準安定遷移で読まれる候補である。
- 反力候補が強い条件で窓が短くなる場合、観測窓の切り方により G1 候補が隠れる。
- C記憶が長い条件で窓が伸びる場合、G1 は C残渣遅延と結びついた候補である。
- 本実験の `valid` は、全ケースで G1 が成立したという意味ではない。準安定窓を多条件で分類できたという意味である。
- これは標準重力の導出ではなく、加速度読出し前の窓長・遅延・反力分離の予備計量である。

## 出力

| 種類 | ファイル |
|---|---|
| JSON | `abc_baseline_stationary_wave_metastable_window_sweep_preliminary_result_v1.json` |
| case CSV | `abc_baseline_stationary_wave_metastable_window_sweep_cases_v1.csv` |
| rows CSV | `abc_baseline_stationary_wave_metastable_window_sweep_rows_v1.csv` |
