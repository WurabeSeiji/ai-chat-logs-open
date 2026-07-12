# AB二体閉鎖位相系における一角度円周位相調和読出しパラメータスイープ予備実験検証メモ v1

## 目的

AB一角度円周位相調和読出し予備実験について、初期偏差、回転周期、読出し漏れ量を変え、結果の安定範囲を調べた。

本スイープでも、観測機 C、標準重力式、標準ばね式、`f_A`, `f_B` は使わない。

## 統合判定

| 量 | 値 |
|---|---:|
| `sweep_configuration_count` | `210` |
| `case_summary_count` | `420` |
| `period_count` | `5` |
| `deviation_count` | `7` |
| `leak_count` | `6` |
| `protocol_count` | `2` |
| `observer_C_used` | `False` |
| `absolute_background_axis_used` | `False` |
| `f_A_or_f_B_used` | `False` |
| `standard_force_law_used` | `False` |
| `max_Q_closed_abs` | `0.0` |
| `max_D_AB_near_protocol_diff` | `0.0` |
| `max_V_AB_protocol_diff` | `0.0` |
| `label_free_protocol_degenerate_all_cases` | `True` |
| `oscillation_detected_all_cases` | `True` |
| `readout_off_decay_max_abs` | `8.255309734203768e-17` |
| `decay_abs_monotonic_by_leak_all_grids` | `True` |
| `normalized_f_error_monotonic_by_leak_all_grids` | `True` |
| `max_normalized_f_AB_projection_error` | `0.0002580660255432312` |
| `leak_1e3_min_normalized_f_AB_projection_error` | `6.245465172379202e-05` |
| `strong_leak_detection_floor` | `5e-05` |
| `strong_leak_perturbs_projection_all_cases` | `True` |
| `parameter_sweep_preliminary_valid` | `True` |

## 漏れ量別サマリー

| readout | leak | max normalized f error | mean normalized f error | max |decay| | min envelope ratio |
|---|---:|---:|---:|---:|---:|
| readout_off | 0.0e+00 | 2.3191544552338146e-16 | 1.7828056158834067e-16 | 8.2553097342037681e-17 | 9.9999999999995037e-01 |
| leak_1e-6 | 1.0e-06 | 2.6104938238942896e-07 | 1.4376634743745425e-07 | 2.0000010002041215e-06 | 9.9693271223251489e-01 |
| leak_1e-5 | 1.0e-05 | 2.6102236499488829e-06 | 1.4373875419886592e-06 | 2.0000100000661670e-05 | 9.6974691528965185e-01 |
| leak_5e-5 | 5.0e-05 | 1.3045115939866196e-05 | 7.1808098223663267e-06 | 1.0000250008349528e-04 | 8.5761170506577911e-01 |
| leak_2e-4 | 2.0e-04 | 5.2090519060146818e-05 | 2.8631549818212146e-05 | 4.0004000533416209e-04 | 5.4093213685773500e-01 |
| leak_1e-3 | 1.0e-03 | 2.5806602554323120e-04 | 1.4075361448771432e-04 | 2.0010006671671585e-03 | 4.6257252117821133e-02 |

## 周期別サマリー

| period_steps | omega_step | max normalized f error | mean normalized f error | max |decay| |
|---:|---:|---:|---:|---:|
| 48 | 1.3089969389957470e-01 | 2.5806602554323120e-04 | 5.4345488929281856e-05 | 2.0010006671671585e-03 |
| 72 | 8.7266462599716474e-02 | 1.7128609259666383e-04 | 3.6108607086114093e-05 | 2.0010006671670670e-03 |
| 96 | 6.5449846949787352e-02 | 1.2776666270611215e-04 | 2.6962737282054558e-05 | 2.0010006671670397e-03 |
| 144 | 4.3633231299858237e-02 | 8.4194639868799927e-05 | 1.7805409837181665e-05 | 2.0010006671669117e-03 |
| 192 | 3.2724923474893676e-02 | 6.2454651723800627e-05 | 1.3233696880282171e-05 | 2.0010006671669629e-03 |

## 解釈

- `Protocol F/B` は全スイープで `D_AB` と `V_AB` では縮退した。
- `readout_off` の減衰は数値丸め範囲に留まった。
- 読出し漏れ量を増やすと、包絡減衰と `f_AB` 射影不整合が単調に増えた。
- したがって、AB一角度版は弱読出しでは安定だが、強い読出し波は補償表示そのものを歪める。
- これは、後続の二角度・三角度実験で観測波を弱く保つ必要があることを示す制御結果である。

## 出力

| 種類 | ファイル |
|---|---|
| JSON | `ab_two_body_one_angle_harmonic_readout_parameter_sweep_preliminary_result_v1.json` |
| case summary CSV | `ab_two_body_one_angle_parameter_sweep_case_summary_v1.csv` |
| protocol comparison CSV | `ab_two_body_one_angle_parameter_sweep_protocol_comparison_v1.csv` |
| leak summary CSV | `ab_two_body_one_angle_parameter_sweep_leak_summary_v1.csv` |
| period summary CSV | `ab_two_body_one_angle_parameter_sweep_period_summary_v1.csv` |
| monotonic checks CSV | `ab_two_body_one_angle_parameter_sweep_monotonic_checks_v1.csv` |
| selected series CSV | `ab_two_body_one_angle_parameter_sweep_selected_series_v1.csv` |
| f error plot | `ab_two_body_one_angle_parameter_sweep_leak_f_error_v1.png` |
| decay plot | `ab_two_body_one_angle_parameter_sweep_leak_decay_v1.png` |
| period plot | `ab_two_body_one_angle_parameter_sweep_period_v1.png` |
| selected series plot | `ab_two_body_one_angle_parameter_sweep_selected_series_v1.png` |
