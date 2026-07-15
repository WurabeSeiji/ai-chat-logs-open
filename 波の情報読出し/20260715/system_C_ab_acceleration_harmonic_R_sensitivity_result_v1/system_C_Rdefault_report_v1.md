# 系統C AB加速度様調和読出しR感度スイープ 予備実験レポート

## 代表結果

| 量 | 値 |
|---|---:|
| `R_star_C` | `0.629` |
| `score_C_min` | `0.9949535051579352` |
| `classification_C` | `flat` |
| `q_out_factor_diagnostic_at_star` | `-0.2580000000000002` |
| `q_out_factor_applied_at_star` | `False` |
| `full_two_channel_scattering_used_rate_at_star` | `1.0` |
| `max_scattering_unitarity_error_at_star` | `6.245004513516506e-17` |
| `distance_to_R_137` | `0.06817787912824735` |
| `distance_to_R_128_nominal` | `0.05767146567112502` |
| `harmonic_valid_rate_nonstrong_at_star` | `1.0` |
| `c1_area_valid_rate_at_star` | `1.0` |
| `c1_calibrated_rate_at_star` | `0.25` |

## 上位ピーク

| rank | kind | R | score_C | q_out_factor_diagnostic |
|---:|---|---:|---:|---:|
| 1 | `local_min` | `0.629` | `0.9949535051579352` | `-0.2580000000000002` |
| 2 | `local_min` | `0.66` | `0.9949637341945585` | `-0.31999999999999995` |
| 3 | `local_min` | `0.690671465671` | `0.9949722605137112` | `-0.381342931342` |
| 4 | `local_min` | `0.613` | `0.995246547097324` | `-0.2260000000000001` |
| 5 | `local_min` | `0.645` | `0.9952609189359016` | `-0.29` |
| 6 | `local_min` | `0.675671465671` | `0.9952661372093654` | `-0.35134293134200006` |
| 7 | `local_min` | `0.706` | `0.9952767090538861` | `-0.4119999999999999` |
| 8 | `local_min` | `0.791` | `0.9953059885862237` | `-0.582` |
| 9 | `local_min` | `0.762` | `0.9953137038803801` | `-0.5240000000000002` |
| 10 | `local_min` | `0.738` | `0.9953199746448292` | `-0.47600000000000015` |
| 11 | `local_min` | `0.84` | `0.995342494782556` | `-0.6799999999999999` |
| 12 | `probe_R137` | `0.697177879128` | `0.9953541922638128` | `-0.39435575825600017` |

## 読み方

本実験は、20260711 のAB二体加速度様調和読出しに部分反射率 `R` を導入し、調和読出しと `c=1` 面積読出しがどの程度 `R` に敏感かを確認する。

`q_out_factor_diagnostic` は透過率と反射率から読まれる診断量であり、演算子として `chi_read` へ掛けていない。実装ではA/B二チャネル散乱行列を入射チャネルへ作用させる。

`score_C` は、調和読出し、c=1面積読出し、射影誤差を合わせた管理用集約値である。物理的なR地形は `projection_penalty` 図で読む。

`R_128_nominal` は高エネルギー側で `1/alpha` が128近傍へ走ることを読むための名目プローブであり、精密測定値そのものではない。

## 出力

| 種類 | ファイル |
|---|---|
| rows | `system_C_Rdefault_rows_v1.csv` |
| condition rows | `system_C_Rdefault_condition_rows_v1.csv` |
| summary | `system_C_Rdefault_summary_v1.csv` |
| best | `system_C_Rdefault_best_v1.csv` |
| peaks | `system_C_Rdefault_peaks_v1.csv` |
| result | `system_C_Rdefault_result_v1.json` |
| scores | `system_C_Rdefault_scores_v1.png` |
| projection landscape overview | `system_C_Rdefault_stability_depth_distribution_overview_v1.png` |
| projection landscape deep | `system_C_Rdefault_stability_depth_distribution_deep_v1.png` |
| peak zoom | `system_C_Rdefault_peak_zoom_v1.png` |
