# 系統C AB加速度様調和読出しR感度スイープ 予備実験レポート

## 代表結果

| 量 | 値 |
|---|---:|
| `R_star_C` | `0.7495` |
| `score_C_min` | `0.9949874006934684` |
| `classification_C` | `flat` |
| `q_out_factor_diagnostic_at_star` | `-0.4990000000000003` |
| `q_out_factor_applied_at_star` | `False` |
| `full_two_channel_scattering_used_rate_at_star` | `1.0` |
| `max_scattering_unitarity_error_at_star` | `7.632783294297951e-17` |
| `distance_to_R_137` | `0.0523221208717527` |
| `distance_to_R_128_nominal` | `0.06282853432887503` |
| `harmonic_valid_rate_nonstrong_at_star` | `1.0` |
| `c1_area_valid_rate_at_star` | `1.0` |
| `c1_calibrated_rate_at_star` | `0.25` |

## 上位ピーク

| rank | kind | R | score_C | q_out_factor_diagnostic |
|---:|---|---:|---:|---:|
| 1 | `local_min` | `0.7495` | `0.9949874006934684` | `-0.4990000000000003` |
| 2 | `local_min` | `0.752` | `0.9953865492025344` | `-0.5040000000000001` |
| 3 | `local_min` | `0.748` | `0.9953872204173` | `-0.49600000000000005` |

## 読み方

本実験は、20260711 のAB二体加速度様調和読出しに部分反射率 `R` を導入し、調和読出しと `c=1` 面積読出しがどの程度 `R` に敏感かを確認する。

`q_out_factor_diagnostic` は透過率と反射率から読まれる診断量であり、演算子として `chi_read` へ掛けていない。実装ではA/B二チャネル散乱行列を入射チャネルへ作用させる。

`score_C` は、調和読出し、c=1面積読出し、射影誤差を合わせた管理用集約値である。物理的なR地形は `projection_penalty` 図で読む。

`R_128_nominal` は高エネルギー側で `1/alpha` が128近傍へ走ることを読むための名目プローブであり、精密測定値そのものではない。

## 出力

| 種類 | ファイル |
|---|---|
| rows | `system_C_R075_dense_rows_v1.csv` |
| condition rows | `system_C_R075_dense_condition_rows_v1.csv` |
| summary | `system_C_R075_dense_summary_v1.csv` |
| best | `system_C_R075_dense_best_v1.csv` |
| peaks | `system_C_R075_dense_peaks_v1.csv` |
| result | `system_C_R075_dense_result_v1.json` |
| scores | `system_C_R075_dense_scores_v1.png` |
| projection landscape overview | `system_C_R075_dense_stability_depth_distribution_overview_v1.png` |
| projection landscape deep | `system_C_R075_dense_stability_depth_distribution_deep_v1.png` |
| peak zoom | `system_C_R075_dense_peak_zoom_v1.png` |
