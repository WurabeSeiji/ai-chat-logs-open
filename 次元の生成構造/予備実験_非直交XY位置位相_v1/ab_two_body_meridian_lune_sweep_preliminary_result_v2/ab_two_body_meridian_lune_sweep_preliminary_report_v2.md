# AB二体・経度±2.5度大円リューネ予備実験 v2

## 1. R²保存の確認

各経度大円について、全ステップで

```text
x² + y² + z² = R²
```

を直接検査した。

- R2_preserved_all_cases: `True`
- max_R2_abs_error: `6.6613381477509392e-16`
- max_R2_rel_error: `6.0744119919478487e-16`
- cycle_returns_all_cases: `True`
- max_cycle_return_error: `7.7477341530709071e-15`

## 2. 二大円と交点

- longitude_minus_deg: `-2.5`
- longitude_plus_deg: `2.5`
- longitude_separation_deg: `5.0`
- nodes_intersect_all_cases: `True`
- max_node_intersection_error: `4.9238120677277192e-16`

二大円は位相0度と180度で交差する。

## 3. 空間的な球面スイープ

- lune_area_detected_all_cases: `True`
- max_lune_area_analytic_rel_error: `2.5491076916181669e-16`
- max_lune_area_cycle_return_error: `1.5449880957918438e-17`
- time_axis_area_used_any: `False`

面積はXT・YT投影ではなく、二つの実三成分大円が挟む球面積として計算した。

## 4. 距離指数

- `max_f_native_mean`: slope=1, alpha=-1, log_rmse=1.426e-14
- `rms_f_native_mean`: slope=1, alpha=-1, log_rmse=6.714e-16
- `max_A_meridian_lune`: slope=2, alpha=-2, log_rmse=1.760e-15
- `max_K_area_times_f`: slope=3, alpha=-3, log_rmse=7.018e-15

- native_inverse_square_detected: `False`
- inverse_term_used_any: `False`
- readout_loss_used_any: `False`

主要なネイティブ加速度フィット：

```json
{
  "candidate": "max_f_native_mean",
  "fit_valid": true,
  "case_count": 7,
  "loglog_slope": 1.0000000000000104,
  "inverse_power_alpha": -1.0000000000000104,
  "loglog_intercept": -5.45329923639949,
  "log_rmse": 1.4262307636984606e-14,
  "min_value": 7.47376780486225e-05,
  "max_value": 0.0044842606829176195,
  "scattering_protocol": "fermionic_reflection_pi",
  "inverse_square_match": false
}
```
