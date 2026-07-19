# AB二体・非直交XY位置位相予備実験 v1

## 実装

既存の `run_ab_two_body_fermionic_reflection_harmonic_readout_v4.py` から、
AB基礎回転、二チャネル散乱行列、面積・実ランク計算をそのまま使用した。

対称関係を `delta=0 deg` とし、第二位置位相を
`delta=-5 deg, +5 deg` だけ外した。二位置位相には総振幅を等分し、
読出し損失、外部軸、逆数、逆面積、逆二乗項は使用していない。

## 判別結果

- case_count: `42`
- zero_offset_rank1_all: `True`
- zero_offset_area_zero_all: `True`
- nonzero_offset_rank2_all: `True`
- nonzero_offset_area_detected_all: `True`
- signed_area_plus_minus_antisymmetry_max_rel_error: `1.7398171804292282e-14`
- max_relation_norm_rel_error: `8.0992159892637980e-16`
- max_scattering_unitarity_error: `1.6653345369377348e-16`
- native_inverse_square_detected: `False`

## fermionic reflection, delta=+5 deg の冪指数

- `max_f_xy_second_difference`: slope=1, alpha=-1, log_rmse=3.188e-14
- `max_abs_A_XY`: slope=2, alpha=-2, log_rmse=5.638e-15
- `max_K_area_times_f`: slope=3, alpha=-3, log_rmse=2.982e-14

## 読み

`delta=0 deg` は一次元へ退化し、XY有向面積を生成しない。
`delta=+/-5 deg` は実ランク2と非零の有向面積を生成し、その符号は反転する。

この最小追加だけで既存の二階差分が逆二乗へ変化したかどうかは、
`native_inverse_square_detected` で判定する。最も指数2に近いネイティブ結果は次である。

```json
{
  "candidate": "max_f_xy_second_difference",
  "fit_valid": true,
  "case_count": 7,
  "loglog_slope": 0.9999999999999917,
  "inverse_power_alpha": -0.9999999999999917,
  "loglog_intercept": -5.454489015482266,
  "log_rmse": 1.339438347704979e-14,
  "min_value": 7.464880959998421e-05,
  "max_value": 0.004478928575998784,
  "scattering_protocol": "pass_through",
  "relation_offset_deg": -5.0,
  "inverse_square_match": false
}
```
