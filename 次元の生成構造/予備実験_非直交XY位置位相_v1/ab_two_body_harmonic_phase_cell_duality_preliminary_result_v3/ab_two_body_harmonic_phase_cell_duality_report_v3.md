# AB二体・調和位相セル双対予備実験 v3

## 実装条件

- fixed_radius_R: `1.0`
- base_period_steps: `12288`
- common_step_count: `720`
- harmonic_orders: `[1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96]`
- phase_cell_width: `L_n = 2 pi / n`
- harmonic_frequency: `omega_n = n omega_1`
- `1/L^2` による除算: `False`
- 正規化: `False`
- 追加空間軸・時間軸・面積項: `False`

## 統合結果

- R2_preserved_all_cases: `True`
- max_R2_abs_error: `2.2204460492503131e-16`
- harmonic_duality_preserved_all_cases: `True`
- max_omega_times_L_abs_error: `0`
- cycles_close_all_cases: `True`
- max_cycle_return_error: `2.4492935982947064e-16`
- max_scattering_unitarity_error: `0`
- max_tangent_acceleration: `7.229644567250593e-15`
- max_pass_vs_fermionic_acceleration_difference: `0`

## 冪回帰

- `pass_through`: acceleration vs L slope = `-1.9999708997597263`
- `fermionic_reflection_pi`: acceleration vs L slope = `-1.9999708997597263`
- `pass_through`: omega vs L slope = `-1.0000000000000011`
- `fermionic_reflection_pi`: omega vs L slope = `-1.0000000000000011`

## 判定範囲

本実験は、閉じた位相円で `L_n = 2 pi / n` と `omega_n = n omega_1` を
同時に用いた場合、固定 `R` の円運動から直接計算した離散二階差分が
`L_n^-2` に近い冪を示すかを検査する。逆二乗項は演算へ入力していない。
高倍音では離散二階差分 `4 R sin^2(omega_n/2)` と連続近似
`R omega_n^2` の差が生じるため、その相対誤差を各ケースへ記録する。
