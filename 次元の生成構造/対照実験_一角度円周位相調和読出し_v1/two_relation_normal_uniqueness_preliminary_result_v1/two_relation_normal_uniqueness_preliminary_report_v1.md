# 二体位相関係からの第三方向一意性予備実験検証メモ v1

## 目的

無改変対照テストを通過した一角度円周位相調和読出しをスイープ核として使い、二つの無名な位相関係から第三方向を追加の特権軸なしに一意に再構成できる条件を検査した。

## 検査方法

二つの関係ベクトルを `u`, `v` とし、両者の位相差を `5, 15, 30, 60 deg` とした。各ケースをランダム直交変換で基底交換し、二体関係が張る行列

```math
A = [u^T; v^T]
```

の零空間を法線候補集合として計算した。候補集合への射影は、

```math
P_perp = I - A^T (A A^T)^{-1} A
```

である。選択に外部軸名は使用していない。

## 一角度核の保存検算

- `step_count`: `720`
- `period_steps`: `96`
- `max_R2_drift`: `2.9056618222611519e-16`
- `max_one_period_repeat_error`: `3.3659000130511772e-16`

## 次元別結果

| test d | normal nullity | unique third line | unique readable directions | XYZ rank | candidate angle min | linear range min | quadratic range max | pass |
|---:|---:|---|---:|---:|---:|---:|---:|---|
| 2 | 0 | False | 2 | 2 | 0.00000000 | 0.00000000 | 0.00000000e+00 | True |
| 3 | 1 | True | 3 | 3 | 0.00000000 | 0.00000000 | 0.00000000e+00 | True |
| 4 | 2 | False | 2 | 3 | 90.00000000 | 2.00000000 | 1.33226763e-15 | True |
| 5 | 3 | False | 2 | 3 | 90.00000000 | 2.00000000 | 1.11022302e-15 | True |
| 6 | 4 | False | 2 | 3 | 90.00000000 | 2.00000000 | 1.11022302e-15 | True |

## 統合判定

| 量 | 値 |
|---|---:|
| `experiment` | `two_relation_normal_uniqueness_preliminary_v1` |
| `trial_case_count` | `1280` |
| `one_angle_control_kernel_used` | `True` |
| `one_angle_control_max_R2_drift` | `2.905661822261152e-16` |
| `one_angle_control_max_one_period_repeat_error` | `3.365900013051177e-16` |
| `two_relation_phase_offset_includes_5deg` | `True` |
| `dimension_3_unique_third_direction_all` | `True` |
| `dimension_4plus_nonunique_normal_family_all` | `True` |
| `dimension_4plus_pair_signatures_degenerate_all` | `True` |
| `dimension_4plus_linear_components_vary_all` | `True` |
| `dimension_4plus_quadratic_norm_invariant_all` | `True` |
| `max_projector_covariance_error` | `1.0224530215383721e-14` |
| `max_label_swap_projector_error` | `1.2251281692220832e-14` |
| `max_mapped_one_angle_R2_drift` | `3.139849491518021e-16` |
| `absolute_background_axis_used_for_selection` | `False` |
| `external_direction_name_used_for_selection` | `False` |
| `ambient_dimension_used_as_test_parameter` | `True` |
| `physical_dimension_selection_derived` | `False` |
| `imaginary_axis_identification_tested` | `False` |
| `preliminary_experiment_valid` | `True` |

## 観測事実

- `d=3` では二関係平面の法線候補空間が1次元となり、符号を同一視した第三方向は全ケースで一意だった。二関係と法線候補を並べた行列のランクは3だった。
- `d>=4` では法線候補空間が2次元以上となり、同じ二体関係読出しを持つ互いに異なる候補が残った。
- `d>=4` の候補族では、任意の診断基底に沿う一次成分は変化したが、候補成分の二乗和は数値精度内で一定だった。
- 法線候補射影はランダム基底交換に対して共変であり、二関係のラベル交換でも変化しなかった。

## 分類

- `d=3` における第三方向の一意性: 本予備実験の数値結果。
- `d>=4` における二体読出しだけからの法線選択不能: 本予備実験の数値結果。
- 一次候補を選べない場合にも二乗和が読めること: 本予備実験の数値結果。
- なぜ完全系が `d=3` の表示を選ぶか: 本実験では未導出。`d` は比較用パラメータとして与えた。
- 選択されない候補を虚数軸と同定すること: 本実験では未検査。

## 出力

- `two_relation_normal_uniqueness_preliminary_result_v1.json`
- `two_relation_normal_uniqueness_trials_v1.csv`
- `two_relation_normal_uniqueness_dimension_summary_v1.csv`
- `two_relation_dimension_uniqueness_v1.png`
- `dimension4_hidden_linear_and_quadratic_readout_v1.png`
- `copied_one_angle_kernel_R2_control_v1.png`
