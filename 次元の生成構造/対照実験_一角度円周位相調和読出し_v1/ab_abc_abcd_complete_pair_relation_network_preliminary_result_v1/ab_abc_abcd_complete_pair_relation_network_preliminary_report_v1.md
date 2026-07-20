# AB・ABC・ABCD完全二体関係波ネットワーク予備実験検証メモ v1

## 目的

個体A・B・C・Dではなく、個体間の全二体関係を物理的な関係波として状態変数に置き、ABCの三関係 `AB`, `BC`, `CA` を三つの軸成分として構成できるかを検査した。

CおよびDを観測器として使用せず、観測減衰と正規化を行わない。

## 状態と閉鎖

関係波の集合を `X_e` とし、各ステップで次を直接計算した。

```math
E = sum_e ((Re X_e)^2 - (Im X_e)^2)
```

```math
F = 2 sum_e (Re X_e)(Im X_e)
```

```math
sum_e X_e^2 = E + i F = R^2
```

比較用に、実数二乗和 `sum_e |X_e|^2` も独立に記録した。

## 作業更新則

関係波同士が端点を共有する場合だけ結合し、初期位相差から実反対称生成子を作った。

```math
K_ef = adjacency(e,f) sin(theta_f - theta_e)
```

実反対称生成子のCayley変換による実直交更新を反復した。この更新は二乗形式と実数二乗和を保存する。

この更新則は保存的な作業仮説であり、第0・第1公理からの導出結果ではない。

## 構成別結果

| system | relation waves | generator rank | nullity | rotation planes | frequencies | kernel drift max | pass |
|---|---:|---:|---:|---:|---:|---:|---|
| AB | 1 | 0--0 | 1--1 | 0--0 | 0--0 | 0.0000000000000000e+00 | True |
| ABC | 3 | 2--2 | 1--1 | 1--1 | 1--1 | 4.0730280738750208e-14 | True |
| ABCD | 6 | 6--6 | 0--0 | 3--3 | 3--3 | 0.0000000000000000e+00 | True |

## 統合判定

| 量 | 値 |
|---|---:|
| `experiment` | `ab_abc_abcd_complete_pair_relation_network_preliminary_v1` |
| `AB_relation_wave_count` | `1` |
| `ABC_relation_wave_count` | `3` |
| `ABCD_relation_wave_count` | `6` |
| `ABC_three_relation_axes_constructed` | `True` |
| `ABC_all_three_relation_waves_active_all_trials` | `True` |
| `AB_single_relation_stationary_all_trials` | `True` |
| `ABCD_six_relation_directions_present` | `True` |
| `AB_generator_rank_range` | `[0, 0]` |
| `ABC_generator_rank_range` | `[2, 2]` |
| `ABC_generator_nullity_range` | `[1, 1]` |
| `ABC_rotation_plane_count_range` | `[1, 1]` |
| `ABC_independent_rotation_frequency_count_range` | `[1, 1]` |
| `ABC_one_plane_plus_one_normal_all_trials` | `True` |
| `ABCD_generator_rank_range` | `[6, 6]` |
| `ABCD_generator_nullity_range` | `[0, 0]` |
| `ABCD_rotation_plane_count_range` | `[3, 3]` |
| `ABCD_independent_rotation_frequency_count_range` | `[3, 3]` |
| `ABCD_unique_invariant_normal_all_trials` | `False` |
| `max_closure_target_error_abs` | `1.9218001973240896e-13` |
| `max_hermitian_amplitude_drift` | `2.41140440948584e-13` |
| `max_label_covariance_error` | `1.4608516699761366e-13` |
| `max_kernel_projection_drift` | `4.073028073875021e-14` |
| `max_kernel_projector_label_covariance_error` | `7.591035761743328e-16` |
| `all_relations_are_physical_waves` | `True` |
| `observer_C_or_D_used` | `False` |
| `normalization_applied` | `False` |
| `absolute_background_axis_used` | `False` |
| `relation_to_spatial_axis_is_model_definition` | `True` |
| `three_spatial_dimensions_derived` | `False` |
| `ABCD_six_to_three_projection_resolved` | `False` |
| `preliminary_experiment_valid` | `True` |

## 観測事実

- ABは一つの関係波だけを持ち、本更新則では連続混合相手がないため定常だった。
- ABCは `AB`, `BC`, `CA` の三つの物理的関係波を持ち、全試行で三波すべてが変動した。
- ABC生成子のランクは `2--2`、零空間次元は `1--1` であり、全試行で一回転平面と一不変法線へ分解した。
- ABCDは六つの物理的関係波を持ち、三軸を超える関係方向が代数的には存在した。
- ABCD生成子のランクは `6--6`、零空間次元は `0--0`、回転平面数は `3--3` だった。
- 全構成で `sum_e X_e^2 = R^2` と実数二乗和は数値精度内で保存された。
- 個体名の置換に対し、生成子、更新行列、軌道は数値精度内で共変だった。
- 零空間射影は軌道上で数値精度内に保存され、名称置換に対して共変だった。

## 分類

- ABC三関係波を三軸成分として置くこと: 本実験のモデル定義。
- ABC三関係波が閉鎖を保存しながら同時に振動できること: 本予備実験の数値結果。
- ABCの反対称生成子が一回転平面と一意な不変法線を持つこと: 本予備実験の数値結果。
- 不変法線を物理的な第三空間方向と同一視すること: 本実験では物理的解釈。
- 三関係波が物理的なXYZ空間と同一であること: 本実験では未導出。
- ABCDの六関係波から観測可能な三軸を一意選択する機構: 本実験では未解決。

## 出力

- `ab_abc_abcd_complete_pair_relation_network_preliminary_result_v1.json`
- `ab_abc_abcd_complete_pair_relation_network_trial_summary_v1.csv`
- `ab_abc_abcd_complete_pair_relation_network_body_summary_v1.csv`
- `ab_abc_abcd_complete_pair_relation_network_selected_series_v1.csv`
- `complete_pair_relation_wave_count_v1.png`
- `generator_plane_normal_structure_v1.png`
- `ABC_three_physical_relation_waves_v1.png`
- `ABC_relation_wave_conservation_v1.png`
- `ABC_one_plane_one_normal_conservation_v1.png`
