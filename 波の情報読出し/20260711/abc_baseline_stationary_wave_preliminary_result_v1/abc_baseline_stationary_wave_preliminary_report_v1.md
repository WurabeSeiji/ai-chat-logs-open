# ABCベースライン定常波C予備実験 v1

## 目的

本予備実験は、サブセル位置位相加速度読出し本実験の前に、ゼロであるべき条件、読出し器影響、C定常波安定性、閉鎖再選別の制御挙動を確認する。

本実験は重力的読出しの検出を主張しない。目的は、後続の Stage I に入る前の数値床と対照条件を確定することである。

## 統合判定

- phase0_no_coupling_valid: `True`
- phase1_readout_only_valid: `True`
- c_baseline_stability_valid: `True`
- closure_reselection_control_valid: `True`
- single_gauge_only_used: `False`
- preliminary_experiment_valid: `True`
- max_phase0_delta_chi_abs: `0.0`
- max_phase1_delta_chi_abs: `0.0`
- max_phase0_acceleration_abs: `0.0`
- max_phase1_acceleration_abs: `0.0`
- max_c_center_error: `0.0`
- max_closure_Q_raw_abs: `4.000453976226206e-06`
- max_closure_Q_closed_abs: `0.0`
- max_seed_preservation_error: `1.1964902862783528e-17`
- gauge_count: `9`

## 予備フェーズ

| phase | purpose | valid | max delta chi | max acceleration |
|---|---|---|---:|---:|
| phase0_no_coupling | 無結合条件で実装ドリフトがないことを確認する | `True` | 0.0000000000000000e+00 | 0.0000000000000000e+00 |
| phase1_readout_only | Cで読むだけでは位置位相加速度が作られないことを確認する | `True` | 0.0000000000000000e+00 | 0.0000000000000000e+00 |

## C定常波安定性

| k_max | max center error | valid |
|---:|---:|---|
| 1 | 0.0000000000000000e+00 | `True` |
| 3 | 0.0000000000000000e+00 | `True` |
| 5 | 0.0000000000000000e+00 | `True` |
| 9 | 0.0000000000000000e+00 | `True` |
| 17 | 0.0000000000000000e+00 | `True` |
| 33 | 0.0000000000000000e+00 | `True` |
| 65 | 0.0000000000000000e+00 | `True` |
| 129 | 0.0000000000000000e+00 | `True` |

## 閉鎖再選別制御

| case | target | Q_raw | Q_closed | delta chi | seed error | valid |
|---|---|---:|---:|---:|---:|---|
| amplitude_only_A | A | 4.0000040001156750e-06 | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 0.0000000000000000e+00 | `True` |
| phase_seed_A | A | 4.0004539762262061e-06 | 0.0000000000000000e+00 | 3.0000000011964900e-08 | 1.1964902862783528e-17 | `True` |
| phase_seed_B | B | 3.1256233744244042e-06 | 0.0000000000000000e+00 | -1.9999999989472886e-08 | 1.0527114065656107e-17 | `True` |

## 解釈

- Phase 0 は無結合条件で、実装ドリフトが数値床に収まるかを確認する。
- Phase 1 は読出し器影響のみで、Cで読むだけでは位置位相加速度が作られないことを確認する。
- C定常波安定性は、倍音上限を変えても左右対称Cが方向ドリフトを持たないことを確認する。
- 閉鎖再選別制御では、振幅だけの閉鎖破れは位置位相を作らず、位相シードを入れた場合のみ Q を閉じた後にも δχ が残ることを確認する。

## 出力

| 種類 | ファイル |
|---|---|
| JSON | `abc_baseline_stationary_wave_preliminary_result_v1.json` |
| timeline CSV | `abc_baseline_stationary_wave_preliminary_timeline_v1.csv` |
| gauge CSV | `abc_baseline_stationary_wave_preliminary_gauge_rows_v1.csv` |
| C stability CSV | `abc_baseline_stationary_wave_preliminary_c_stability_v1.csv` |
| closure CSV | `abc_baseline_stationary_wave_preliminary_closure_controls_v1.csv` |
