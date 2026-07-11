# ABCベースライン定常波C 遷移プロトコル予備実験 v1

## 目的

同じ最終 `C_memory` に到達する場合でも、quench, adiabatic ramp, stable restart, overshoot により加速度候補がどう変わるかを調べる。

準安定傾斜起源なら、加速度候補は最終 `C_memory` の大きさではなく、`ΔC_memory` の大きさと符号に従う。

## 統合判定

- protocol_count: `9`
- single_gauge_only_used: `False`
- final_C_memory_target: `3.8e-06`
- all_protocols_same_final_C_except_down_quench: `True`
- ramp_max_Ra_strictly_decreases_with_duration: `True`
- quench_to_slowest_max_Ra_ratio: `31.999999999998824`
- restart_stable_acceleration_near_zero: `True`
- down_quench_negative_only_A: `True`
- overshoot_has_positive_and_negative_A: `True`
- post_settle_acceleration_near_zero_all_protocols: `True`
- max_balance_error_abs: `1.0339757656912846e-23`
- transition_protocol_preliminary_valid: `True`

## プロトコル別サマリー

| protocol | kind | duration | max |ΔC| | max |R*a| | post max |R*a| | sign A | max error |
|---|---|---:|---:|---:|---:|---|---:|
| ramp_duration_1 | ramp | 1 | 3.8000000000000000e-06 | 2.3170731707317074e-09 | 1.0339757656912846e-23 | +1/-0 | 1.0339757656912846e-23 |
| ramp_duration_2 | ramp | 2 | 1.9000000000000000e-06 | 1.1585365853658537e-09 | 1.0339757656912846e-23 | +2/-0 | 1.0339757656912846e-23 |
| ramp_duration_4 | ramp | 4 | 9.5000000000000022e-07 | 5.7926829268292727e-10 | 1.0339757656912846e-23 | +4/-0 | 1.0339757656912846e-23 |
| ramp_duration_8 | ramp | 8 | 4.7500000000000011e-07 | 2.8963414634146374e-10 | 6.6174449004242214e-24 | +8/-0 | 6.6174449004242214e-24 |
| ramp_duration_16 | ramp | 16 | 2.3750000000000027e-07 | 1.4481707317073290e-10 | 6.6174449004242214e-24 | +16/-0 | 6.6174449004242214e-24 |
| ramp_duration_32 | ramp | 32 | 1.1875000000000035e-07 | 7.2408536585368520e-11 | 6.6174449004242214e-24 | +32/-0 | 6.6174449004242214e-24 |
| restart_stable_final_C | stable_restart | 0 | 0.0000000000000000e+00 | 1.0339757656912846e-23 | 1.0339757656912846e-23 | +0/-0 | 1.0339757656912846e-23 |
| down_quench_from_final_C | down_quench | 1 | 3.8000000000000000e-06 | 2.3170731707317074e-09 | 0.0000000000000000e+00 | +0/-1 | 8.2718061255302767e-25 |
| overshoot_relax_to_final_C | overshoot_relax | 16 | 1.3300000000000000e-06 | 8.1097560975609754e-10 | 6.6174449004242214e-24 | +4/-12 | 6.6174449004242214e-24 |

## 解釈

- 同じ最終 `C_memory` でも ramp が遅いほど `|R*a|` が小さくなるなら、加速度候補は定常レベルではなく遷移速度に支配される。
- stable restart で非ゼロ `C_memory` があるにもかかわらず加速度候補が消えるなら、安定状態そのものは力候補を生まない。
- down quench で符号が反転し、overshoot relax で正負の両符号が出るなら、候補はポテンシャルの絶対量ではなく、準安定傾斜の向きに従う。
- これは標準重力の導出ではなく、引力的読出し候補が定常場ではなく遷移プロトコルに依存するかを調べる予備試験である。

## 出力

| 種類 | ファイル |
|---|---|
| JSON | `abc_baseline_stationary_wave_transition_protocol_preliminary_result_v1.json` |
| protocol CSV | `abc_baseline_stationary_wave_transition_protocol_cases_v1.csv` |
| rows CSV | `abc_baseline_stationary_wave_transition_protocol_rows_v1.csv` |
