# ABCベースライン定常波C Stage III 別読出し照合予備実験 v1

## 目的

準安定傾斜から読まれる加速度候補が、単一の観測方法に依存していないかを検査する。

同じ遷移プロトコルに対して、位置位相二階差分、`p_read` 相当の速度差分、C傾斜による G1 予測、Cからの積分再構成を照合する。

## 統合判定

- protocol_count: `9`
- single_gauge_only_used: `False`
- all_chi_p_readouts_match: `True`
- all_chi_g1_readouts_match: `True`
- all_p_g1_readouts_match: `True`
- all_active_signs_consistent: `True`
- all_integrated_errors_near_zero: `True`
- all_R_weighted_integrated_balances_near_zero: `True`
- max_chi_minus_g1_error_abs: `1.0339757656912846e-23`
- max_integrated_error_abs: `0.0`
- stage3_cross_readout_preliminary_valid: `True`

## プロトコル別サマリー

| protocol | kind | active | max chi | max p | max G1 | chi-p err | chi-G1 err | signs | integrated balance |
|---|---|---:|---:|---:|---:|---:|---:|---|---:|
| ramp_duration_1 | ramp | 2 | 2.3170731707317074e-09 | 2.3170731707317074e-09 | 2.3170731707317074e-09 | 1.0339757656912846e-23 | 1.0339757656912846e-23 | `True` | 2.6469779601696886e-23 |
| ramp_duration_2 | ramp | 4 | 1.1585365853658537e-09 | 1.1585365853658537e-09 | 1.1585365853658537e-09 | 1.0339757656912846e-23 | 1.0339757656912846e-23 | `True` | 1.3234889800848443e-23 |
| ramp_duration_4 | ramp | 8 | 5.7926829268292727e-10 | 5.7926829268292706e-10 | 5.7926829268292696e-10 | 1.0339757656912846e-23 | 1.0339757656912846e-23 | `True` | 1.3234889800848443e-23 |
| ramp_duration_8 | ramp | 16 | 2.8963414634146374e-10 | 2.8963414634146374e-10 | 2.8963414634146348e-10 | 6.6174449004242214e-24 | 6.6174449004242214e-24 | `True` | 0.0000000000000000e+00 |
| ramp_duration_16 | ramp | 32 | 1.4481707317073290e-10 | 1.4481707317073208e-10 | 1.4481707317073187e-10 | 6.6174449004242214e-24 | 6.6174449004242214e-24 | `True` | 1.3234889800848443e-23 |
| ramp_duration_32 | ramp | 64 | 7.2408536585368520e-11 | 7.2408536585366258e-11 | 7.2408536585366064e-11 | 6.6174449004242214e-24 | 6.6174449004242214e-24 | `True` | 2.6469779601696886e-23 |
| restart_stable_final_C | stable_restart | 0 | 1.0339757656912846e-23 | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 1.0339757656912846e-23 | 1.0339757656912846e-23 | `True` | 3.9704669402545328e-23 |
| down_quench_from_final_C | down_quench | 2 | 2.3170731707317074e-09 | 2.3170731707317074e-09 | 2.3170731707317074e-09 | 8.2718061255302767e-25 | 8.2718061255302767e-25 | `True` | 0.0000000000000000e+00 |
| overshoot_relax_to_final_C | overshoot_relax | 32 | 8.1097560975609754e-10 | 8.1097560975609764e-10 | 8.1097560975609754e-10 | 6.6174449004242214e-24 | 6.6174449004242214e-24 | `True` | 2.6469779601696886e-23 |

## 解釈

- 位置位相二階差分と `p_read` 差分が一致する場合、加速度候補は位置読出しだけの人工物ではない。
- それらが C傾斜による G1 予測とも一致する場合、時間位相側の補償読出しとも整合する。
- Cからの積分再構成が最終変位と一致し、R重み付き積分バランスが保たれる場合、局所差分だけの偶然ではない。
- 本実験は標準重力の導出ではなく、準安定傾斜候補の観測方法依存性を下げるための Stage III 予備照合である。

## 出力

| 種類 | ファイル |
|---|---|
| JSON | `abc_baseline_stationary_wave_stage3_cross_readout_preliminary_result_v1.json` |
| protocol CSV | `abc_baseline_stationary_wave_stage3_cross_readout_cases_v1.csv` |
| rows CSV | `abc_baseline_stationary_wave_stage3_cross_readout_rows_v1.csv` |
| integrated CSV | `abc_baseline_stationary_wave_stage3_cross_readout_integrated_v1.csv` |
