# ABCベースライン定常波C Stage II frozen fanout 予備実験 v1

## 目的

同一状態スナップショットから、相対位相、位置位相、C残渣、反力候補、beat、閉鎖残差を同時に読む。

目的は、C媒介接近候補、R2反力、R3反力が混在したとき、合成が単純和に近いか、反力差し引き後の残差として C媒介接近が残るかを確認することである。

## 統合判定

- case_count: `8`
- composition_check_count: `5`
- all_composition_checks_valid: `True`
- all_R_balances_valid: `True`
- all_gauge_stability_valid: `True`
- all_Q_closed_valid: `True`
- single_gauge_only_used: `False`
- stage2_frozen_fanout_preliminary_valid: `True`
- C_mediated_distance_change: `-1.3414465233063666e-07`
- R2_distance_change: `0.00014981441759243097`
- R3_distance_change: `9.363072540669037e-05`
- combined_distance_change: `0.00024332848252417438`
- residual_after_R2_R3: `-1.166604749469613e-07`
- reaction_dominance_ratio_abs: `1814.7957355697142`
- beat_alias_reverse_all_cases: `True`

## ケース別サマリー

| case | C | R2 | R3 | distance change | R balance | Q_raw max | beat reverse |
|---|---|---|---|---:|---:|---:|---|
| C_mediated_only_persistent | `True` | `False` | `False` | -1.3414465233063666e-07 | -7.2672779896458799e-20 | 3.5873841187439844e-09 | `True` |
| C_mediated_only_reset | `True` | `False` | `False` | -1.4031966755201353e-08 | -7.2686014786259648e-20 | 2.9233264549036767e-10 | `True` |
| R2_only | `False` | `True` | `False` | 1.4981441759243097e-04 | -2.7105054312137611e-20 | 3.8117932330638805e-09 | `True` |
| R3_only | `False` | `False` | `True` | 9.3630725406690374e-05 | -1.0164395367051604e-19 | 2.3821850484109543e-09 | `True` |
| C_plus_R2 | `True` | `True` | `False` | 1.4968024876083064e-04 | -6.7762635780344027e-20 | 7.3998776036025568e-09 | `True` |
| C_plus_R3 | `True` | `False` | `True` | 9.3496565642503526e-05 | -7.4538899358378430e-20 | 5.9700068183052459e-09 | `True` |
| C_plus_R2_plus_R3 | `True` | `True` | `True` | 2.4332848252417438e-04 | -8.1315162936412833e-20 | 9.7834908868101551e-09 | `True` |
| C_plus_R2_plus_R3_mirrored | `True` | `True` | `True` | 2.4332848252417438e-04 | -5.4210108624275222e-20 | 9.7834908868101551e-09 | `True` |

## 合成・残差検査

| metric | value | tolerance | valid |
|---|---:|---:|---|
| C_plus_R2_minus_sum | 2.4179269697555128e-11 | 4.9999999999999998e-08 | `True` |
| C_plus_R3_minus_sum | 1.5111856210836550e-11 | 4.9999999999999998e-08 | `True` |
| C_plus_R2_plus_R3_minus_sum | 1.7484177383675359e-08 | 4.9999999999999998e-08 | `True` |
| residual_after_R2_R3_minus_C | 1.7484177383675359e-08 | 4.9999999999999998e-08 | `True` |
| mirrored_combined_distance_match | 0.0000000000000000e+00 | 4.9999999999999998e-08 | `True` |

## 解釈

- `combined ≈ C + R2 + R3` が成立する場合、Stage I 候補の粗い線形分離が可能である。
- `combined - R2 - R3 ≈ C` が成立し、かつ C 残差が負の距離変化を保つ場合、反力候補を差し引いた後にも C媒介接近候補が残る。
- ただし本予備実験では R2/R3 が C媒介接近より大きい。したがって、観測方法依存性を避けるには、G1 へ進む前に同時多読出しの記録を保持する必要がある。
- beat 系列は負のサブ巻数残差により逆向きに読まれる。これは G3 候補として、距離位相残差とは別に保持する。

## 出力

| 種類 | ファイル |
|---|---|
| JSON | `abc_baseline_stationary_wave_stage2_frozen_fanout_preliminary_result_v1.json` |
| case CSV | `abc_baseline_stationary_wave_stage2_frozen_fanout_cases_v1.csv` |
| fanout CSV | `abc_baseline_stationary_wave_stage2_frozen_fanout_rows_v1.csv` |
