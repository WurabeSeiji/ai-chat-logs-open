# AB二体閉鎖位相系における一角度円周位相調和読出し予備実験検証メモ v1

## 目的

観測機 C を置かない AB 二体系で、ラベルなし二弧相対位相 `D_AB`、対称偏差 `V_AB`、AB 合成補償 `f_AB` を読む予備実験を行った。

本実験では、標準重力式、標準クーロン式、標準ばね式を使わない。

一角度の閉鎖補助平面における複素回転

```text
z(s+1) = lambda * exp(i Omega) * z(s)
```

を用い、`Protocol F/B` と読出し波条件の違いがラベルなし読出しへどう現れるかを検査した。

## 統合判定

| 量 | 値 |
|---|---:|
| `case_count` | `32` |
| `initial_case_count` | `4` |
| `protocol_count` | `2` |
| `readout_mode_count` | `4` |
| `observer_C_used` | `False` |
| `single_gauge_only_used` | `False` |
| `absolute_background_axis_used` | `False` |
| `f_A_or_f_B_used` | `False` |
| `standard_force_law_used` | `False` |
| `max_Q_closed_abs` | `0.0` |
| `max_Q_raw_abs` | `2.4369393582936685e-05` |
| `max_f_AB_projection_consistency_error` | `9.089175649484855e-06` |
| `max_f_AB_projection_consistency_error_nonstrong` | `2.2803188765678194e-06` |
| `max_f_AB_projection_consistency_error_strong` | `9.089175649484855e-06` |
| `f_AB_projection_consistent_nonstrong_modes` | `True` |
| `strong_readout_perturbs_f_AB_projection` | `True` |
| `max_D_AB_near_protocol_diff` | `0.0` |
| `max_V_AB_protocol_diff` | `0.0` |
| `oscillation_detected_all_cases` | `True` |
| `label_free_protocol_degenerate_all_cases` | `True` |
| `readout_decay_monotonic_all_cases` | `True` |
| `readout_off_decay_max_abs` | `6.490828430734141e-17` |
| `readout_strong_decay_min_abs` | `0.0004000400053341305` |
| `ab_one_angle_harmonic_readout_preliminary_valid` | `True` |

## ケース別サマリー

| case | protocol | readout | decay_rate | envelope final/initial | sign changes | max f error |
|---|---|---|---:|---:|---:|---:|
| near_pi_02deg | Protocol_F | readout_off | -6.4908284307341414e-17 | 9.9999999999995803e-01 | 15 | 6.3696877633523385e-18 |
| near_pi_02deg | Protocol_F | readout_weak | -2.0000100000660538e-05 | 9.8570311315135206e-01 | 15 | 4.5649268955116595e-08 |
| near_pi_02deg | Protocol_F | readout_normal | -1.0000250008349211e-04 | 9.3052922080116551e-01 | 15 | 2.2803188765711786e-07 |
| near_pi_02deg | Protocol_F | readout_strong | -4.0004000533413277e-04 | 7.4973999653669299e-01 | 15 | 9.0891756494848099e-07 |
| near_pi_02deg | Protocol_B | readout_off | -6.4908284307341414e-17 | 9.9999999999995803e-01 | 15 | 6.3696877633523385e-18 |
| near_pi_02deg | Protocol_B | readout_weak | -2.0000100000660538e-05 | 9.8570311315135206e-01 | 15 | 4.5649268955116595e-08 |
| near_pi_02deg | Protocol_B | readout_normal | -1.0000250008349211e-04 | 9.3052922080116551e-01 | 15 | 2.2803188765711786e-07 |
| near_pi_02deg | Protocol_B | readout_strong | -4.0004000533413277e-04 | 7.4973999653669299e-01 | 15 | 9.0891756494848099e-07 |
| near_pi_05deg | Protocol_F | readout_off | -5.4726592651287856e-17 | 9.9999999999996181e-01 | 15 | 1.3552527156068805e-17 |
| near_pi_05deg | Protocol_F | readout_weak | -2.0000100000659265e-05 | 9.8570311315135062e-01 | 15 | 1.1412317238787395e-07 |
| near_pi_05deg | Protocol_F | readout_normal | -1.0000250008348830e-04 | 9.3052922080116796e-01 | 15 | 5.7007971914195485e-07 |
| near_pi_05deg | Protocol_F | readout_strong | -4.0004000533413336e-04 | 7.4973999653669221e-01 | 15 | 2.2722939123712136e-06 |
| near_pi_05deg | Protocol_B | readout_off | -5.4726592651287856e-17 | 9.9999999999996181e-01 | 15 | 1.3552527156068805e-17 |
| near_pi_05deg | Protocol_B | readout_weak | -2.0000100000659265e-05 | 9.8570311315135062e-01 | 15 | 1.1412317238787395e-07 |
| near_pi_05deg | Protocol_B | readout_normal | -1.0000250008348830e-04 | 9.3052922080116796e-01 | 15 | 5.7007971914195485e-07 |
| near_pi_05deg | Protocol_B | readout_strong | -4.0004000533413336e-04 | 7.4973999653669221e-01 | 15 | 2.2722939123712136e-06 |
| near_pi_10deg | Protocol_F | readout_off | -5.5999304108294559e-17 | 9.9999999999996181e-01 | 15 | 2.7105054312137611e-17 |
| near_pi_10deg | Protocol_F | readout_weak | -2.0000100000659265e-05 | 9.8570311315135062e-01 | 15 | 2.2824634477574789e-07 |
| near_pi_10deg | Protocol_F | readout_normal | -1.0000250008349021e-04 | 9.3052922080116796e-01 | 15 | 1.1401594382839097e-06 |
| near_pi_10deg | Protocol_F | readout_strong | -4.0004000533413081e-04 | 7.4973999653669221e-01 | 15 | 4.5445878247424273e-06 |
| near_pi_10deg | Protocol_B | readout_off | -5.5999304108294559e-17 | 9.9999999999996181e-01 | 15 | 2.7105054312137611e-17 |
| near_pi_10deg | Protocol_B | readout_weak | -2.0000100000659265e-05 | 9.8570311315135062e-01 | 15 | 2.2824634477574789e-07 |
| near_pi_10deg | Protocol_B | readout_normal | -1.0000250008349021e-04 | 9.3052922080116796e-01 | 15 | 1.1401594382839097e-06 |
| near_pi_10deg | Protocol_B | readout_strong | -4.0004000533413081e-04 | 7.4973999653669221e-01 | 15 | 4.5445878247424273e-06 |
| near_pi_20deg | Protocol_F | readout_off | -5.3772059058532839e-17 | 9.9999999999996181e-01 | 15 | 5.4210108624275222e-17 |
| near_pi_20deg | Protocol_F | readout_weak | -2.0000100000659265e-05 | 9.8570311315135062e-01 | 15 | 4.5649268955149579e-07 |
| near_pi_20deg | Protocol_F | readout_normal | -1.0000250008348606e-04 | 9.3052922080116796e-01 | 15 | 2.2803188765678194e-06 |
| near_pi_20deg | Protocol_F | readout_strong | -4.0004000533413049e-04 | 7.4973999653669221e-01 | 15 | 9.0891756494848546e-06 |
| near_pi_20deg | Protocol_B | readout_off | -5.3772059058532839e-17 | 9.9999999999996181e-01 | 15 | 5.4210108624275222e-17 |
| near_pi_20deg | Protocol_B | readout_weak | -2.0000100000659265e-05 | 9.8570311315135062e-01 | 15 | 4.5649268955149579e-07 |
| near_pi_20deg | Protocol_B | readout_normal | -1.0000250008348606e-04 | 9.3052922080116796e-01 | 15 | 2.2803188765678194e-06 |
| near_pi_20deg | Protocol_B | readout_strong | -4.0004000533413049e-04 | 7.4973999653669221e-01 | 15 | 9.0891756494848546e-06 |

## Protocol F/B 比較

| case | readout | max D diff | max V diff | display diff | label-free degenerate |
|---|---|---:|---:|---:|---|
| near_pi_02deg | readout_off | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 6.9813170079773043e-02 | True |
| near_pi_02deg | readout_weak | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 6.9779667631853157e-02 | True |
| near_pi_02deg | readout_normal | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 6.9645815193869906e-02 | True |
| near_pi_02deg | readout_strong | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 6.9146103979029219e-02 | True |
| near_pi_05deg | readout_off | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 1.7453292519943270e-01 | True |
| near_pi_05deg | readout_weak | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 1.7444916907963279e-01 | True |
| near_pi_05deg | readout_normal | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 1.7411453798467469e-01 | True |
| near_pi_05deg | readout_strong | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 1.7286525994757307e-01 | True |
| near_pi_10deg | readout_off | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 3.4906585039886540e-01 | True |
| near_pi_10deg | readout_weak | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 3.4889833815926558e-01 | True |
| near_pi_10deg | readout_normal | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 3.4822907596934938e-01 | True |
| near_pi_10deg | readout_strong | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 3.4573051989514614e-01 | True |
| near_pi_20deg | readout_off | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 6.9813170079773079e-01 | True |
| near_pi_20deg | readout_weak | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 6.9779667631853115e-01 | True |
| near_pi_20deg | readout_normal | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 6.9645815193869876e-01 | True |
| near_pi_20deg | readout_strong | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 6.9146103979029228e-01 | True |

## 読出し波停止反証テスト

| case | protocol | off | weak | normal | strong | monotonic |
|---|---|---:|---:|---:|---:|---|
| near_pi_02deg | Protocol_F | -6.4908284307341414e-17 | -2.0000100000660538e-05 | -1.0000250008349211e-04 | -4.0004000533413277e-04 | True |
| near_pi_02deg | Protocol_B | -6.4908284307341414e-17 | -2.0000100000660538e-05 | -1.0000250008349211e-04 | -4.0004000533413277e-04 | True |
| near_pi_05deg | Protocol_F | -5.4726592651287856e-17 | -2.0000100000659265e-05 | -1.0000250008348830e-04 | -4.0004000533413336e-04 | True |
| near_pi_05deg | Protocol_B | -5.4726592651287856e-17 | -2.0000100000659265e-05 | -1.0000250008348830e-04 | -4.0004000533413336e-04 | True |
| near_pi_10deg | Protocol_F | -5.5999304108294559e-17 | -2.0000100000659265e-05 | -1.0000250008349021e-04 | -4.0004000533413081e-04 | True |
| near_pi_10deg | Protocol_B | -5.5999304108294559e-17 | -2.0000100000659265e-05 | -1.0000250008349021e-04 | -4.0004000533413081e-04 | True |
| near_pi_20deg | Protocol_F | -5.3772059058532839e-17 | -2.0000100000659265e-05 | -1.0000250008348606e-04 | -4.0004000533413049e-04 | True |
| near_pi_20deg | Protocol_B | -5.3772059058532839e-17 | -2.0000100000659265e-05 | -1.0000250008348606e-04 | -4.0004000533413049e-04 | True |

## 解釈

- `Protocol F/B` は内部表示としては異なるが、`D_AB` と `V_AB` では縮退した。
- `readout_off` では包絡減衰が数値丸め範囲に留まり、読出し波を強くするほど減衰率が大きくなった。
- これは、読出し波が長期振幅へ影響しうるという反証テストの検出系として機能する。
- ただし、この予備実験は複素回転写像の検査であり、調和読出しが第一原理から自発的に出現したことの証明ではない。
- 逆二乗型は本実験の対象外であり、二角度以上の位置位相自由度拡張で検査する。

## 出力

| 種類 | ファイル |
|---|---|
| JSON | `ab_two_body_one_angle_harmonic_readout_preliminary_result_v1.json` |
| series CSV | `ab_two_body_one_angle_harmonic_readout_series_v1.csv` |
| case summary CSV | `ab_two_body_one_angle_harmonic_readout_case_summary_v1.csv` |
| protocol comparison CSV | `ab_two_body_one_angle_harmonic_readout_protocol_comparison_v1.csv` |
| readout decay CSV | `ab_two_body_one_angle_harmonic_readout_readout_decay_v1.csv` |
| readout mode plot | `ab_two_body_one_angle_harmonic_readout_readout_mode_comparison_v1.png` |
| protocol plot | `ab_two_body_one_angle_harmonic_readout_protocol_comparison_v1.png` |
| decay plot | `ab_two_body_one_angle_harmonic_readout_envelope_decay_v1.png` |
