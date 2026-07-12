# AB二体閉鎖位相系における調和読出しとc=1面積スイープ予備実験総括 v1

**日付:** 2026-07-12  
**著者:** 木原 範昭  
**位置づけ:** 波の情報読出しシリーズ・AB二体予備実験群の総括  

---

## 0. 結論

本総括は、AB二体閉鎖位相系に対して行った次の予備実験群をまとめる。

1. 一角度円周位相調和読出し
2. 一角度円周位相調和読出しパラメータスイープ
3. `c=1` 内部較正と `chi-tau` 面積スイープ
4. `c=1` 内部較正パラメータスイープ
5. `chi-tau` 面積逆数補償診断
6. `chi-tau` ネイティブ逆面積拡張スイープ

結論は次である。

```text
前回のABC多ゲージ干渉読出し実験では、
背景空間を先に置かないABC閉鎖位相系という最小モデルにおいて、
測定機自身も系の中の複素位相波として定義し、
干渉読出しから質量様、運動量様、エネルギー様に見える保存量を構成できることを示した。

本実験群では、さらに踏み込み、
AB二体の閉鎖位相関係だけから、加速度様に見える調和読出しが可能であることを確認した。

ただし、その加速度様読出しが位置位相差、すなわち距離に対して
逆比例または逆二乗比例で変化することは、AB二体系だけでは判別できなかった。
```

本実験群の中心成果は、単なる調和振動の表示ではない。

外部の標準力、標準ばね式、重力式、クーロン式を置かず、
`f_A`, `f_B` という個体別力も置かず、
ラベルなし二体関係 `f_AB` のみを使って、加速度様に見える読出しを作れた点にある。

一方で、本実験群は次の限界も明確にした。

AB二体系だけでは、測定機自身が測定対象の上に乗っている。

そのため、相対距離の変化量は読めても、その変化が

```text
L_AB
1 / L_AB
1 / L_AB^2
```

のどれに従うかを独立に判定するゲージがない。

加速度様読出しによって間隔が変化しても、その間隔を測るゲージ自身も同じ閉鎖位相系の変化を受けるためである。

さらに、空間方向の位置位相差だけでなく、時間軸方向の位相差も導入し、内部 `c=1` 相当の較正を置いて時間発展に近い効果を再現しても、結論は変わらなかった。

したがって、距離指数の判別は、別途独立した計量ゲージを持つABC三体実験へ送るべき課題として残った。

---

## 1. 実験一覧

| No. | 実験 | 目的 | 主結果 |
|---:|---|---|---|
| 1 | 一角度調和読出し | `D_AB`, `V_AB`, `f_AB` がラベルなしで読めるか | 読める |
| 2 | 一角度パラメータスイープ | 初期偏差、周期、読出し漏れへの頑健性 | 全ケースで調和振動を検出 |
| 3 | `c=1` 面積スイープ | `s` と `tau_read` を分離し、`chi-tau` 面が立つか | 独立 `tau` で面積検出 |
| 4 | `c=1` パラメータスイープ | `c=1` が面積成立の十分条件か | 十分条件ではない |
| 5 | 逆面積補償診断 | `1/A_chi_tau` が native に出るか | 未検出 |
| 6 | native 逆面積拡張スイープ | 広い条件で `alpha≈2` が出るか | native では未検出 |

---

## 2. 一角度円周位相調和読出し

### 2.1 実験の意味

AB二体だけを閉じた位相系として読み、観測機 `C` を置かない条件で、

```text
f_A
f_B
```

を独立力として導入せず、

```text
f_AB
```

のみを関係性補償として読むことを検査した。

この実験では、標準重力式、標準クーロン式、標準ばね式は使わない。

### 2.2 主結果

主要な数値結果は次である。

| 量 | 値 |
|---|---:|
| `case_count` | `32` |
| `observer_C_used` | `False` |
| `f_A_or_f_B_used` | `False` |
| `max_Q_closed_abs` | `0.0` |
| `oscillation_detected_all_cases` | `True` |
| `label_free_protocol_degenerate_all_cases` | `True` |
| `readout_decay_monotonic_all_cases` | `True` |
| `readout_off_decay_max_abs` | `6.49e-17` |
| `readout_strong_decay_min_abs` | `4.0004e-4` |

成立したことは次である。

```text
Protocol F/B は内部表示として異なるが、
D_AB と V_AB では完全に縮退した。
```

これは、ラベルなし読出しとしては、反跳型と通過型を区別できないことを意味する。

また、

```text
readout_off では減衰なし
readout が強いほど包絡減衰が増加
```

が確認された。

これは、読出し波が閉鎖位相系から情報を外へ読むとき、系の包絡に影響しうることを示す反証テストである。

### 2.3 図

| AB二体幾何 | 主要観測総括 |
|---|---|
| <img src="AB二体問題の図化_fAB_v1.png" width="420"> | <img src="ab_two_body_one_angle_harmonic_readout_observation_figures_v1/ab_two_body_one_angle_observation_summary_v1.png" width="520"> |

| 調和状態 | 読出しリーク応答 |
|---|---|
| <img src="ab_two_body_one_angle_harmonic_readout_observation_figures_v1/ab_two_body_one_angle_harmonic_state_v1.png" width="520"> | <img src="ab_two_body_one_angle_harmonic_readout_observation_figures_v1/ab_two_body_one_angle_readout_leak_response_v1.png" width="520"> |

| Protocol 縮退 | 位相差スケーリング |
|---|---|
| <img src="ab_two_body_one_angle_harmonic_readout_observation_figures_v1/ab_two_body_one_angle_protocol_degeneracy_v1.png" width="520"> | <img src="ab_two_body_one_angle_harmonic_readout_observation_figures_v1/ab_two_body_one_angle_phase_difference_scaling_v1.png" width="520"> |

---

## 3. 一角度パラメータスイープ

### 3.1 主結果

| 量 | 値 |
|---|---:|
| `sweep_configuration_count` | `210` |
| `case_summary_count` | `420` |
| `period_count` | `5` |
| `deviation_count` | `7` |
| `leak_count` | `6` |
| `max_Q_closed_abs` | `0.0` |
| `label_free_protocol_degenerate_all_cases` | `True` |
| `oscillation_detected_all_cases` | `True` |
| `readout_off_decay_max_abs` | `8.26e-17` |
| `decay_abs_monotonic_by_leak_all_grids` | `True` |
| `max_normalized_f_AB_projection_error` | `2.58e-4` |

本スイープにより、一角度版は弱読出しでは安定である一方、強い読出し波は `f_AB` 射影を歪めることが確認された。

| 漏れ量と射影誤差 | 漏れ量と減衰 |
|---|---|
| <img src="ab_two_body_one_angle_harmonic_readout_parameter_sweep_preliminary_result_v1/ab_two_body_one_angle_parameter_sweep_leak_f_error_v1.png" width="520"> | <img src="ab_two_body_one_angle_harmonic_readout_parameter_sweep_preliminary_result_v1/ab_two_body_one_angle_parameter_sweep_leak_decay_v1.png" width="520"> |

| 周期スイープ | 選択時系列 |
|---|---|
| <img src="ab_two_body_one_angle_harmonic_readout_parameter_sweep_preliminary_result_v1/ab_two_body_one_angle_parameter_sweep_period_v1.png" width="520"> | <img src="ab_two_body_one_angle_harmonic_readout_parameter_sweep_preliminary_result_v1/ab_two_body_one_angle_parameter_sweep_selected_series_v1.png" width="520"> |

---

## 4. c=1内部較正とchi-tau面積スイープ

### 4.1 実験の意味

一角度調和読出しでは、位相差に応じた距離減衰は読めなかった。

そこで、実験ステップ `s` を時間そのものとして扱わず、独立な時間位相候補

```text
tau_read
```

を導入した。

目的は、

```text
chi_read と tau_read が独立な面を作るか
```

を検査することである。

### 4.2 主結果

| 量 | 値 |
|---|---:|
| `case_summary_count` | `288` |
| `power_candidate_count` | `48` |
| `max_Q_closed_abs` | `0.0` |
| `disabled_max_area` | `0.0` |
| `locked_max_area` | `7.11e-15` |
| `independent_min_area` | `0.0024367633602631385` |
| `c1_readout_off_max_epsilon_c_abs` | `2.33e-15` |
| `c1_area_sweep_detected_all_cases` | `True` |
| `tau_is_step_used_any` | `False` |
| `external_c_used_any` | `False` |

成立したことは次である。

```text
tau disabled は面積を作らない。
tau locked は面積を作らない。
tau independent は chi-tau 面積を作る。
```

したがって、時間位相を独立に読めば `chi-tau` 面は立つ。

しかし、それだけでは距離減衰は出ない。

### 4.3 図

| `chi-tau` 面 | 面積スイープ |
|---|---|
| <img src="ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1/ab_two_body_c1_internal_calibration_chi_tau_surface_v1.png" width="520"> | <img src="ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1/ab_two_body_c1_internal_calibration_area_sweep_v1.png" width="520"> |

| `c=1` 較正誤差 | 冪候補 |
|---|---|
| <img src="ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1/ab_two_body_c1_internal_calibration_error_v1.png" width="520"> | <img src="ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1/ab_two_body_c1_internal_calibration_power_candidate_v1.png" width="520"> |

| 読出しリーク |
|---|
| <img src="ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1/ab_two_body_c1_internal_calibration_readout_leak_v1.png" width="620"> |

---

## 5. c=1内部較正パラメータスイープ

### 5.1 主結果

| 量 | 値 |
|---|---:|
| `sweep_case_count` | `756` |
| `readout_off_case_count` | `252` |
| `rank2_readout_off_count` | `246` |
| `c1_surface_like_readout_off_count` | `13` |
| `c1_locked_like_readout_off_count` | `8` |
| `min_c_error_readout_off` | `1.11e-15` |
| `max_area_readout_off` | `0.19685536479742288` |

この結果により、

```text
c=1
rank_chi_tau = 2
A_chi_tau != 0
```

の三条件を同時に要求する必要が明確になった。

`c=1` に近いだけでは、`chi-tau` 面が成立したとは言えない。

| c 誤差 heatmap | 面積 heatmap |
|---|---|
| <img src="ab_two_body_c1_internal_calibration_parameter_sweep_preliminary_result_v1/ab_two_body_c1_parameter_sweep_c_error_heatmap_v1.png" width="520"> | <img src="ab_two_body_c1_internal_calibration_parameter_sweep_preliminary_result_v1/ab_two_body_c1_parameter_sweep_area_heatmap_v1.png" width="520"> |

| 位相差応答 | 読出しリーク |
|---|---|
| <img src="ab_two_body_c1_internal_calibration_parameter_sweep_preliminary_result_v1/ab_two_body_c1_parameter_sweep_phase_response_v1.png" width="520"> | <img src="ab_two_body_c1_internal_calibration_parameter_sweep_preliminary_result_v1/ab_two_body_c1_parameter_sweep_readout_leak_v1.png" width="520"> |

---

## 6. 逆面積補償診断

### 6.1 目的

`chi-tau` 面積が成立した後、

```text
1 / A_chi_tau
```

が閉鎖補償側に自然に残るかを検査した。

重要なのは、`1/A_chi_tau` を後処理で作ることではない。

それは構成済み対照であり、native readout とは区別する。

### 6.2 主結果

| 量 | 値 |
|---|---:|
| `diagnostic_case_count` | `288` |
| `area_valid_case_count` | `144` |
| `fit_count` | `576` |
| `native_fit_count` | `132` |
| `derived_fit_count` | `84` |
| `native_positive2_count` | `0` |
| `constructed_reciprocal_positive2_count` | `2` |

結果は明確である。

```text
1 / A_chi_tau を作れば alpha≈2 になる。
しかし native readout には alpha≈2 は出ていない。
```

| 構成済み逆面積対照 | native 候補 |
|---|---|
| <img src="ab_two_body_chi_tau_inverse_area_compensation_diagnostic_preliminary_result_v1/ab_two_body_chi_tau_inverse_area_constructed_control_v1.png" width="520"> | <img src="ab_two_body_chi_tau_inverse_area_compensation_diagnostic_preliminary_result_v1/ab_two_body_chi_tau_inverse_area_native_candidates_v1.png" width="520"> |

| alpha 比較 |
|---|
| <img src="ab_two_body_chi_tau_inverse_area_compensation_diagnostic_preliminary_result_v1/ab_two_body_chi_tau_inverse_area_alpha_comparison_v1.png" width="620"> |

---

## 7. native逆面積拡張スイープ

### 7.1 主結果

| 量 | 値 |
|---|---:|
| `sweep_case_count` | `1323` |
| `area_valid_case_count` | `1260` |
| `c1_surface_like_case_count` | `126` |
| `fit_count` | `4158` |
| `native_fit_count` | `1056` |
| `native_positive2_count` | `0` |
| `c1_native_positive2_count` | `0` |
| `constructed_reciprocal_positive2_count` | `198` |

広いパラメータ範囲でも、

```text
native inverse-area scaling は未検出
```

であった。

一方で、構成済み対照では `alpha≈2` が出る。

したがって、今回の境界は次である。

```text
逆二乗は作れる。
しかし、まだ読まれていない。
```

| 参照曲線 | alpha scan |
|---|---|
| <img src="ab_two_body_chi_tau_native_inverse_area_extended_sweep_preliminary_result_v1/ab_two_body_native_inverse_area_extended_reference_curves_v1.png" width="520"> | <img src="ab_two_body_chi_tau_native_inverse_area_extended_sweep_preliminary_result_v1/ab_two_body_native_inverse_area_extended_alpha_scan_v1.png" width="520"> |

---

## 8. 総合判定

### 8.1 成立したこと

| 項目 | 判定 |
|---|---|
| AB二体で `f_AB` を単一関係性として読む | 成立 |
| `f_A`, `f_B` を置かずに調和振動を読む | 成立 |
| Protocol F/B のラベルなし縮退を確認 | 成立 |
| 読出し波停止で減衰が消えることを確認 | 成立 |
| 独立 `tau_read` により `chi-tau` 面積を読む | 成立 |
| `c=1` が必要条件候補であることを確認 | 成立 |
| 後処理 `1/A_chi_tau` が `alpha≈2` を示す | 成立 |

### 8.2 成立していないこと

| 項目 | 判定 |
|---|---|
| AB二体だけで距離指数を独立計量する | 未成立 |
| `c=1` だけで時間位相面が成立する | 未成立 |
| native な逆比例型読出し | 未検出 |
| native な逆二乗型読出し | 未検出 |
| `chi-tau` 面積が自動で `1/A` 補償へ変換される | 未検出 |
| 標準重力対応 | 未成立 |

---

## 9. 解釈

AB二体系の最大の成果は、次の二点である。

```text
1. AB二体の関係 f_AB だけで、加速度様に見える調和読出しを構成できる。
2. AB二体だけでは、距離指数を読む独立ゲージが足りない。
```

第一の点は積極的な成果である。

前回のABC多ゲージ干渉読出しでは、質量様、運動量様、エネルギー様の読出しを干渉相関から再構成した。

本実験では、それに続いて、加速度様読出しがAB二体の調和関係から構成できることを示した。

これは、外部から標準力を入れた結果ではない。

二体関係 `f_AB` の閉鎖読出しとして現れた調和変位である。

第二の点は境界条件である。

AB二体の中では、A と B は互いを通じてしか読めない。

そのため、相対距離の変化は読めても、その変化率を

```text
比例
逆比例
逆二乗
```

のどれとして計量すべきかを固定する独立ゲージがない。

これは無名性の失敗ではなく、むしろ無名性の帰結である。

読めないものを、後処理で読めたことにしてはならない。

`chi-tau` 面積スイープを導入した実験では、独立時間位相を置くと面積読出しは成立した。

しかし、それは距離指数の読出しを自動的に与えなかった。

また、後処理で `1/A_chi_tau` を作れば逆二乗型は当然に現れるが、これは native readout ではない。

したがって、本総括の判定は次である。

```text
加速度様読出しは確認された。
距離指数読出しは、AB二体だけでは判別不能である。
```

---

## 10. 次の実験への接続

この総括から導かれる次の段階は、ABC三体系である。

ただし、第三波 `C` は外部の絶対観測器ではない。

`C` 自身も閉鎖位相系の一部であり、

```text
f_AC
f_BC
f_ABC
```

を生成する。

したがって、ABC実験では、

```text
C を独立計量ゲージとして使えるか
C が AB 主読出しを汚染しないか
f_ABC を代表時間として使えるか
f_AB, f_AC, f_BC を円周方向候補として分離できるか
```

を検査する必要がある。

---

# 参考文献

## 自己引用

1. 木原範昭「無名等振幅複合波モデル基本公理系 v4」Version DOI: `10.5281/zenodo.21316620`, Concept DOI: `10.5281/zenodo.21315735`, 2026.
2. 木原範昭「ABC閉鎖位相系における多ゲージ干渉読出し保存量の構成実験」Version DOI: `10.5281/zenodo.21308050`, Concept DOI: `10.5281/zenodo.21308049`, 2026.
3. 木原範昭, [AB二体閉鎖位相系におけるラベルなし二弧相対位相と調和読出しに関する定義補足.md](AB二体閉鎖位相系におけるラベルなし二弧相対位相と調和読出しに関する定義補足.md), 2026.
4. 木原範昭, [AB二体閉鎖位相系における一角度円周位相調和読出し実験仕様書 v1.md](AB二体閉鎖位相系における一角度円周位相調和読出し実験仕様書%20v1.md), 2026.
5. 木原範昭, [AB二体閉鎖位相系におけるc=1内部較正と空間位相・時間位相面積スイープ読出し実験仕様書 v1.md](AB二体閉鎖位相系におけるc=1内部較正と空間位相・時間位相面積スイープ読出し実験仕様書%20v1.md), 2026.
6. 木原範昭, [閉鎖複素位相波における自己項の内部閉鎖とN体外部読出し分離に関する定義補足.md](閉鎖複素位相波における自己項の内部閉鎖とN体外部読出し分離に関する定義補足.md), 2026.

## 外部参考文献

外部参考文献は、本稿の導出根拠ではなく、波動読出し、位相、観測確率に関する標準的背景を示すために最小限に用いる。

7. Max Born, "Zur Quantenmechanik der Stossvorgaenge", *Zeitschrift fuer Physik* 37, 863-867, 1926. DOI: `10.1007/BF01397477`.
8. Y. Aharonov and D. Bohm, "Significance of electromagnetic potentials in the quantum theory", *Physical Review* 115, 485-491, 1959. DOI: `10.1103/PhysRev.115.485`.
