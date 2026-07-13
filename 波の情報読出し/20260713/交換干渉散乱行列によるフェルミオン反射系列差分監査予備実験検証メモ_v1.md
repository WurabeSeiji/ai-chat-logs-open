# 交換干渉散乱行列によるフェルミオン反射系列差分監査 v1

## 目的

20260710 および 20260711 のフェルミオン反射・弾性衝突読出し系列について、旧 JSON の同一条件を読み直し、交換干渉散乱行列版で旧主読出しが再現されるかを監査した。

ここでの主判定は、名前の毛ではない。旧有効ケースに対して、二つの出力チャネルが旧ターゲットの運動量位相 `q` と状態形状を再現するかで判定する。

## 総合判定

| 指標 | 値 |
|---|---:|
| 監査対象実験数 | `19` |
| 監査対象ケース数 | `316` |
| 旧有効ケース数 | `255` |
| 完全反射ターゲットの旧有効ケース数 | `189` |
| 非完全反射または一般化ターゲットの旧有効ケース数 | `66` |
| 散乱行列完全反射で旧主読出しを再現したケース数 | `189` |
| 確率比のみなら最適近似できるケース数 | `192` |
| V2再計算検討が必要な旧有効ケース数 | `66` |
| 全旧有効ケースを再現 | `false` |

## 実験別監査

| series | experiment | legacy valid | complete target | non-complete target | scattering reproduced | best probability p-fit | needs V2 | all reproduced |
|---|---|---:|---:|---:|---:|---:|---:|---|
| 20260710 | elastic_collision_simulation | 1 | 1 | 0 | 1 | 1 | 0 | `true` |
| 20260710 | fermionic_interference_reflection | 1 | 1 | 0 | 1 | 1 | 0 | `true` |
| 20260710 | elastic_collision_multi_collision | 1 | 1 | 0 | 1 | 1 | 0 | `true` |
| 20260710 | elastic_collision_control_maps | 2 | 2 | 0 | 2 | 2 | 0 | `true` |
| 20260710 | elastic_collision_cell_resolution_sweep | 57 | 57 | 0 | 57 | 57 | 0 | `true` |
| 20260710 | elastic_collision_label_robustness | 20 | 20 | 0 | 20 | 20 | 0 | `true` |
| 20260710 | elastic_collision_eta_resolution_sweep | 59 | 59 | 0 | 59 | 59 | 0 | `true` |
| 20260710 | elastic_collision_observer_sweep | 7 | 7 | 0 | 7 | 7 | 0 | `true` |
| 20260710 | elastic_collision_observation_perturbation | 5 | 5 | 0 | 5 | 5 | 0 | `true` |
| 20260710 | elastic_collision_asymmetry_sweep | 7 | 7 | 0 | 7 | 7 | 0 | `true` |
| 20260711 | abc_multigauge_interference_readout | 1 | 1 | 0 | 1 | 1 | 0 | `true` |
| 20260711 | abc_multigauge_interference_readout_multi_collision | 1 | 1 | 0 | 1 | 1 | 0 | `true` |
| 20260711 | abc_multigauge_interference_readout_asymmetric_amplitude_sweep | 8 | 8 | 0 | 8 | 8 | 0 | `true` |
| 20260711 | abc_multigauge_interference_readout_robustness_sweep | 5 | 5 | 0 | 5 | 5 | 0 | `true` |
| 20260711 | abc_multigauge_generalized_elastic_collision_readout | 8 | 1 | 7 | 1 | 1 | 7 | `false` |
| 20260711 | abc_multigauge_generalized_elastic_collision_velocity_sweep | 9 | 1 | 8 | 1 | 3 | 8 | `false` |
| 20260711 | abc_multigauge_generalized_elastic_collision_multi_collision | 4 | 1 | 3 | 1 | 1 | 3 | `false` |
| 20260711 | abc_multigauge_generalized_elastic_collision_noise_robustness | 47 | 11 | 36 | 11 | 11 | 36 | `false` |
| 20260711 | abc_multigauge_generalized_elastic_collision_extreme_R_sweep | 12 | 0 | 12 | 0 | 1 | 12 | `false` |

## 最大差分例

| experiment | case | target class | q target A | q target B | p minus | p plus | max p error | max copy distance | best R | best max p error |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| abc_multigauge_generalized_elastic_collision_velocity_sweep | c08_A1.50_B1.00_u1.80_v-0.20 | non_complete_or_generalized | 0.569231 | 2.56923 | -10.4671 | -8.73712 | 11.3064 | 0.645958 | 0 | 0.769231 |
| abc_multigauge_generalized_elastic_collision_multi_collision | c04_A1.50_B1.00_u1.80_v-0.20 | non_complete_or_generalized | 0.569231 | 2.56923 | -10.4671 | -8.73712 | 11.3064 | 0.645958 | 0 | 0.769231 |
| abc_multigauge_generalized_elastic_collision_noise_robustness | c04_A1.50_B1.00_u1.80_v-0.20 | non_complete_or_generalized | 0.569231 | 2.56923 | -10.4671 | -8.73712 | 11.3064 | 0.645958 | 0 | 0.769231 |
| abc_multigauge_generalized_elastic_collision_noise_robustness | c04_A1.50_B1.00_u1.80_v-0.20 | non_complete_or_generalized | 0.569231 | 2.56923 | -10.4671 | -8.73712 | 11.3064 | 0.645958 | 0 | 0.769231 |
| abc_multigauge_generalized_elastic_collision_noise_robustness | c04_A1.50_B1.00_u1.80_v-0.20 | non_complete_or_generalized | 0.569231 | 2.56923 | -10.4671 | -8.73712 | 11.3064 | 0.645958 | 0 | 0.769231 |
| abc_multigauge_generalized_elastic_collision_noise_robustness | c04_A1.50_B1.00_u1.80_v-0.20 | non_complete_or_generalized | 0.569231 | 2.56923 | -10.4671 | -8.73712 | 11.3064 | 0.645958 | 0 | 0.769231 |
| abc_multigauge_generalized_elastic_collision_noise_robustness | c04_A1.50_B1.00_u1.80_v-0.20 | non_complete_or_generalized | 0.569231 | 2.56923 | -10.4671 | -8.73712 | 11.3064 | 0.645958 | 0 | 0.769231 |
| abc_multigauge_generalized_elastic_collision_noise_robustness | c04_A1.50_B1.00_u1.80_v-0.20 | non_complete_or_generalized | 0.569231 | 2.56923 | -10.4671 | -8.73712 | 11.3064 | 0.645958 | 0 | 0.769231 |
| abc_multigauge_generalized_elastic_collision_noise_robustness | c04_A1.50_B1.00_u1.80_v-0.20 | non_complete_or_generalized | 0.569231 | 2.56923 | -10.4671 | -8.73712 | 11.3064 | 0.645958 | 0 | 0.769231 |
| abc_multigauge_generalized_elastic_collision_noise_robustness | c04_A1.50_B1.00_u1.80_v-0.20 | non_complete_or_generalized | 0.569231 | 2.56923 | -10.4671 | -8.73712 | 11.3064 | 0.645958 | 0 | 0.769231 |
| abc_multigauge_generalized_elastic_collision_noise_robustness | c04_A1.50_B1.00_u1.80_v-0.20 | non_complete_or_generalized | 0.569231 | 2.56923 | -10.4671 | -8.73712 | 11.3064 | 0.645958 | 0 | 0.769231 |
| abc_multigauge_generalized_elastic_collision_noise_robustness | c04_A1.50_B1.00_u1.80_v-0.20 | non_complete_or_generalized | 0.569231 | 2.56923 | -10.4671 | -8.73712 | 11.3064 | 0.645958 | 0 | 0.769231 |
| abc_multigauge_generalized_elastic_collision_noise_robustness | c04_A1.50_B1.00_u1.80_v-0.20 | non_complete_or_generalized | 0.569231 | 2.56923 | -10.4671 | -8.73712 | 11.3064 | 0.645958 | 0 | 0.769231 |
| abc_multigauge_generalized_elastic_collision_noise_robustness | c04_A1.50_B1.00_u1.80_v-0.20 | non_complete_or_generalized | 0.569231 | 2.56923 | -10.4671 | -8.73712 | 11.3064 | 0.645958 | 0 | 0.769231 |
| abc_multigauge_generalized_elastic_collision_velocity_sweep | c07_A2.00_B1.00_u0.80_v-1.50 | non_complete_or_generalized | -0.12 | 2.18 | -9.60216 | 0.717987 | 9.48216 | 1.09002 | 1 | 0.68 |
| abc_multigauge_generalized_elastic_collision_velocity_sweep | c03_A1.00_B1.00_u0.50_v-1.70 | non_complete_or_generalized | -1.7 | 0.5 | -0.5 | 9.97951 | 9.47951 | 1.2277 | 0 | 0 |

## 判定

完全反射ターゲット、すなわち旧ターゲットが `q_A -> -q_A`, `q_B -> -q_B` であるケースでは、交換干渉散乱行列の `R=1` 極限が旧主読出しを再現する。

一方、一般化弾性衝突や非完全反射ターゲットでは、単純な `R=1` 反射だけでは旧ターゲットを再現しない。さらに、`R:T` の確率比だけを調整しても再現できないケースは、散乱確率ではなく、出力キャリア位相そのものを更新する一般化散乱写像が必要である。

したがって、20260710/20260711 側を V2 化する必要性は一律ではない。完全反射の主結論だけに依存する論文は、散乱行列版でも主読出しが保持される。一般化速度、非等振幅、ABC 多ゲージ一般化弾性衝突に依存する論文は、V2 再計算の対象である。

## 出力

| 種別 | ファイル |
|---|---|
| JSON | `exchange_scattering_matrix_fermionic_reflection_legacy_audit_result_v1.json` |
| case CSV | `exchange_scattering_matrix_fermionic_reflection_legacy_audit_cases_v1.csv` |
| experiment CSV | `exchange_scattering_matrix_fermionic_reflection_legacy_audit_experiments_v1.csv` |
| report | `exchange_scattering_matrix_fermionic_reflection_legacy_audit_report_v1.md` |
