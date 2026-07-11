# ABCベースライン定常波Cにおける反力候補粗計量予備実験 v1

## 目的

Stage I の前段として、引力様に見える距離位相縮小と混同し得る反力候補を、`R1`, `R2`, `R3` に分けて粗く計量する。

本実験は、標準力や標準重力を導入しない。読出し器依存の共通モード、C定常波圧、A-B直接干渉が、それぞれどの符号と桁で現れるかを調べる。

## 統合判定

- candidate_count: `3`
- case_count: `6`
- R1_valid: `True`
- R2_valid: `True`
- R3_valid: `True`
- all_candidates_valid: `True`
- single_gauge_only_used: `False`
- reaction_candidates_preliminary_valid: `True`
- max_R_weighted_delta_balance_for_R2_R3: `1.0164395367051604e-19`
- max_R_weighted_acceleration_balance_for_R2_R3: `3.8963515573697816e-20`
- max_gauge_delta_std: `6.429344348871916e-19`

## ケース別判定

| candidate | case | class | A delta | B delta | distance change | R balance | R a balance | valid |
|---|---|---|---:|---:|---:|---:|---:|---|
| R1 | C_reference_plus | reference_common_mode | 1.9999999999971638e-08 | 1.9999999999971638e-08 | 0.0000000000000000e+00 | 5.1249999999927320e-08 | 0.0000000000000000e+00 | `True` |
| R1 | C_reference_minus | reference_common_mode | -2.0000000000028366e-08 | -2.0000000000028366e-08 | 0.0000000000000000e+00 | -5.1250000000072685e-08 | 0.0000000000000000e+00 | `True` |
| R2 | C_pressure_repulsion | repulsive_distance_growth | -9.1350254629537857e-05 | 5.8464162962904213e-05 | 1.4981441759243097e-04 | -2.7105054312137611e-20 | 3.6422416731934915e-20 | `True` |
| R2 | C_pressure_repulsion_mirrored | repulsive_distance_growth | 9.1350254629537817e-05 | -5.8464162962904267e-05 | 1.4981441759243097e-04 | -9.4867690092481638e-20 | 3.8963515573697816e-20 | `True` |
| R3 | AB_direct_repulsion | repulsive_distance_growth | -5.7091905735827904e-05 | 3.6538819670929792e-05 | 9.3630725406690374e-05 | -1.0164395367051604e-19 | 2.5834504891256160e-20 | `True` |
| R3 | AB_direct_repulsion_mirrored | repulsive_distance_growth | 5.7091905735827850e-05 | -3.6538819670929846e-05 | 9.3630725406690374e-05 | -3.3881317890172014e-20 | 2.4140438996747560e-20 | `True` |

## 解釈

- `R1` は C 参照に依存した共通モードとして分類する。A/B が同符号に読まれ、距離位相を変えず、R重み付き外部収支を満たさないなら、外部並進候補ではなく読出し器反力または読出し器バイアスとして保留する。
- `R2` は C 定常波圧による離反候補として分類する。距離位相が増大し、R重み付き収支が小さいなら、引力様候補とは別の反力指紋として扱う。
- `R3` は A-B 直接干渉による離反候補として分類する。C媒介を切っても距離位相が増大するなら、C残渣由来の接近候補とは別系列として扱う。
- この予備実験で反力候補の桁が C媒介 persistent 応答と同程度なら、次は Stage II の同時多読出しへ進む必要がある。

## 出力

| 種類 | ファイル |
|---|---|
| JSON | `abc_baseline_stationary_wave_reaction_candidates_preliminary_result_v1.json` |
| case CSV | `abc_baseline_stationary_wave_reaction_candidates_cases_v1.csv` |
| gauge CSV | `abc_baseline_stationary_wave_reaction_candidates_gauge_rows_v1.csv` |
