# ABCベースライン定常波C 準安定傾斜の幾何スケーリング予備実験 v1

## 目的

準安定傾斜から生じる加速度候補が、距離位相 `sin(Δχ)` と R重みに従うかを調べる。

ここでは定常場の大きさではなく、同じ ramp duration の `ΔC_memory` が、距離位相と R重みによりどう変わるかを読む。

## 統合判定

- case_count: `15`
- family_count: `2`
- single_gauge_only_used: `False`
- ramp_duration: `8`
- all_family_correlations_near_one: `True`
- all_family_relative_errors_small: `True`
- all_family_signs_valid: `True`
- max_balance_error_abs: `1.3234889800848443e-23`
- metastable_geometry_scaling_preliminary_valid: `True`

## family summary

| family | cases | corr(actual, expected) | rel error | norm std | signs | balance |
|---|---:|---:|---:|---:|---|---:|
| distance_phase | 8 | 1.0 | 3.2373070891891948e-15 | 2.7101212931035675e-24 | `True` | 1.3234889800848443e-23 |
| R_weight | 7 | 1.0 | 6.2112882302027380e-16 | 5.4712854753734605e-25 | `True` | 1.3234889800848443e-23 |

## case summary

| case | family | distance | R_B/R_A | max |R*a| | expected | normalized | sign ok |
|---|---|---:|---:|---:|---:|---:|---|
| distance_-1.00 | distance_phase | -1.000000 | 1.562500 | 4.8884533527749182e-10 | 3.1286101457759391e-01 | 1.5625000000000043e-09 | `True` |
| distance_-0.75 | distance_phase | -0.750000 | 1.562500 | 3.9599217821849469e-10 | 2.5343499405983572e-01 | 1.5625000000000055e-09 | `True` |
| distance_-0.50 | distance_phase | -0.500000 | 1.562500 | 2.7851814547481728e-10 | 1.7825161310388274e-01 | 1.5625000000000028e-09 | `True` |
| distance_-0.25 | distance_phase | -0.250000 | 1.562500 | 1.4372720342623236e-10 | 9.1985410192788122e-02 | 1.5625000000000101e-09 | `True` |
| distance_+0.25 | distance_phase | 0.250000 | 1.562500 | 1.4372720342623236e-10 | 9.1985410192788122e-02 | 1.5625000000000101e-09 | `True` |
| distance_+0.50 | distance_phase | 0.500000 | 1.562500 | 2.7851814547481728e-10 | 1.7825161310388274e-01 | 1.5625000000000028e-09 | `True` |
| distance_+0.75 | distance_phase | 0.750000 | 1.562500 | 3.9599217821849469e-10 | 2.5343499405983572e-01 | 1.5625000000000055e-09 | `True` |
| distance_+1.00 | distance_phase | 1.000000 | 1.562500 | 4.8884533527749182e-10 | 3.1286101457759391e-01 | 1.5625000000000043e-09 | `True` |
| R_ratio_0.25 | R_weight | 0.500000 | 0.250000 | 2.9964096162762731e-11 | 1.9177021544168123e-02 | 1.5625000000000020e-09 | `True` |
| R_ratio_0.5 | R_weight | 0.500000 | 0.500000 | 8.3233600452118664e-11 | 5.3269504289355889e-02 | 1.5625000000000016e-09 | `True` |
| R_ratio_1 | R_weight | 0.500000 | 1.000000 | 1.8727560101726717e-10 | 1.1985638465105075e-01 | 1.5625000000000030e-09 | `True` |
| R_ratio_1.562 | R_weight | 0.500000 | 1.562500 | 2.7851814547481728e-10 | 1.7825161310388274e-01 | 1.5625000000000028e-09 | `True` |
| R_ratio_2 | R_weight | 0.500000 | 2.000000 | 3.3293440180847465e-10 | 2.1307801715742355e-01 | 1.5625000000000016e-09 | `True` |
| R_ratio_4 | R_weight | 0.500000 | 4.000000 | 4.7942553860420370e-10 | 3.0683234470668996e-01 | 1.5625000000000020e-09 | `True` |
| R_ratio_8 | R_weight | 0.500000 | 8.000000 | 5.9188338099284420e-10 | 3.7880536383541963e-01 | 1.5625000000000026e-09 | `True` |

## 解釈

- 距離位相スイープで `|R*a|` が `|sin(Δχ)|` に比例するなら、準安定傾斜には幾何的な距離位相依存がある。
- R比スイープで `|R*a|` が `R_red^2 |sin(Δχ)|` に比例するなら、Cへの書き込みとCからの戻りの両方に R重みが入っている。
- `A` と `B` の符号が反対で、`A` の符号が `sin(Δχ)` と対応するなら、ラベル固定の押し込みではなく相対位相の向きに従っている。
- これは標準重力の距離法則ではない。準安定傾斜候補に距離位相・R重みのスケーリングがあるかを見る予備試験である。

## 出力

| 種類 | ファイル |
|---|---|
| JSON | `abc_baseline_stationary_wave_metastable_geometry_scaling_preliminary_result_v1.json` |
| case CSV | `abc_baseline_stationary_wave_metastable_geometry_scaling_cases_v1.csv` |
| rows CSV | `abc_baseline_stationary_wave_metastable_geometry_scaling_rows_v1.csv` |
