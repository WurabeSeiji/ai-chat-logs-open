# ABCベースライン定常波C媒介応答の符号対称性予備実験 v1

## 目的

C媒介 persistent 応答が、A/B というラベルに固定された押し引きではなく、位置位相差の符号に従うかを確認する。

これは、引力的に見える距離位相縮小が、単なる実装上のラベル依存反力ではないかを切り分けるための予備テストである。

## 統合判定

- case_count: `4`
- all_cases_valid: `True`
- single_gauge_only_used: `False`
- c_mediated_symmetry_preliminary_valid: `True`
- max_R_weighted_delta_balance: `7.277865901486559e-20`
- max_R_weighted_acceleration_balance: `4.011825970882184e-23`
- max_gauge_delta_std: `6.414269100766159e-19`

## ケース別判定

| case | d_AB sign | A delta | B delta | distance change | expected signs A/B/d | observed signs A/B/d | sign valid | valid |
|---|---:|---:|---:|---:|---|---|---|---|
| normal | 5.0000000000000000e-01 | 8.1795519756403468e-08 | -5.2349132644144730e-08 | -1.3414465233063666e-07 | 1/-1/-1 | 1/-1/-1 | `True` | `True` |
| mirrored_positions | -5.0000000000000000e-01 | -8.1795519756460206e-08 | 5.2349132644087992e-08 | -1.3414465233063666e-07 | -1/1/-1 | -1/1/-1 | `True` | `True` |
| inverted_C_source | 5.0000000000000000e-01 | -8.1795532454306435e-08 | 5.2349140770709541e-08 | 1.3414467314731837e-07 | -1/1/1 | -1/1/1 | `True` | `True` |
| inverted_C_return | 5.0000000000000000e-01 | -8.1795532454306435e-08 | 5.2349140770709541e-08 | 1.3414467314731837e-07 | -1/1/1 | -1/1/1 | `True` | `True` |

## 解釈

- normal と mirrored_positions の両方で距離位相が縮む場合、応答は少なくともラベル固定の片側押しではない。
- inverted_C_source と inverted_C_return では距離位相が広がることを要求する。これにより、C 残渣の符号枝を反転したときに応答符号も反転するかを見る。
- R重み付き変位収支と加速度収支が小さいことを同時に要求する。
- この予備テストが失敗する場合、C媒介応答は Stage I 候補ではなく、実装ラベル依存または読出し枝依存として保留する。

## 出力

| 種類 | ファイル |
|---|---|
| JSON | `abc_baseline_stationary_wave_c_mediated_symmetry_preliminary_result_v1.json` |
| case CSV | `abc_baseline_stationary_wave_c_mediated_symmetry_cases_v1.csv` |
