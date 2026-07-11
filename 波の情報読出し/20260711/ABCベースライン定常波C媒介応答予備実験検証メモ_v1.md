# ABCベースライン定常波C媒介応答予備実験 v1

## 目的

`A,B -> C` で C が歪むが `A,B` は動かない条件と、`C -> A,B` に戻した場合の応答を分離して測る。

本予備実験は、重力的読出しの成立ではなく、C媒介経路が読出し器擾乱、frozen応答、persistent戻入応答に分けられるかを調べる。

## 統合判定

- phase2_c_deformation_only_valid: `True`
- phase3_c_return_frozen_valid: `True`
- phase3_c_return_persistent_valid: `True`
- single_gauge_only_used: `False`
- c_mediated_response_preliminary_valid: `True`
- max_phase2_C_memory_abs: `3.5873847375164163e-06`
- max_frozen_C_memory_abs: `2.923326454903677e-07`
- max_persistent_C_memory_abs: `3.5873841187439837e-06`
- persistent_distance_change: `-1.3414465233063666e-07`
- persistent_R_delta_balance: `-7.26727798964588e-20`
- persistent_R_acceleration_balance_max: `4.011825970882184e-23`
- gauge_count: `7`

## フェーズ別サマリー

| phase | case | valid | max A delta | max B delta | distance change | R delta balance | R a balance |
|---|---|---|---:|---:|---:|---:|---:|
| phase2_c_deformation_only | A_B_to_C_no_return | `True` | 2.8365488260742951e-20 | 2.8365488260742951e-20 | 0.0000000000000000e+00 | -7.2686563668153807e-20 | 0.0000000000000000e+00 |
| phase3_c_return_frozen | C_return_without_reembedding | `True` | 1.7825161307551725e-10 | 1.1408103241485042e-10 | 0.0000000000000000e+00 | -7.2686531774142493e-20 | 0.0000000000000000e+00 |
| phase3_c_return_persistent | C_return_with_persistent_reembedding | `True` | 8.1795519756403468e-08 | 5.2349132644144730e-08 | -1.3414465233063666e-07 | -7.2672779896458799e-20 | 4.0118259708821842e-23 |

## 解釈

- Phase 2 では `A,B -> C` により C の記憶変数は非ゼロになるが、`C -> A,B` を切るため A/B の `δχ` はゼロでなければならない。
- frozen return では C からの瞬間応答を読むが、状態へ戻入しないため蓄積や加速度は出ない。
- persistent return では C 残渣を次ステップへ残し、readout reembedding を行うため、距離位相の縮小、R重み付き収支、微小加速度候補を同時に見る。
- ここでの応答は本命効果ではなく、Stage I 粗計量へ進むための C媒介経路の床確認である。

## 出力

| 種類 | ファイル |
|---|---|
| JSON | `abc_baseline_stationary_wave_c_mediated_response_preliminary_result_v1.json` |
| timeline CSV | `abc_baseline_stationary_wave_c_mediated_response_timeline_v1.csv` |
| gauge CSV | `abc_baseline_stationary_wave_c_mediated_response_gauge_rows_v1.csv` |
| C memory CSV | `abc_baseline_stationary_wave_c_mediated_response_c_memory_v1.csv` |
