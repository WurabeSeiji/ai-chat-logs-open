# 再構成スクリプト検証記録（2026-08-28）

## 判定の原則

- **RECOVERED_ORIGINAL**: 元 package / Drive folder に残っていた原本スクリプト。
- **RECONSTRUCTED_EQUIVALENT**: 原本スクリプトは残っていないが、保存された基礎データ・列定義・数値出力から生成則を逆算し、旧 CSV と比較した等価生成器。
- **NOT_RECOVERED**: 歴史的な実行値を裏付ける記録が無く、推測で埋めていないもの。

## B-1 N=5

`B1_N5/reconstruct_N5_missing_outputs.py` は engine を seed `40265722` で再実行し、5000 step の軌道から7 CSVを生成する。

旧出力との比較:

- `N5_all_steps_a_b_a2_b2_ab.csv`: 最大数値差 `8.88e-16`
- `N5_final_four_group_pattern_by_step.csv`: 最大数値差 `9.02e-16`
- `N5_key_step_ab2_ab_classes.csv`: 最大数値差 `6.11e-16`
- `N5_relation_class_counts_by_step.csv`: 一致
- `N5_step5000_four_group_summary.csv`: 最大数値差 `1.11e-16`
- `N5_time_separation_milestones.csv`: 一致
- `N5_inflation_vs_ordering_timeseries.csv`: 最大数値差 `9.02e-16`

4群誤差が以後ずっと閾値未満になる step も再現:

- `1e-4`: 2627
- `1e-6`: 3791
- `1e-8`: 4923

7枚の PNG は失われた描画原本を数値データから再構成した **科学内容等価図**。PNG のピクセル SHA 一致を主張しない。

## B-2 N=16

`B2_N16/reconstruct_N16_missing_outputs.py` は、既存の `N16_selected_snapshots_long.csv`, `N16_global_summary.csv`, `decompact_N16_geometry_summary.csv` から4ファイルを生成する。

- `N16_step5000_final_edges_with_normalized_components.csv`: 文字列レベル一致
- `N16_triplet_cluster_counts_selected_steps.csv`: 文字列レベル一致
- `N16_time_separation_milestones.csv`: 文字列レベル一致
- `N16_phase_structure_summary.json`: 全キー・全値一致

クラスタリングは「最初に現れたメンバーを固定中心とする greedy max-component distance」であることを旧出力全行から確認した。

## B-3 N=3/N=4

`B3_N3_N4/reconstruct_N3_N4_phase_outputs.py` の6 CSV（N3 3本 + N4 3本）および `N3_N4_comparison_summary.csv` は旧 CSV と `1e-14` 以下で一致。文字列列も一致。

6枚の PNG は科学内容等価図として再構成。

## B-4 comparison summary 6件

`B4_comparison/build_pair_comparison_summary.py` で以下6ファイルを生成し、**全6件を旧 CSV と文字列レベル一致**で確認。

- N3_N4
- N6_N7
- N8_N9
- N10_N11
- N12_N13
- N14_N15

## B-5 非自明ゼロ閉包

`B5_nontrivial/rebuild_nontrivial_zero_closure_analysis.py`:

- small-subset classification 69行: 一致
- summary 14行: 一致
- assessment 14行: 一致
- N5 exact covers 12通り: 一致
- N5 best-pair time evolution 5001行: 一致

`B5_nontrivial/mitm_subset_search.py` で N=14,k=6 を再探索し、旧結果と同じ
`1-7,2-14,4-10,5-9,5-14,10-13`, residual `9.249993443238053e-07`
を再取得した。runtime 列はマシン依存なので一致対象外。

## B-6 N=14..16 全 cardinality search

回収済み原本:

- `anneal_subsets.cpp`
- `exact_k234.cpp`
- `search_subsets.cpp`
- orphaned historical `mitm56.py`, `refine_subset.py`

再構成:

- `aggregate_search_results.py`: preserved method-level resultsを入力した検証で `N14_N16_all_cardinalities_best_results.csv` と `N14_N16_summary.csv` が旧 CSV と一致。
- `build_time_evolution_and_n5_covers.py`: N14 の engine 軌道を再実行して検証し、N14/N5 time evolution と N5 exact covers が旧 CSV と一致。N14 time evolution の最大数値差 `9.97e-17`。

### 唯一の未回収点

歴史的な `anneal_subsets N k input steps seed` の **steps/seed 全21ケース**は package / 利用可能な保存記録から回収できなかった。

そのため:

- `historical_anneal_args_NOT_RECOVERED.md` に欠損を明記。
- `reproduction_anneal_args.csv` に今後用の新しい固定値を与えた。
- この新しい値を歴史的引数とは主張しない。
- `run_all_reconstructed.py` はデフォルトで B6 stochastic search を実行せず、`--run-new-stochastic-search` 指定時だけ新しい固定値で実行する。

## B-7 partial zero closure

監査の結果、

`N3_N16_zero_subset_exact_covers.csv`
と
`N3_N16_zero_triple_exact_covers.csv`

は SHA256 も内容も完全同一だった。元解析スクリプトは後者だけを書いていたため、前者は計算の別枝ではなく filename alias と判定した。

`B7_partial/write_zero_subset_alias.py` および patched original を同梱。byte-identical を確認済み。

## B-8 K/sigma

ユーザー判断どおり廃止枝として再構成対象外。
