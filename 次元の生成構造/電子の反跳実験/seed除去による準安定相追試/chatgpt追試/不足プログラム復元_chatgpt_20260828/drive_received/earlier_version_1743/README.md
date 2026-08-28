# 論文 v1 欠落生成プログラム回収・再構成パッケージ (2026-08-28)

目的: 公開再現パッケージで出力だけ残り、生成 `.py` / 実行条件が欠落していた B-1〜B-7 を、(A) 回収できた原本と (B) 保存出力から数式を逆算して検証した等価生成器に分けて収録する。K/σ 正規化枝 (B-8) は依頼どおり対象外。

## 重要な区分

- `recovered_original/`: Drive 上の旧実験フォルダから回収した**原本スクリプト**。改変していない。
- `reconstructed_validated/`: 欠落していた一回限りの派生生成処理を、保存 CSV/JSON と照合して再構成したスクリプト。
- `validation/VALIDATION_REPORT.md`: 旧保存出力との機械比較結果。
- `HISTORICAL_ANNEAL_ARGS_UNRECOVERABLE.md`: B-6 の historical `steps,seed` を回収できなかった事実と、再現用の代替条件を区別して記録。

## B-1 N=5 complex simplex

`B1_generate_N5_complete_analysis.py`

入力: `N5_phase_by_edge_5000steps.csv`, `N5_raw_K_raw_observables.csv`。
要求された 7 CSV と 7 PNG を生成する。4群 A+/A-/B+/B- は保存最終状態の辺群を固定し、各群中心の `(a^2/r^2,b^2/r^2,ab/r^2)` 最大偏差と 90° complement error から `four_group_error` を計算する。旧保存値の「以後ずっと 1e-4 未満 = step 2627」「以後ずっと 1e-8 未満 = step 4923」を数値再現した。

## B-2 N=16 phase outputs

`B2_generate_N16_phase_outputs.py`

入力: `N16_selected_snapshots_long.csv`, `N16_global_summary.csv`, `decompact_N16_geometry_summary.csv`。要求された 4 ファイルを生成し、旧保存値と一致する。

## B-3 N=3/N=4 phase analysis

`B3_generate_N3_N4_phase_outputs.py`

入力: `N3_all_steps_long.csv`, `N4_all_steps_long.csv`。要求された 6 CSV と 6 PNG を生成する。N=3 は pairwise circular phase distance の 1/3 収束、N=4 は対辺 (12,34),(13,24),(14,23) の三クラスを評価する。

## B-4 pair comparison summaries

`B4_generate_pair_comparison_summaries.py`

N3/N4, N6/N7, N8/N9, N10/N11, N12/N13, N14/N15 の 6 比較 CSV を、各パッケージの `N*_summary.json` / `N*_triplet_cluster_counts.csv` から生成する。6件すべて旧保存 CSV と一致確認済み。

## B-5 nontrivial zero closure

`B5_generate_nontrivial_zero_closure_outputs.py`

入力: `SOURCE_N3..16_step5000_final_edges.csv`。小部分集合は保存版と同じ exhaustive 判定 (tol=1e-6, star-span test)、大きい k は MITM。N=5 exact cover は 12 通りを旧保存順まで再現する。時系列 2 本と比較図を作る場合は `--n5-history` と `--n14-history` を指定する。

## B-6 N14-N16 complete search

原本: `exact_k234.cpp`, `anneal_subsets.cpp`, `search_subsets.cpp`, `mitm56.py`, `refine_subset.py`。

再構成:
- `B6_mitm_generic.py`
- `B6_aggregate_search_results.py`
- `B6_run_searches_reconstructed.sh`
- `B6_generate_shared_time_evolution.py`
- `anneal_args_reconstructed.csv`

k=2..4 は recovered original の完全総当たり、k=5,6 は MITM、k>=7 は recovered original annealer を使う。**historical annealing の `steps,seed` は保存されておらず回収不能**なので、再実行用テーブルは明示的に `reconstructed_not_historical` とした。旧結果を再現したと偽装しない。

共有時系列は `B6_generate_shared_time_evolution.py` が `N5_nontrivial_pair_exact_covers.csv`, `N5_best_nontrivial_pair_time_evolution.csv`, `N14_best_k6_candidate_time_evolution.csv`, `N14_candidate_vs_N5_benchmark_time_evolution.png` を生成する。

## B-7 partial zero closure exact covers

`B7_generate_zero_subset_exact_covers.py`

旧 `analyze_partial_zero_closures_N3_N16.py` が `N3_N16_zero_triple_exact_covers.csv` として出していた同一表を、要求名 `N3_N16_zero_subset_exact_covers.csv` で出力する。17行8列が旧表と完全一致確認済み。

## 図について

CSV/JSON の数値生成則を優先して復元した。再構成 PNG はその再生成 CSV から作る等価図であり、旧一回限り作図コードが残っていない図について**ピクセル同一**とは主張しない。

## 最小実行例

```bash
python3 reconstructed_validated/B1_generate_N5_complete_analysis.py --input-dir N5_complex_simplex_complete_analysis_20260826 --output-dir out/B1
python3 reconstructed_validated/B2_generate_N16_phase_outputs.py --input-dir N16_complex_simplex_complete_analysis_20260826 --output-dir out/B2
python3 reconstructed_validated/B3_generate_N3_N4_phase_outputs.py --input-dir N3_N4_complex_simplex_complete_analysis_20260826 --output-dir out/B3
python3 reconstructed_validated/B5_generate_nontrivial_zero_closure_outputs.py --source-dir N3_N16_nontrivial_zero_closure_analysis_20260826 --output-dir out/B5
python3 reconstructed_validated/B7_generate_zero_subset_exact_covers.py --source-dir N3_N16_partial_zero_closure_analysis_20260826 --output out/B7/N3_N16_zero_subset_exact_covers.csv
```
