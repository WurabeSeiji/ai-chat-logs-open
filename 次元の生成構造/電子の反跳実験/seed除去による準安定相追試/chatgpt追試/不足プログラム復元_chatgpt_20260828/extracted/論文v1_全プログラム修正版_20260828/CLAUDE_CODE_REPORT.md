# Claude Code 引継ぎ報告書 — 論文v1 全再現テスト不足プログラム復元

作成日: 2026-08-28  
対象: `論文v1_全再現テスト_20260828` の NOT_REGENERATED 項目のうち B-1〜B-7  
方針: **既存出力の再送ではなく、出力を再生成するプログラムを回収・再構成する。**  
B-8 K/σ 枝: ユーザー指示により対象外。

---

## 1. Google Drive 上の基準位置

Google Drive API の親フォルダ列をたどって確認した論理パスは次のとおり。

```text
マイドライブ/
└── OneDrive/
    └── GitHub/
        └── ai-chat-logs-open/
            └── 次元の生成構造/
                └── 電子の反跳実験/
                    └── seed除去による準安定相追試/
                        └── chatgpt追試/
                            └── 論文v1_全再現テスト_20260828/
```

**Google Drive 論理パス**

```text
マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/論文v1_全再現テスト_20260828
```

**Drive folder ID**

```text
14k8LWkk7Ht0EA4Yv-r4wByuySWFymfQF
```

**Drive URL**

```text
https://drive.google.com/drive/folders/14k8LWkk7Ht0EA4Yv-r4wByuySWFymfQF
```

Google Drive for desktop のローカル実体パスは環境依存なので固定文字列を推測しないこと。Claude Code 側では、上記の `マイドライブ/OneDrive/...` を実際の同期ルートに接続する。

例:

```bash
GD_ROOT="<Google Drive for desktop のマイドライブ実体パス>"
SOURCE_ROOT="$GD_ROOT/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/論文v1_全再現テスト_20260828"
```

---

## 2. 今回作成した復元パッケージ

パッケージルート:

```text
論文v1_全プログラム修正版_20260828/
```

主要ファイル:

```text
README.md
VALIDATION.md
CLAUDE_CODE_REPORT.md
coverage_manifest.csv
requirements.txt
CHECKSUMS.sha256
copy_A_aliases.py
run_all_reconstructed.py
recovered_original/
B1_N5/
B2_N16/
B3_N3_N4/
B4_comparison/
B5_nontrivial/
B6_all_cardinalities/
B7_partial/
```

ZIP:

```text
論文v1_全プログラム修正版_20260828.zip
```

この ZIP は旧結果ファイルを集め直したものではない。**不足していた生成器、回収できた原本、検証用ドキュメントをまとめたもの**である。

---

## 3. Claude Code が最初に読むファイル

順番は次のとおり。

```text
1. CLAUDE_CODE_REPORT.md
2. README.md
3. VALIDATION.md
4. coverage_manifest.csv
5. B6_all_cardinalities/historical_anneal_args_NOT_RECOVERED.md
6. B6_all_cardinalities/anneal_args_audit.csv
```

特に B-6 について、欠損した歴史的引数を推測で埋めないこと。

---

## 4. 回収・再構成した基礎 generator

`recovered_original/base_generators/` に以下を収録した。

```text
run_n_scaling_lowrank_v1_no_sigma_norm.py
run_N3_N4_complete_analysis.py
run_N6_N7_complete_analysis.py
run_N8_N9_complete_analysis.py
run_N10_N11_complete_analysis.py
run_N12_N13_complete_analysis.py
run_N14_N15_complete_analysis.py
run_N16_complex_simplex_physics.py
analyze_and_plot_N16.py
```

N5 provenance:

```text
recovered_original/N5_provenance/analyze_N5_complex_simplex_structure.py
recovered_original/N5_provenance/plot_N5_inflation_vs_ordering.py
recovered_original/N5_provenance/run_N5_physical_phase_step_test.py
```

これらを基礎に B-1〜B-7 の不足出力を生成する。

---

## 5. B-1 N=5

生成器:

```text
B1_N5/reconstruct_N5_missing_outputs.py
```

engine を N=5, seed=`40265722`, 5000 step で再実行し、保存されていた N5 軌道と照合した。

数値一致:

```text
N5_all_steps_a_b_a2_b2_ab.csv             max error 8.88e-16
N5_final_four_group_pattern_by_step.csv    max error 9.02e-16
N5_key_step_ab2_ab_classes.csv             max error 6.11e-16
N5_relation_class_counts_by_step.csv       exact
N5_step5000_four_group_summary.csv         max error 1.11e-16
N5_time_separation_milestones.csv          exact
N5_inflation_vs_ordering_timeseries.csv    max error 9.02e-16
```

4群誤差が以後ずっと閾値未満になる step:

```text
1e-4 -> 2627
1e-6 -> 3791
1e-8 -> 4923
```

7図も再生成する。ただし失われた旧描画スクリプトそのものではないため、PNG pixel SHA 一致ではなく**科学内容等価**である。

---

## 6. B-2 N=16

生成器:

```text
B2_N16/reconstruct_N16_missing_outputs.py
```

入力:

```text
N16_selected_snapshots_long.csv
N16_global_summary.csv
decompact_N16_geometry_summary.csv
```

再生成対象:

```text
N16_phase_structure_summary.json
N16_step5000_final_edges_with_normalized_components.csv
N16_time_separation_milestones.csv
N16_triplet_cluster_counts_selected_steps.csv
```

旧保存結果と比較し、CSV 3件は文字列レベル一致、JSON は全 key/value 一致を確認済み。

---

## 7. B-3 N=3/N=4 phase outputs

生成器:

```text
B3_N3_N4/reconstruct_N3_N4_phase_outputs.py
```

N3/N4 の 6 CSV と 6図、および `N3_N4_comparison_summary.csv` を生成する。

数値 CSV は旧結果と `1e-14` 以下で一致。図は科学内容等価。

---

## 8. B-4 comparison summaries

生成器:

```text
B4_comparison/build_pair_comparison_summary.py
B4_comparison/run_all_six_comparisons.py
```

再生成:

```text
N3_N4_comparison_summary.csv
N6_N7_comparison_summary.csv
N8_N9_comparison_summary.csv
N10_N11_comparison_summary.csv
N12_N13_comparison_summary.csv
N14_N15_comparison_summary.csv
```

6件とも旧 CSV と文字列レベル一致を確認済み。

---

## 9. B-5 非自明ゼロ閉包

主要生成器:

```text
B5_nontrivial/mitm_subset_search.py
B5_nontrivial/run_mitm_surveys.py
B5_nontrivial/rebuild_nontrivial_zero_closure_analysis.py
B5_nontrivial/exact_k234.cpp
```

検証済み:

```text
small-subset classification 69行     exact
summary 14行                         exact
assessment 14行                      exact
N5 exact covers 12通り               exact
N5 best-pair time evolution 5001行   exact
```

N=14,k=6 MITM も旧結果と同じ subset を再取得:

```text
1-7,2-14,4-10,5-9,5-14,10-13
residual = 9.249993443238053e-07
```

runtime はマシン依存なので一致判定対象外。

---

## 10. B-6 N=14..16 全 cardinality search

### 回収できた原本

```text
B6_all_cardinalities/anneal_subsets.cpp
B6_all_cardinalities/exact_k234.cpp
B6_all_cardinalities/search_subsets.cpp
B6_all_cardinalities/mitm56_historical_orphan.py
B6_all_cardinalities/refine_subset_historical_orphan.py
```

historical orphan 2本は古い絶対パスを含むため、そのまま再現入口には使わず provenance として保持した。

### 再現用に整備したスクリプト

```text
B6_all_cardinalities/mitm_repro.py
B6_all_cardinalities/refine_subset_repro.py
B6_all_cardinalities/run_search_methods.py
B6_all_cardinalities/aggregate_search_results.py
B6_all_cardinalities/build_time_evolution_and_n5_covers.py
```

### 重要: 唯一の歴史的未回収事項

旧 `anneal_subsets N k input steps seed` の **steps と seed が 21ケース分回収できていない**。

対象:

```text
N14: k=7..12   6 cases
N15: k=7..13   7 cases
N16: k=7..14   8 cases
TOTAL          21 cases
```

したがって、Claude Code は旧実験値を推測してはならない。

監査:

```text
B6_all_cardinalities/historical_anneal_args_NOT_RECOVERED.md
B6_all_cardinalities/anneal_args_audit.csv
```

新しい再現実験用の固定値は別ファイル:

```text
B6_all_cardinalities/reproduction_anneal_args.csv
```

この値は **NEW_REPRODUCTION_NOT_HISTORICAL** であり、旧実験の再発見値ではない。

`run_all_reconstructed.py` はデフォルトでは B6 stochastic search を実行しない。`--run-new-stochastic-search` を明示したときだけ新しい固定値を使う。

---

## 11. B-7 partial zero closure

監査により次の2ファイルは byte-identical と判定した。

```text
N3_N16_zero_subset_exact_covers.csv
N3_N16_zero_triple_exact_covers.csv
```

生成器:

```text
B7_partial/write_zero_subset_alias.py
```

また、元解析スクリプトの patched 版も同梱:

```text
B7_partial/analyze_partial_zero_closures_N3_N16_patched.py
```

これは別物理計算の枝ではなく filename alias である。

---

## 12. A群コピー処理

コピー専用:

```text
copy_A_aliases.py
```

`SOURCE_N{3..16}_step5000_final_edges.csv` 等を必要な解析フォルダへコピーする。

N5 については完全 package に `N5_step5000_final_edges.csv` が無かったため、`N5_all_steps_a_b_a2_b2_ab.csv` の step=5000 と K5 edge order から必要な source を構成する。

---

## 13. 全体実行入口

通常の復元実行:

```bash
cd "<論文v1_全プログラム修正版_20260828 の展開先>"

python run_all_reconstructed.py \
  --root "$SOURCE_ROOT" \
  --decompact-results "<decompactification results path>" \
  --n5-raw-k-source "<N5 raw K observables CSV>"
```

このモードでは B6 k>=7 annealing は走らない。

新しい固定 seed/steps による B6 再現実験を行う場合のみ:

```bash
python run_all_reconstructed.py \
  --root "$SOURCE_ROOT" \
  --decompact-results "<decompactification results path>" \
  --n5-raw-k-source "<N5 raw K observables CSV>" \
  --run-new-stochastic-search
```

**旧 steps/seed を再現したことにはならない。**

---

## 14. 依存環境

Python:

```text
numpy
pandas
matplotlib
scipy
```

B5/B6 exact / annealing search:

```text
C++17 compiler (g++)
```

全 Python script は syntax compile 済み。B6 C++ source も compile 確認済み。

---

## 15. Claude Code への禁止事項

1. B6 の historical `steps/seed` を推測で補わない。
2. `reproduction_anneal_args.csv` を historical data と呼ばない。
3. B1/B3 の再構成図について pixel/SHA identical と主張しない。
4. B8 K/σ 枝を勝手に復活させない。
5. 旧 CSV が存在するからといって、それをコピーして「再生成」と扱わない。
6. 原本 generator が残っている箇所は `recovered_original/` を優先して読む。
7. Drive のローカル同期絶対パスを推測せず、上記 Drive 論理パスを実環境で解決する。

---

## 16. 検証資料

詳細は:

```text
VALIDATION.md
coverage_manifest.csv
CHECKSUMS.sha256
```

を参照すること。

今回の復元範囲における結論は、**B-1〜B-5 および B-7 は生成ロジックを再構成・検証済み、B-6 は生成系自体を回収・再構成できたが historical anneal steps/seed だけが未回収**、である。
