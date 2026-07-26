# 第7論文 N=5 イベント定義校正報告（Stage A1）

## 1. 範囲と総合状態

- 実施内容: Stage A0で完全再現されたN=5 CSV 3件の後処理だけ
- 原本力学コードの実行: なし
- 新しい軌道生成: なし
- 第8論文本実験、高精度親、Series 1〜3、N=40、N=300: 未実行
- 単一イベント時刻の採用: なし
- 仮説H1/H2/H0の判定: なし
- **総合状態: `CALIBRATION_DATA_COMPLETE`**

本報告は候補集合の校正資料であり、指数成長開始、指数成長終了、rank_Q=4持続開始のいずれについても最終時刻を決定しない。

## 2. 入力SHA-256

| 入力 | 行数 | 絶対パス | SHA-256 | Stage A0照合 |
|---|---:|---|---|---|
| `fcurve` | 21168 | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_reproduction/reproduced/metastable_series_result_v1/fcurve_N00005_delta1e-15_seed0.csv` | `9220c5f3c1f570c8a52ea24a3cdd95568354cea0943d9bee7d8ed20316d3a9d0` | 一致 |
| `q_svd` | 995 | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_reproduction/reproduced/exact_lowN_eigenspectrum_v2/raw/N00005_dimension_saturation_v2/q_svd_N00005.csv` | `7c16a364c6cc9145293c2625dfe4ebb1f9962655d212679188215e8fad5e5155` | 一致 |
| `paper7_long` | 2201 | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_reproduction/reproduced/exact_lowN_eigenspectrum_v2/paper7_longtime/raw/N00005/paper7_long_timeseries.csv` | `efeaf9dab753c057ad0c6109b9e4a8919f8d8db1249da186658bfed9fda784e3` | 一致 |

照合はStage A0報告書の `REPRODUCED_EXACTLY` とCSV別 `exact` 行、Stage A0 CSV比較記録の `reproduced_sha256`、実ファイルSHA-256の四者で行った。

## 3. 使用列と時間軸

- `fcurve`: `tau`, `f`
- `q_svd`: `step`, `time`, `relative_time`, `q1`, `q2`, `q3`, `q4`, `rank_q`
- `paper7_long_timeseries`: `step`, `time`, `crossing_flag`, `splitting_fraction`

- 共通軸: `absolute_step`
- f範囲: `[0, 21167]`、1 step連続: `True`
- q範囲: `[0, 51167]`、`step-relative_time`: `[1167]`
- 長時間CSV範囲: `[0, 55000]`
- 既存crossing: `1167`
- qの未保存stepは補間していない。rank持続長は連続する保存レコード数であり、実step spanも候補表に併記した。

## 4. 全候補パラメータ

- 回帰窓: `[11, 21, 41, 81, 161, 321]`
- R²閾値: `[0.9, 0.95, 0.98, 0.99, 0.995, 0.999]`
- 成長区間最小持続長: `[10, 20, 40, 80, 160]`
- 成長終了条件: `{'A': {'description': 'slope <= 0', 'median_slope_multiplier': 0.0}, 'B': {'description': 'slope <= 0.1 * interval_median_slope', 'median_slope_multiplier': 0.1}, 'C': {'description': 'slope <= 0.01 * interval_median_slope', 'median_slope_multiplier': 0.01}}`
- 成長終了持続長: `[10, 20, 40, 80]`
- rank相対閾値: `[1e-06, 1e-07, 1e-08, 1e-09, 1e-10, 1e-11, 1e-12]`
- rank=4持続長（保存レコード数）: `[1, 5, 10, 20, 40, 80, 160]`

回帰対象は自然対数 `log_f` である。中心窓の通常最小二乗を用い、残差標準偏差は `sqrt(SSE/(window-2))` とした。`f<=0` は正数へ置換せず、対数不能点はNaNとした。

## 5. 各窓の成長区間候補数

| window | 候補行数 |
|---:|---:|
| 11 | 1765 |
| 21 | 1545 |
| 41 | 1175 |
| 81 | 642 |
| 161 | 177 |
| 321 | 116 |
| **全窓** | **5420** |

同じ最大連続区間が複数のR²閾値・最小持続長を満たす場合も、各パラメータ候補行を削除していない。

## 6. 候補分布

| 候補集合 | 件数 | 最小step | Q25 | 中央値 | Q75 | 最大step | 異なるstep数 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 成長区間開始 | 5420 | 7 | 5420 | 10705 | 17011.75 | 21039 | 1840 |
| 成長区間終了端 | 5420 | 1038 | 5444.25 | 10764 | 17066 | 21083 | 1872 |
| 成長終了候補（found） | 59623 | 1339 | 5337 | 10774 | 16305 | 21090 | 679 |
| rank=4持続開始候補（found） | 49 | 15 | 650 | 840 | 940 | 985 | 8 |
| rank開始 - 成長終了（全組合せ） | 2921527 | -21075 | -15620 | -9892 | -4497 | -354 | 4782 |

- 成長終了候補表全行: `65040`
- 成長終了 `found`: `59623`
- 成長終了 `not_found`: `5417`
- rank候補表全行: `49`
- rank候補 `found`: `49`
- 成長終了×rank開始の全ペア: `2921527`

## 7. パラメータに対して反復する候補群

以下は、異なる候補パラメータで同じstepが何回現れたかの頻度である。頻度順は採用順位ではない。

### 成長区間開始

| 頻度順位 | step | 出現候補行数 |
|---:|---:|---:|
| 1 | 1560 | 12 |
| 2 | 14 | 10 |
| 3 | 25 | 10 |
| 4 | 160 | 10 |
| 5 | 1554 | 10 |
| 6 | 4353 | 10 |
| 7 | 4075 | 9 |
| 8 | 6350 | 9 |
| 9 | 9922 | 9 |
| 10 | 2051 | 8 |

### 成長終了

| 頻度順位 | step | 出現候補行数 |
|---:|---:|---:|
| 1 | 8161 | 1253 |
| 2 | 4591 | 1085 |
| 3 | 1769 | 812 |
| 4 | 4592 | 760 |
| 5 | 3201 | 747 |
| 6 | 10916 | 625 |
| 7 | 2249 | 560 |
| 8 | 10777 | 509 |
| 9 | 6556 | 483 |
| 10 | 8160 | 480 |

### rank=4持続開始

| 頻度順位 | step | 出現候補行数 |
|---:|---:|---:|
| 1 | 840 | 16 |
| 2 | 940 | 7 |
| 3 | 985 | 7 |
| 4 | 885 | 6 |
| 5 | 15 | 4 |
| 6 | 475 | 4 |
| 7 | 650 | 4 |
| 8 | 265 | 1 |

これらは「パラメータに対して安定して反復する候補群」を抽出するための記述統計に限られ、どのstepも採用していない。

## 8. パラメータに強く依存する候補群

| パラメータ群 | 件数 | 最小 | 中央値 | 最大 | 異なるstep数 |
|---|---:|---:|---:|---:|---:|
| 成長開始 window=11 | 1765 | 7 | 11938 | 21030 | 606 |
| 成長開始 window=21 | 1545 | 12 | 11337 | 21037 | 603 |
| 成長開始 window=41 | 1175 | 22 | 11182 | 21036 | 505 |
| 成長開始 window=81 | 642 | 41 | 9922 | 21039 | 301 |
| 成長開始 window=161 | 177 | 80 | 4221 | 10563 | 56 |
| 成長開始 window=321 | 116 | 160 | 5022 | 10590 | 29 |
| 成長終了 condition=A | 19400 | 1416 | 10777 | 21090 | 243 |
| 成長終了 condition=B | 20823 | 1339 | 10774 | 21087 | 426 |
| 成長終了 condition=C | 19400 | 1407 | 10651 | 21089 | 283 |
| rank開始 threshold=1e-06 | 7 | 985 | 985 | 985 | 1 |
| rank開始 threshold=1e-07 | 7 | 940 | 940 | 940 | 1 |
| rank開始 threshold=1e-08 | 7 | 265 | 885 | 885 | 2 |
| rank開始 threshold=1e-09 | 7 | 15 | 840 | 840 | 4 |
| rank開始 threshold=1e-10 | 7 | 15 | 840 | 840 | 4 |
| rank開始 threshold=1e-11 | 7 | 15 | 840 | 840 | 4 |
| rank開始 threshold=1e-12 | 7 | 15 | 840 | 840 | 4 |

候補範囲・異なるstep数が広い群は、窓、R²、区間持続長、終了条件、終了持続長、rank閾値、rank保存レコード持続長への依存を人間が確認すべき群である。ここでも許容範囲や代表値は設定していない。

## 9. 解析不能または曖昧な範囲

- 各回帰窓の半窓より端側では中心回帰を定義できず、該当列はNaNである。
- `fcurve` はstep 21167で終了するため、それ以後の成長slope・成長終了候補は解析できない。
- q系列は可変間隔で保存されている。連続保存レコード間の未観測stepでrank=4が維持されたかは、この入力だけでは判定できない。
- rank持続長を「保存レコード数」ではなく「物理step数」とする最終仕様は未決定である。
- `paper7_long_timeseries` は25 step間隔であり、crossing=1167そのものは保存点ではない。保存点のflag整合だけを検査した。
- 成長区間候補が複数ある系列について、どの区間を主区間と呼ぶかは未決定である。
- 成長終了候補が見つからないパラメータ組合せも削除せず `not_found` として保存した。
- 成長終了候補とrank開始候補の対応付け規則がないため、全直積だけを保存した。

## 10. q比と既存rank定義の検証

- 既存定義: `rank_Q = count(q_j > 1e-8 q1), j=1,...,4`
- 既存`rank_q`との不一致行数: `0`
- 比較用rank閾値は頑健性確認だけであり、既存定義を変更していない。
- `q3/q1`, `q4/q1`, `min(q3,q4)/q1`, `q3-q4`, `q3/q4` を保存した。
- ゼロ除算はNaNであり、比に固定イベント閾値を置いていない。

## 11. 出力表

| 表 | CSV | Markdown |
|---|---|---|
| `growth_intervals_all_candidates` | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/processed/growth_intervals_all_candidates.csv` (921.01 KiB) | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/processed/growth_intervals_all_candidates.md` (1.04 MiB) |
| `growth_end_all_candidates` | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/processed/growth_end_all_candidates.csv` (9.25 MiB) | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/processed/growth_end_all_candidates.md` (11.42 MiB) |
| `rank4_onset_all_candidates` | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/processed/rank4_onset_all_candidates.csv` (5.17 KiB) | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/processed/rank4_onset_all_candidates.md` (6.27 KiB) |
| `growth_end_vs_rank4_onset_all_pairs` | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/processed/growth_end_vs_rank4_onset_all_pairs.csv` (203.02 MiB) | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/processed/growth_end_vs_rank4_onset_all_pairs.md` (261.53 MiB) |
| `candidate_summary` | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/processed/candidate_summary.csv` (5.68 KiB) | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/processed/candidate_summary.md` (7.05 KiB) |

候補数が多い表も削除・間引きしていない。全ペア表は候補IDを使って正規化し、個別パラメータは成長終了表とrank開始表から参照できる。

## 12. 出力図

- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/figures/figure01_f_and_log10.png`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/figures/figure01_f_and_log10.svg`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/figures/figure02_regression_slopes.png`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/figures/figure02_regression_slopes.svg`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/figures/figure03_regression_r2.png`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/figures/figure03_regression_r2.svg`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/figures/figure04_growth_intervals_by_window.png`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/figures/figure04_growth_intervals_by_window.svg`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/figures/figure05_growth_end_candidate_distribution.png`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/figures/figure05_growth_end_candidate_distribution.svg`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/figures/figure06_q1_q4_and_rank_q.png`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/figures/figure06_q1_q4_and_rank_q.svg`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/figures/figure07_q_ratios.png`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/figures/figure07_q_ratios.svg`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/figures/figure08_rank4_onset_heatmap.png`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/figures/figure08_rank4_onset_heatmap.svg`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/figures/figure09_growth_end_vs_rank4_difference.png`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/figures/figure09_growth_end_vs_rank4_difference.svg`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/figures/figure10_crossing_1167_zoom.png`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_event_calibration/figures/figure10_crossing_1167_zoom.svg`

全図は候補を同等に表示し、採用候補の色・線幅・注釈による強調を行っていない。

## 13. 人間が決定すべき事項

1. 指数成長回帰に採用する窓幅、R²閾値、最小持続長。
2. 複数の連続成長区間がある場合の対象区間。
3. 成長終了条件A/B/Cと終了持続長。
4. `rank_Q=4` 頑健性確認で重視する相対閾値。
5. rank持続を保存レコード数で定義するか、別の毎stepデータを要求するか。
6. 成長終了候補とrank開始候補を対応付ける規則。
7. 既存crossing=1167と新たに校正するイベント群の関係。

上記事項が承認されるまで、単一のイベント時刻を固定してはならない。

## 14. 実行環境と停止

- Python: `3.9.6 (default, Jan  9 2026, 11:03:41) 
[Clang 17.0.0 (clang-1700.6.4.2)]`
- NumPy: `2.0.2`
- OS: `macOS-26.3.1-arm64-arm-64bit`
- 報告生成日時（UTC）: `2026-07-26T05:37:30.076845+00:00`
- 報告生成時間（秒）: `0.528008`
- **最終状態: `CALIBRATION_DATA_COMPLETE`**
- Stage A1はここで停止する。Stage A2へ進まない。
