# Reproduction Package: 二段階seed除去による三方向生成の因果分離（第8論文）

論文「N体関係波閉鎖系における三方向生成の時間構造——二段階seed除去による因果分離」の再現パッケージ。

## 構成

```
programs/
  originals_paper7/                第7論文原本（read-only 再利用、無変更）
  preliminary_seed_ablation/       予備実験: 条件A/B/D 全N比較（表2, 図1, 図7）
  stage_A0_paper7_reproduction/    第7論文 N=5 軌道のビット一致再現
  stage_A1b_transition_anatomy/    移行解剖: rank_q/占有/f の時刻分離（図5, 表4）
  stage_A2a_seedless_N5/           完全無seed N=5 毎step観測（図2, 図3, 図4）
  stage_A2c_direction_lineage_N5/  seedあり方向系譜（§7.2, §7.3）
  stage_A2d_seedless_direction_lineage_N5/  無seed方向系譜（表3, 図6）
  seedless_natural_and_long_horizon/  無seed自然軌道・t=110000 長時間（図8, 図9）
  supporting_analysis/             単調性検査・親状態構造・残余占有・N=5長時間判定（§4.1, 表1, 表5）
figures/                           論文掲載図 9 点
results/
  stage_reports/                   各実験の観測報告書（数値のみ・解釈なし）
  processed_tables/                first-passage・decade増幅率・overlap 等の表
  data_summaries/                  実行メタデータ・SHA-256・対照結果 JSON
```

## 図の対応

| 論文図 | ファイル |
|:--|:--|
| 図1 | figures/fig01_f_compare_N00005.png |
| 図2 | figures/figure02_seedless_log10_f.png |
| 図3 | figures/figure05_decade_growth_rate_comparison.png |
| 図4 | figures/figure03_seeded_vs_seedless_absolute_step.png |
| 図5 | figures/figure06_q_ratios_and_rank_q_0_3000.png |
| 図6 | figures/figure01_seedless_early_vs_late_lineage.png |
| 図7 | figures/fig04_metastable_B_vs_D_N00005.png |
| 図8 | figures/figure_long_horizon_one_exp_residual_x20.png |
| 図9 | figures/figure3_compare_N5_N40_N300.png |

## 原本の固定（SHA-256）

第7論文原本は編集せず read-only import で再利用した。各実行は実行前に SHA-256 照合を行う。

| file | sha256 |
|:--|:--|
| run_n_scaling_lowrank_v1.py | ba0fc19b03caf06d16e97e2cd5da499ed2ed8e95288f6dbf2062bedcaf11176d |
| run_plane_flow_exact_v1.py | 9cf28ca8c0d2ad8fac2f0f6dae045248695247c5809c21ccb2069ef91a94ab76 |
| run_plane_flow_approx_v1.py | a9d247a8070d849fe989e35e00320e470968354614424c9bdceda57132d9f0fa |
| run_n300_dimension_saturation_v2.py | 229938a66631057426f187ed80b17de08cfcb9107cfe509c30f5bbdcca3a03e6 |
| run_paper7_5color_timeseries.py | fe5c7cbc33437890a5f50944cbbae1594e5f647739d4955402b578f515658503 |
| run_paper7_transverse.py | ac1073bea329971de3ff4c2fd1588d926029a8502c21e8cc01f406acb86ad60b |

## 再生成手順（概要）

各 programs/ サブフォルダのスクリプトは、原本パスの SHA-256 照合 → 実行 → 報告書生成の順で動く。決定論（`numpy.random.default_rng` 固定 seed）のため、同一環境（Python 3.9.6 / NumPy 2.0.2）で全数値が再現する。主要な実行時間: 予備実験 N=300 が約 40 分、その他は各数秒〜数分、長時間対照 N=300 が約 28 分、五色分解 N=300 が約 14 分。

- 予備実験: `run_preliminary_seed_ablation_v1.py`（条件 A/B/D、N=5/40/300）
- Stage A2a: `verify_sources.py` → `run_seedless.py` → 比較・作図・報告
- 方向系譜: `replay_and_extract_bases.py` → `analyze_direction_lineage.py`（A2c）、`analyze_saved_seedless_lineage.py`（A2d）
- 長時間: `run_seedless_natural_long_horizon_v1.py run N`（N=5 は `run_long_horizon_N5_driver_v1.py`）
- 五色分解: `run_seedless_natural_figures3_4_v1.py run N` → `figures`

## データ方針

生 CSV 時系列は決定論的に再生成可能なため同梱しない（本シリーズの規約）。各 CSV の SHA-256 は results/data_summaries/ の実行メタ JSON に記録されており、再生成結果の同一性を検証できる。

## 数値健全性

規格化誤差 ≤1.4×10⁻¹⁴、零二乗閉鎖 ≤1.8×10⁻¹⁰（条件Dの注入直後、他は ≤1.3×10⁻¹³）、射影閉鎖誤差 0、占有和誤差 0。独立 2 実行のビット一致（Stage A2a）、第7論文軌道のビット一致再現（Stage A0）、N=300 五色分解の既存条件 A との共通 27 列全行一致を確認済み。
