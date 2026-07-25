# N=5 branch追跡・初期床判定・縮退感度 修正版 観測報告（観測値のみ・解釈なし）

`exact_lowN_eigenspectrum_v2_branch修正指示.md` への対応。固有分解本体（eigh(iK)・全固有値・
瞬時クラスタ）は維持し、**後処理のみ**修正した。**N=5 のみ。N=40・N=300 未着手。物理解釈なし。**

## 1. 修正点

1. branch追跡を **Hungarian 一対一**（`scipy.optimize.linear_sum_assignment`）へ。各現在クラスタが
   独立に最大重なり過去を取る方式を廃止。同時刻の branch ID 重複を assert で禁止。
2. 分裂・合流を **lineage edge** として別保存（continuation / split_candidate / merge_candidate /
   birth / death / ambiguous）。
3. `initial_floor_flag`（旧：前時刻対応失敗＝床）を廃止。**初期時刻部分空間との直接重なり**で
   `initial_origin_status`（initial_nonfloor / initial_floor / mixed / undetermined）を定義。
4. 併合閾値 **1e-10〜1e-14** の感度表。
5. crossing 表記：真 crossing=step1167 を区別（v2 主表は 5step 格子のため直前点 step1165 を使用）。

## 2. 計算条件

N=5, M=10, seed=40265722, δ=1e-15。crossing（真）= step **1167**。観測終端 51167。
主分解 merge_tol=1e-12、eigh(1j*K)。branch受入主閾値 0.7（診断併記 0.5/0.7/0.9/0.99）、
lineage 重なり閾値 0.3、初期床 ρ=σ/(ε‖K‖)<1e3（併記 1e2/1e4/1e6）。

## 3. branch 一対一追跡（§2）

| 量 | 値 |
|:--|--:|
| 同時刻 branch ID 重複 | **0**（全時刻 assert 通過） |
| 受入(≥0.7)なしで新規付与 (unmatched) | 78 |
| 総 branch 数 | 81 |

最終スペクトル σ_j/σ_1 = **1.000000, 0.251458, 0.250026, 0.249974, 0.248533**（v2 と一致）。
出力：`raw/N00005/branch_tracking_bijective.csv`（source/target step・cluster_id・branch_id・
overlap・assignment_cost・accepted・tracking_status）。図 FigA（同時刻 branch 重複なし）、
FigF（採用対応の最小重なり）。

## 4. lineage（§3）

| relation_type | 総数 |
|:--|--:|
| split_candidate | 8 |
| merge_candidate | 6 |
| ambiguous | 134 |
| continuation / birth / death | `cluster_lineage.csv` に全保存 |

判定規則：source が閾値以上の target を複数持てば split_candidate、target が複数 source を持てば
merge_candidate、両方なら ambiguous、対応 target なしは death、対応 source なしは birth。
これらは幾何的追跡ラベルであり物理解釈を付けない。出力：`raw/N00005/cluster_lineage.csv`。図 FigB。

## 5. 初期時刻部分空間との重なり（§4）

初期 t0 の全固有値の ρ_j=|σ_j|/(ε‖K‖) から床/非床部分空間を構成（ρ<1e3 を床）。t0 では
ρ が 1〜10（床）と 1e15（非床）に明瞭分離し、床空間は閾値 1e2〜1e6 で不変。各クラスタ×時刻の
初期床/非床空間重なりと `initial_origin_status` を保存。

| initial_origin_status | 件数（クラスタ×時刻） |
|:--|--:|
| mixed | 3864 |
| initial_nonfloor | 679 |
| initial_floor | 129 |
| undetermined | 21 |

出力：`raw/N00005/initial_space_origin.csv`（overlap_with_initial_floor_space,
overlap_with_initial_nonfloor_space, initial_origin_status, floor_threshold_label）。図 FigC, FigD。

## 6. 併合閾値感度（§5）

代表時刻（initial / nearest_before_crossing_step1167 / true_crossing_step1167 / post_crossing /
plateau_start / final）で merge_tol ∈ {1e-10,…,1e-14} を比較。

| 量 | 挙動（全代表時刻・全閾値） |
|:--|:--|
| cluster_count | initial=3、他=5 で**閾値不変** |
| dominant_delta | 各時刻で閾値不変 |
| dominant_occupation | 各時刻で閾値不変 |
| q1..q4 | 閾値不変 |
| max_intercluster_overlap | 全代表時刻・全閾値で ≤1.9×10⁻¹⁵ |
| closure_error | 倍精度範囲 |

代表時刻では主要観測量が 1e-10〜1e-14 で**定性・定量とも安定**（FigE は全量平坦）。
出力：`raw/N00005/merge_tolerance_sensitivity.csv`。図 FigE。
（近縮退による非縮退平面間重なりのスパイクは crossing 近傍の非代表時刻で生じ、主分解は §11 の
diagnostics_timeseries に既出。）

## 7. crossing 表記（§6）

実際の crossing は step=**1167**。v2 主表（`run_exact_lowN_eigenspectrum_v2`）は 5step 格子のため
直前点 step=1165 を 'crossing' 表記に用いていた。本修正では真 crossing step=1167 をサンプルに追加し、
`representative_crossing_step`=1167（=真 crossing）で後処理した。`diagnostics/N00005_branch_revision.json`
に `crossing_step` と `representative_crossing_step` を明記。

## 8. 図（§8, A〜F）

figures/N00005/：figA（一対一 branch σ/σ1・同時刻重複なし）、figB（lineage 件数）、
figC（初期床/非床空間重なり branch別）、figD（初期床重なり≥0.99 の branch の σ/σ1）、
figE（併合閾値感度）、figF（採用 branch 対応の最小重なり）。生データ CSV から再生成可能
（`code/run_branch_revision_v2.py 5`）。

## 9. 再検収（§11）状況

| # | 条件 | 状態 |
|--:|:--|:--|
| 1 | 同時刻 branch ID 重複ゼロ | ✓（0, assert） |
| 2 | Hungarian 一対一 | ✓ |
| 3 | 分裂・合流 lineage を別保存 | ✓ |
| 4 | 初期床を初期時刻部分空間との直接比較で定義 | ✓ |
| 5 | 旧 initial_floor_flag 定義を廃止 | ✓ |
| 6 | 併合閾値 1e-10〜1e-14 感度表 | ✓ |
| 7 | 主要観測量の閾値安定性 | ✓（代表時刻で全量不変, FigE） |
| 8 | crossing と近傍の表記区別 | ✓ |
| 9 | branch別全図が修正版生データから再生成可 | ✓ |
| 10 | 固有値スペクトルが v2 と一致 | ✓（1, 0.251458, 0.250026, 0.249974, 0.248533） |
| 11 | 数値診断が倍精度範囲 | ✓ |
| 12 | 観測報告に物理解釈なし | ✓ |

## 10. 未確定を解消した項目

v2 で branch ID を用いた図・branch別δ・初期床由来 branch・birth/death・分裂合流・
`initial_floor_flag` 集計・最終表 branch 列は、本修正で Hungarian 一対一 + lineage + 初期空間重なりへ
置換した。維持項目（全固有値スペクトル・正負対検証・固有対残差・瞬時クラスタ生データ・瞬時占有・
親平面との瞬時比較・q特異値・全体数値診断）は不変。

## 11. 出力一覧

- `raw/N00005/branch_tracking_bijective.csv`, `cluster_lineage.csv`, `initial_space_origin.csv`,
  `merge_tolerance_sensitivity.csv`
- `diagnostics/N00005_branch_revision.json`
- `figures/N00005/figA〜figF`

**N=40 は人間検収後に同一コードで適用。N=300 未着手。**
