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

## 5. 初期時刻部分空間との重なり（§4, initial_origin 修正済）

初期 t0 の全固有値の ρ_j=|σ_j|/(ε‖K‖) から床/非床部分空間を構成（ρ<1e3 を床）。t0 では
ρ が 1〜10（床）と 1e15（非床）に明瞭分離し、床空間は閾値 1e2〜1e6 で不変。
重なりは**クラスタ基底の全列**を用いた

  O(U,B) = ‖U†B‖_F² / dim(B) = Tr(Π_U Π_B)/rank(Π_B)

で計算（旧実装は先頭2列のみを使い多次元クラスタで過小評価していた。廃止）。

**重なり閉鎖（O_floor + O_nonfloor = 1）**：初期床・非床が全空間を直交分解するため理論上 1。

| 量 | 値 |
|:--|--:|
| max\|O_floor+O_nonfloor−1\| | 1.55×10⁻¹⁵ |
| median\|…−1\| | 2.22×10⁻¹⁶ |
| 誤差>10⁻¹² の件数 | 0 |
| 誤差>10⁻¹⁰ の件数 | 0 |
| 合格（<10⁻¹²） | ✓ |

**initial_origin_status（全列 Frobenius・閾値 0.99）**：

| status | 件数（クラスタ×時刻） |
|:--|--:|
| mixed | 3852 |
| initial_nonfloor | 679 |
| initial_floor | 162 |
| undetermined | 0 |

**origin 閾値感度**：

| 閾値 | initial_floor | initial_nonfloor | mixed | undetermined |
|--:|--:|--:|--:|--:|
| 0.90 | 186 | 747 | 3760 | 0 |
| 0.95 | 177 | 705 | 3811 | 0 |
| 0.99 | 162 | 679 | 3852 | 0 |
| 0.999 | 143 | 650 | 3852 | 48 |

出力：`raw/N00005/initial_space_origin.csv`（step, time, cluster_id, branch_id, cluster_dimension,
sigma_min/max/representative, overlap_with_initial_floor_space, overlap_with_initial_nonfloor_space,
origin_overlap_sum, origin_overlap_closure_error, initial_origin_status, origin_threshold,
floor_threshold_label）、`diagnostics/N00005_initial_origin_revision.json`、
`tables/N00005/initial_origin_summary.csv`。図 FigC（床/非床重なり）, FigC(closure)（和≈1）, FigD
（O_floor≥0.99 の branch）。

**旧結果（未確定として置換）**：旧 `overlap_space`（先頭2列）による initial_origin_status
（mixed 3864 / nonfloor 679 / floor 129 / undetermined 21）と旧 FigC/FigD は廃止。旧
`initial_floor_flag`（前時刻対応失敗＝床）由来の集計・図（fig09/fig10）は `deprecated/` へ隔離。

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

- `raw/N00005/branch_tracking_bijective.csv`, `cluster_lineage.csv`, `initial_space_origin.csv`
  （全列 Frobenius・新列付）, `merge_tolerance_sensitivity.csv`
- `diagnostics/N00005_branch_revision.json`, `diagnostics/N00005_initial_origin_revision.json`
- `tables/N00005/initial_origin_summary.csv`
- `figures/N00005/figA〜figF`, `figC_origin_overlap_closure.png`
- `deprecated/`：`fig09_floor_branches_obsolete.png`, `fig10_delta_branches_obsolete.png`,
  `cluster_tracking_obsolete.csv`（旧方式・解析使用禁止）

## 12. initial_origin 修正の検収（§12）

| # | 条件 | 状態 |
|--:|:--|:--|
| 1 | overlap_space が全列使用 | ✓（‖U†B‖_F²/dim(B)） |
| 2 | 初期床＋非床が全空間を構成 | ✓（U の直交補で構成） |
| 3 | 全時刻・全クラスタで \|O_floor+O_nonfloor−1\|<10⁻¹² | ✓（max 1.55×10⁻¹⁵） |
| 4 | initial_origin_status 全面再生成 | ✓ |
| 5 | 旧 initial_floor_flag を解析から排除 | ✓（deprecated 隔離） |
| 6,7 | FigC / FigD 再生成 | ✓（＋closure 図） |
| 8 | origin 閾値感度保存 | ✓（0.90/0.95/0.99/0.999） |
| 9 | 旧方式成果物を削除/隔離 | ✓（deprecated/, _obsolete） |
| 10 | README に旧方式使用禁止 | ✓ |
| 11 | 報告書の件数・表・図を修正版へ | ✓（本§5, §11） |
| 12 | 固有値スペクトル・branch追跡 不変 | ✓（1,0.251458,…／unmatched=78,総81 不変） |
| 13,14 | N=40・N=300 未着手 | ✓ |
| 15 | 物理解釈なし | ✓ |

**N=40 は人間検収後に同一コードで適用。N=300 未着手。**
