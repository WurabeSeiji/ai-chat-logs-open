# N=40 厳密全固有値・全固有方向 観測報告（観測値のみ・解釈なし）

N=5 で検収済みの `exact_lowN_eigenspectrum_v2` を**実装変更なし**で N=40 へ適用した結果。
本書は数値のみ。物理解釈・N依存の断定・1/N則・時間軸/Q軸/フェルミオン生成・N=300外挿・論文6評価は
書かない（実行指示 §8）。**N=300 未着手。**

## 1. 計算条件（N=5 と同一）

N=40, M=780, seed=40300722, δ=1e-15。crossing（真, f>0.05 最初）= step **2011**。観測終端 52011。
記録時刻数 1164。固有分解 `eigh(1j*K)`。併合閾値 1e-12、回転/核相対床 σ/(ε‖K‖)。
軌道生成エンジン・seed規則・零閉鎖・正規化・Cayley更新・演算順序・crossing定義・crossing後50000step・
時刻サンプリング・branch追跡(Hungarian一対一)・lineage・初期床判定・縮退クラスタ処理・origin overlap・
数値診断・図生成方法は N=5 と同一。f は原本走行とビット一致（`max_f_ref_dev`=0.0）。

## 2. 数値診断（全時刻最大）

| 量 | N=40 |
|:--|--:|
| \|K+Kᵀ\| | 0.0 |
| 固有対残差 \|Hu−μu\| | 2.4×10⁻¹³ |
| \|U†U−I\| | 8.5×10⁻¹⁴ |
| 正負対誤差 | 9.2×10⁻¹⁴ |
| 非縮退平面間直交 max‖B_C^T B_D‖₂ | 2.9×10⁻¹⁵（≤10⁻¹²） |
| 射影冪等 \|Π²−Π\| | 4.3×10⁻¹⁸ |
| 状態分解閉鎖 \|\|Z\|²−ΣE_C−E_ker\| | 5.7×10⁻¹⁵ |
| origin 重なり閉鎖 max\|O_floor+O_nonfloor−1\| | 7.3×10⁻¹⁵（≤10⁻¹²） |
| f 参照一致 | 0.0 |

branch ID 同時刻重複＝0（Hungarian 一対一, assert 通過）。

## 3. 主判定（実行指示 §7）

### 3.1 新しい平面（q₃, q₄）

$Q(t)=[B_0\mid B_{\mathrm{dom}}(t)]$ の特異値：

| 時刻 | q₁ | q₂ | q₃ | q₄ |
|:--|--:|--:|--:|--:|
| 最終 (52011) | 1.3785 | 1.3728 | **0.3396** | **0.3159** |
| 準安定域 (7011) | — | — | 0.3382 | 0.3114 |

$q_3>0,\ q_4>0$。

### 3.2 独立方向数

数値閾値 10⁻⁸ で $q_1,q_2,q_3,q_4$ すべて > 10⁻⁸ ⇒ **rank[B_0∣B_dom] = 4**。

### 3.3 瞬時支配平面の初期空間成分（最終時刻）

| 量 | 値 |
|:--|--:|
| O_floor（初期床空間） | 0.1341 |
| O_nonfloor（初期非床空間） | 0.8659 |
| O_floor + O_nonfloor | 1.000000（誤差 ≤7.3×10⁻¹⁵） |

### 3.4 状態占有（最終時刻・準安定域）

| 量 | 値 |
|:--|--:|
| 瞬時支配平面占有 E_dom/\|Z\|² | 1.000000 |
| 最大非支配クラスタ占有 | 8.27×10⁻¹⁸ |
| 瞬時核占有 | 閉鎖より 10⁻¹⁵ 規模 |
| 固定初期親基底で見た支配平面外成分（f=1−E_P1, 準安定代表） | ≈0.20（第6論文 §5.4 と一致） |

固定初期基底の結果（f, 平面別ノルム）と瞬時固有基底の結果（E_C, q, O_floor）は別量として区別。

## 4. 全固有値・スペクトル

全時刻・全 M=780 固有値を閾値削除なしで保存（`raw/N00040/eigenvalues.csv`）。σ_j/σ_1, Nσ_j/σ_1,
N²σ_j/σ_1, σ/(ε‖K‖) を保存。最終時刻の σ_j/σ_1 は 1.000000 と、0.4622〜0.4872 の 39 枝の帯。
図：fig01（比 線形）, fig02（対数）, fig01C（ヒートマップ）, fig01D（branch追跡）, fig03（N倍）,
fig04（N²倍）, fig08（全代表時刻スペクトル）。

## 5. branch・lineage・origin・感度

| 量 | 値 |
|:--|--:|
| 同時刻 branch ID 重複 | 0 |
| unmatched（受入<0.7 で新規） | 108 |
| 総 branch 数 | 148 |
| lineage: split_candidate / merge_candidate / ambiguous | 4 / 4 / 172 |
| 併合感度 max_intercluster_overlap（1e-10〜1e-14） | 全て 2.2×10⁻¹⁵（安定） |
| origin status（0.99）: floor / nonfloor / mixed / undetermined | 0 / 14009 / 32591 / 0 |

crossing 真=2011、代表=2011（真 crossing をサンプルに追加）。

## 6. N=5 と N=40 の比較（実行指示 §7.5）

| 量 | N=5 | N=40 |
|:--|--:|--:|
| q₃ | 0.8140 | 0.3396 |
| q₄ | 0.6893 | 0.3159 |
| rank[B₀∣B_dom]（閾値10⁻⁸） | 4 | 4 |
| 支配平面 O_floor | 0.4721 | 0.1341 |
| 支配平面 O_nonfloor | 0.5279 | 0.8659 |
| 支配占有 E_dom/\|Z\|² | 1.000000 | 1.000000 |
| 最大非支配占有 | 4.55×10⁻¹⁴ | 8.27×10⁻¹⁸ |
| 非支配 σ_j/σ_1 帯 | 0.2485〜0.2515 | 0.4622〜0.4872 |

比較図：`figures/comparison/fig18_q3q4_compare.png`, `fig19_dominant_origin_compare.png`,
`fig20_dominant_occupation_compare.png`, `comparison_table_N5_N40.json`。
（N依存の傾向・法則は数値確認前に断定しない。N=300 未着手。）

## 7. 図（実行指示 §5, 1〜20）

figures/N00040/：fig01A/B（σ比 線形/対数）, fig01C（ヒートマップ）, fig01D（branch）, fig02A/2B相当
（N/N²）, fig03（δ）, fig05（床比）, fig06（正負対）, fig07（縮退クラスタ）, fig08（q）, fig12（占有）,
fig13（平面間直交）, figA（一対一branch）, figB（lineage）, figC（初期床/非床）, figC(closure),
figD（初期床由来）, figE（併合感度）, figF（対応重なり）。
figures/comparison/：fig18（q3q4）, fig19（支配origin）, fig20（支配占有）。
（本 v2 の図番号は N=5 と同一形式。実行指示 §5 の 20 図に対応。）

## 8. 完了条件（実行指示 §9）

| # | 条件 | 状態 |
|--:|:--|:--|
| 1 | 全時刻で完全固有分解 | ✓（eigh(iK), 1164時刻） |
| 2 | 全固有値保存 | ✓（閾値削除なし, 780×1164） |
| 3 | branch重複ゼロ | ✓（0, assert） |
| 4 | lineage保存 | ✓ |
| 5 | origin closure誤差<10⁻¹² | ✓（7.3×10⁻¹⁵） |
| 6 | 状態分解閉鎖 倍精度 | ✓（5.7×10⁻¹⁵） |
| 7 | q₁〜q₄ 全時系列保存 | ✓ |
| 8 | 代表時刻の全モード表 | ✓（6時刻） |
| 9 | 必須図生成 | ✓ |
| 10 | N=5との比較表 | ✓ |
| 11 | N=300 未着手 | ✓ |
| 12 | 物理解釈なし | ✓ |

## 9. 出力

- `raw/N00040/`：eigenvalues.csv（全固有値）, clusters.csv, delta_targets.csv, q_svd.csv,
  branch_tracking_bijective.csv, cluster_lineage.csv, initial_space_origin.csv,
  merge_tolerance_sensitivity.csv, diagnostics_timeseries.csv（*.csv は Drive 上, git 追跡外）
- `binary/N00040/`：timeseries.npz, decomp_*.npz（代表6時刻）
- `tables/N00040/`：fulltable_*, clusters_*, initial_origin_summary.csv
- `diagnostics/`：N00040.json, N00040_branch_revision.json, N00040_initial_origin_revision.json
- `figures/N00040/`：fig01〜fig13, figA〜figF ；`figures/comparison/`：fig18〜20 + 比較表
