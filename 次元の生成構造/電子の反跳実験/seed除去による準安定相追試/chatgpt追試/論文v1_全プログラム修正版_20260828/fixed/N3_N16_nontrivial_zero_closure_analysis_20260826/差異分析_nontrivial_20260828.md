# 差異分析——N3_N16_nontrivial_zero_closure_analysis_20260826（原本 / 修正版 / baseline）

**作成日:** 2026-08-28　全表は `results/three_way_comparison.md`。

## 適用した修正（全パッケージ共通、判断 1〜6）

- A1 親の振幅正規化 `v/‖v‖` を除去（親の固有モード反復は位相のみ K のまま＝S3）
- A2/A3 初期状態は親そのもの `Z = v`（外部 seed δg も正規化も無し）
- A4 動力学の生成子を振幅込み `K_ij = Im(z̄ᵢ zⱼ)`（`set_state`）。baseline は `KMODE=phase` で位相のみ K
- A5 `zero_closure_generic` の正規化除去（未使用）／A6(b) σ は実際に回している K の厳密スペクトルで読む（冪反復は実行しない）
- R1 Cayley → `exp((2π/144)·K)`（`linear_rotation_step`、`cayley_step` は定義のみ残す）／R2 K/σ 枝廃止／R3 validate 書換え、γ は dphi=2π/n_den
- S1 `zero_closure_kernel_seed` を呼ばない／S2 `sigma_max_power`・`wp=rng.normal` を呼ばない（定義は残す）

**既知の問題（親探索）**：A1 により `_eigenmode_residual` が単位ノルム前提のため残差が意味を失い、収束判定が効かず、3 リスタートの選択も無効。N=6, 7, 10, 11 では**収束していない親**（真残差 0.47 / 0.39 / 0.081 / 1.4e-3）が選ばれている。summary.json の `parent_residual`（0.49〜73）は意味のない数。残差判定のスケール不変化は次回の判断項目。


## 読み取り

- 同梱スクリプトは `reproduced_small_subset_classification.csv` を書くのみ。評価表（B-5）は生成プログラム未同梱。SOURCE は修正版の最終状態に差替え済。

## 再生成した図

`N3_N16_trivial_vs_nontrivial_small_closures.png`, `N5_nontrivial_pair_closure_graph.png`, `N5_vs_N14_nontrivial_closure_time_evolution.png`

## このパッケージへの変更箇所

```
（パス変更・入力差替えのみ）
```
