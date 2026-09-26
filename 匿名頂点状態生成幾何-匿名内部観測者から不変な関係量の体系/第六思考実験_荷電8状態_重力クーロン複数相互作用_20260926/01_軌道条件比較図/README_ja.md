# 論文6・5条件軌道比較図 v1

## 目的
論文6で使用する荷電二体系5条件について、条件による軌道差を同一スケールで比較する。
新しい物理実験は行わず、既存の解析基準式と、保存済み厳格8状態実験の個別軌道図を用いる。

## 生成図
1. `figure01_five_case_reference_orbits_same_scale.svg`
   - 5条件の独立解析基準軌道。
   - 全パネルで x,y 範囲を同一に固定。
   - `r=50 -> 20` の全軌道。
2. `figure02_five_case_strict_generator_overlay_panels.svg`
   - 各条件の保存済み `analytic reference vs strict 8-state generator` 軌道図を5パネルへ統合。
   - 元SVGは変更しない。
3. `figure03_five_case_radius_vs_cycles.svg`
   - 軌道半径 r と累積周回数の関係。軌道差の定量的補助図。
4. `figure04_five_case_first_10_orbits_same_scale.svg`
   - 各条件の最初の10周を同一空間スケールで比較。
5. `paper6_orbit_condition_case_summary.csv`
   - lambda_A, lambda_B, C, D, A, B, r=50->20 の周回数と時間。

## 5条件
- n1 attractive: lambda_A=+0.30, lambda_B=-0.30
- n1 repulsive: lambda_A=+0.30, lambda_B=+0.30
- n3 attractive: lambda_A=+0.90, lambda_B=-0.90
- n3 repulsive: lambda_A=+0.90, lambda_B=+0.90
- n4 attractive: lambda_A=+1.20, lambda_B=-1.20

## 再現
`generate_paper6_orbit_condition_comparison_v1.py` を使用する。

解析基準図だけなら外部入力は不要。厳格生成器 overlay 統合図を作る場合は、既存5ケースの `figure01_orbit_overlay.svg` をローカルへコピーし、次のファイル名で `strict_sources/` に配置する。

- `n1_attractive_figure01_orbit_overlay.svg`
- `n1_repulsive_figure01_orbit_overlay.svg`
- `n3_attractive_figure01_orbit_overlay.svg`
- `n3_repulsive_figure01_orbit_overlay.svg`
- `n4_attractive_figure01_orbit_overlay.svg`

実行例:

```bash
python generate_paper6_orbit_condition_comparison_v1.py \
  --outdir figures \
  --strict-svg-dir strict_sources
```

解析側は既存の `07_荷電二体系_解析近似基準解_20260925` と同じ LO 閉形式を使用し、厳格8状態生成器へは一切フィードバックしない。
