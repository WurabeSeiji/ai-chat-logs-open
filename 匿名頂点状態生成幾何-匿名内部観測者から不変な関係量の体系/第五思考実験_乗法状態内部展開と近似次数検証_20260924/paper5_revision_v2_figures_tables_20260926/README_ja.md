# 第五思考実験・改訂版 図表再生成パッケージ v2

日付: 2026-09-26

目的: 第五論文の主張を「見かけの2状態 → 5永続状態 → 外部νの発見 → 6永続状態」へ更新し、保存済み監査データだけから論文用図表を再生成する。

**このパッケージは物理計算を再実行しない。** 既存の監査 CSV/JSON を入力として、図表・索引・SHA256 を生成する。

## 実行

```bash
python scripts/generate_paper5_revision_v2.py
```

依存: Python 3.11+, numpy, pandas, matplotlib, tabulate。

## 出力

- `figures/figure01_state_audit_2_to_5_to_6_v2.svg`
- `figures/figure02_state_closure_sync_scope_v2.svg`
- `figures/figure03_methodB_strict_12case_heatmap_v2.svg`
- `figures/figure04_angular_order_convergence_v2.svg`
- `figures/figure05_integration_order_convergence_v2.svg`
- `figures/figure06_integration_order_9case_summary_v2.svg`
- `tables/table01_state_audit_summary_v2.csv`
- `tables/table02_approximation_and_state_scope_v2.csv`
- `tables/table03_key_numerical_results_v2.csv`
- `tables/tables_for_paper_v2_ja.md`
- `FIGURE_INDEX_v2_ja.md`
- `EXPERIMENT_INDEX_v2.csv`
- `SOURCE_PROVENANCE_v2.json`
- `SHA256SUMS_v2.txt`

## 主張範囲

本パッケージが反映するのは、状態閉包・同期更新・読み出し非帰還・νの第6永続状態への内部化である。6状態の数学的最小性、11相機械の自然界での必然性、既存B方式全体の後付け最厳格「完全乗法内部演算」適合は主張しない。
