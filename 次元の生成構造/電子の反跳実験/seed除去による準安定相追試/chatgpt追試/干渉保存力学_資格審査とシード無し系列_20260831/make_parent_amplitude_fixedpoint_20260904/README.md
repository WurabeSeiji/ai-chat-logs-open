# make_parent amplitude fixed-point analysis (2026-09-04)

このフォルダは、`make_parent` の段階1が選ぶ振幅スケールを独立に再現・監査するための専用環境です。

## 原則
- `source_snapshot/original_engine.py` は Google Drive 上の原本スナップショットを無変更で保存。
- 原本コードは変更しない。
- 診断は `program/analyze_stage1_amplitude.py` から外部再現する。
- seed、反復数、beta を固定して再現可能にする。

## 構成
- `source_snapshot/`: 原本コードのスナップショット
- `program/`: 診断プログラム
- `data/input_saved_parents/`: N=3..6 の保存済み parent_v.csv（比較用）
- `results/`: 再現結果、反復末尾トレース、exact候補、分析メモ
- `MANIFEST.sha256`: 全ファイル SHA256

## 実行
```bash
python3 program/analyze_stage1_amplitude.py
```

## 主要出力
- `results/stage1_amplitude_N3_N16.csv`
- `results/stage1_last20_iterations_N3_N16.csv`
- `results/exact_candidates_N3_N6.json`
- `results/analysis_summary.md`

## 注意
現行 `_eigenmode_residual` は `v` の振幅正規化削除以前の式を残しているため、非正規化 `v` に対して Rayleigh quotient の分母を欠いています。本環境は原本を修正せず、現行残差とスケール不変残差を並記します。
