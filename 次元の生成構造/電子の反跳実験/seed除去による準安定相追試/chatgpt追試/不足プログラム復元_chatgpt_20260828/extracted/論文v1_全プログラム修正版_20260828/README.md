# 論文v1 全プログラム修正版 2026-08-28

目的: `論文v1_全再現テスト_20260828/results/comparison.json` で NOT_REGENERATED とされた出力のうち、依頼文 B-1〜B-7 を再生成できるプログラムを回収・再構成する。

**既存出力の再送ではなく、生成スクリプトを収録している。** 旧出力はこの ZIP に同梱していない。

## 構成

- `recovered_original/`: Drive / package に残っていた原本。基礎 N-series generator と engine、N5/N16 provenance script。
- `B1_N5/`: N5 4群・a/b/a2/b2/ab・milestone・7図の等価再構成。
- `B2_N16/`: N16 phase structure 4ファイルの厳密 postprocessor。
- `B3_N3_N4/`: N3/N4 1/3 phase ordering CSV + 6図。
- `B4_comparison/`: 6組の comparison summary generator。
- `B5_nontrivial/`: small-subset classification, MITM, assessment/summary, N5 exact covers/time evolution, 3図。
- `B6_all_cardinalities/`: exact/MITM/annealing search, aggregator, time evolution/covers, plots。旧 anneal steps/seed 未回収を明示。
- `B7_partial/`: `zero_subset_exact_covers` の alias 修正。
- `copy_A_aliases.py`: SOURCE CSV 等のコピー専用処理。
- `run_all_reconstructed.py`: 全体 orchestrator。
- `VALIDATION.md`: 旧出力との比較結果。
- `coverage_manifest.csv`: 依頼された各出力と生成器の対応。

## 推奨実行順

展開済みの original package 群を同一 root に置いた上で:

```bash
python run_all_reconstructed.py --root /path/to/extracted_packages \
  --decompact-results /path/to/decompactification/results \
  --n5-raw-k-source /path/to/N5_raw_K_raw_observables.csv
```

この既定実行では B6 の k>=7 stochastic search は実行しない。理由は、歴史的 `steps,seed` を回収できなかったため。

新しい明示的な固定値で今後の再現実験を開始する場合のみ:

```bash
python run_all_reconstructed.py --root /path/to/extracted_packages \
  --decompact-results /path/to/decompactification/results \
  --n5-raw-k-source /path/to/N5_raw_K_raw_observables.csv \
  --run-new-stochastic-search
```

この場合の引数は `B6_all_cardinalities/reproduction_anneal_args.csv`。**historical values ではない。**

## 依存

Python: numpy, pandas, matplotlib, scipy。B5/B6 の exact/annealing search には C++17 compiler (`g++`) が必要。

## K/sigma 枝

依頼どおり B-8 は再構成しない。
