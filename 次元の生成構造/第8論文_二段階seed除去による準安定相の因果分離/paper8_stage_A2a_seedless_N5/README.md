# Stage A2a: N=5、明示seedなしfloat64自然床系列

このパッケージは、既存の第7論文と同じ `N=5` / `float64` 力学を用い、親状態生成後の初期状態を厳密に `Z0 = v` とした2回の決定論的実行を行う。第8論文の仮説判定、下位階層の認定、高精度化、Delta掃引、`N=40`、`N=300` は行わない。

## 固定事項

- 親: `make_parent(sys, np.random.default_rng(40265722), iters=1200, tol=1e-12)`
- 初期状態: `Z0 = v.copy()`（明示的摂動なし、再正規化なし）
- 実行: `step=0..5000`、同じPRNG seedで2回
- 禁止された `zero_closure_kernel_seed` は呼ばない
- warm-start用実ベクトル `wp` は親生成後に同じPRNGから生成する。これは状態 `Z` に加算されず、Cayley更新の `sigma_max_power` にだけ渡される。
- q測定は5 stepごと、状態占有は25 stepごと。指定f水準の初回通過stepでは追加測定する。
- 測定値および支配平面は、状態更新にも `sigma_max_power` のwarm-startにも戻さない。
- seedあり参照に実保存値がないstepでは、Stage A1bの直前・直後の実保存値を併記する。時間補間は行わない。

## 実行順序

```bash
python3 verify_sources.py
python3 run_seedless.py
python3 compare_execs.py
python3 compare_with_seeded_reference.py
python3 make_figures.py
python3 make_report.py
```

各スクリプトは前工程の成功記録を確認する。原本・Stage A0・Stage A1b成果物は読み取り専用であり、このディレクトリ外へ書き込まない。

## 出力

- `raw/`: 2実行の毎step f、q、状態占有、初回到達、支配平面
- `processed/`: 必須比較表、時間平行移動比較、数値健全性
- `figures/`: 図1〜14のPNG/SVG
- `reports/stage_A2a_seedless_N5_report.md`: 最終報告書
- `logs/`: 工程ログと固定入力検証記録

最終報告後は停止し、後続Stageへ進まない。
