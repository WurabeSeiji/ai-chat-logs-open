# Stage A0: 第7論文 N=5 再現パッケージ

## 目的

このパッケージは、第7論文の既存 N=5 結果を既存コード・既存定義のまま再生成し、保存済み成果物と比較するためだけの隔離ラッパーである。第8論文の新実験ではない。

次のものは実装しない。

- 新しい物理量またはイベント定義
- `q5`〜`q8`
- `flux`
- `eta_noise`
- 方向連続性指標
- 80/128/256 bit 高精度計算
- Series 1〜3
- N=40、N=300

既存原本コードと既存成果物は読み取り専用で扱う。全再生成物は `reproduced/`、比較結果は `comparison/`、ログは `logs/`、報告書は `reports/` の下にだけ保存する。

## 固定された実行順序

リポジトリのルートから、次の順序で実行する。

```bash
PYTHONDONTWRITEBYTECODE=1 python3 次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_reproduction/verify_sources.py
PYTHONDONTWRITEBYTECODE=1 python3 次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_reproduction/run_reproduction.py
PYTHONDONTWRITEBYTECODE=1 python3 次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_reproduction/compare_outputs.py
PYTHONDONTWRITEBYTECODE=1 python3 次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_reproduction/make_reproduction_report.py
```

各工程は直前工程の成功記録を検査する。失敗時は停止し、出力を削除・修正・再計算しない。

## N=5 固定

- すべての計算関数へ渡す `n` は整数 `5` だけである。
- `make_paper7_figures.py` と `make_saturation_comparison.py` の原本は編集せず、ラッパー内でモジュール変数 `NS` を `[5]` に固定してから呼ぶ。
- N=40、N=300 の既存成果物は読まない。

## 比較規則

CSVでは列名・行数・整数列・step列・NaN位置を検査し、浮動小数列ごとに最大絶対誤差と最大相対誤差を出す。相対誤差は、基準値が非ゼロなら `abs(reproduced - expected) / abs(expected)`、両方ゼロなら `0`、基準値だけゼロなら無限大とする。これは再現比較上の計算規約であり、物理イベント定義ではない。

- `exact`: 整数・文字列・step列が完全一致し、浮動小数値が IEEE 754 binary64 として bitwise 一致する。
- `numerical_match`: 列名・行数・整数・文字列・NaN位置が一致し、全浮動小数列で最大絶対誤差 `<= 1e-12` かつ最大相対誤差 `<= 1e-10`。
- `mismatch`: 上記を満たさない。

PNGは生成有無、画像サイズ、SHA-256を比較する。PNGのSHA-256不一致だけでは総合再現失敗にしない。

## 出力保護

`run_reproduction.py` は `reproduced/` 内に既存ファイルが1個でもあれば停止する。これにより再実行時の上書きを防止する。既存原本成果物へ向く出力変数は、計算開始前に `reproduced/` 配下へ差し替える。

## 最終成果物

最終報告書は次に生成される。

```text
reports/paper7_N5_reproduction_report.md
```

総合判定は `REPRODUCED_EXACTLY`、`REPRODUCED_NUMERICALLY`、`NOT_REPRODUCED` のいずれかだけである。
