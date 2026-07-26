# Stage A1: 第7論文 N=5 イベント定義校正

## 目的

Stage A0で完全再現されたN=5 CSVだけを読み取り、人間がイベント定義を検討するための候補表と可視化を作る。単一の指数成長開始・終了時刻、`rank_Q=4` 持続開始時刻は採用しない。

このパッケージは後処理専用であり、原本力学コードのimport・実行、新しい軌道生成、第8論文本実験、高精度計算、Series 1〜3、N=40、N=300を行わない。

## 入力

読み取るファイルは、兄弟ディレクトリ `paper7_N5_reproduction/reproduced/` 内の次の3件だけである。

- `fcurve_N00005_delta1e-15_seed0.csv`
- `q_svd_N00005.csv`
- `paper7_long_timeseries.csv`

`verify_inputs.py` は、Stage A0報告書の総合判定とCSV別 `exact` 判定、Stage A0比較記録内の再現SHA-256、実ファイルのSHA-256を照合する。

## 時間軸

- `fcurve.tau` は絶対stepであり、0から21167まで1 step刻みである。
- `q_svd.step` と `q_svd.time` は等しく、全保存行で `step - relative_time = 1167` である。
- `paper7_long_timeseries.step` と `time` は等しく、25 step間隔である。
- 3入力は同じ絶対stepを使うが、保存間隔と終端は異なる。未保存stepのq値は補間しない。

したがって、rank候補の持続長 `1, 5, 10, 20, 40, 80, 160` は「連続するq保存レコード数」として計算し、開始・確認stepと実際の観測spanも併記する。これは最終イベント定義ではなく、人間判断のための候補規約である。

## f(t)の後処理

`f <= 0` は変更せず、対数とその差分を `NaN` にする。

- `log_f = ln(f)`
- `log10_f = log10(f)`
- `log_f_diff_1(t) = log_f(t) - log_f(t-1)`
- `log_f_central_diff(t) = [log_f(t+1)-log_f(t-1)]/[step(t+1)-step(t-1)]`

回帰は各奇数窓の中心を代表時刻とし、`log_f = intercept + slope × step` の通常最小二乗を使う。`R² = 1-SSE/SST`、回帰残差標準偏差は `sqrt(SSE/(window-2))` とする。完全窓を持たない端点、非有限値を含む窓、`SST=0` は該当量を `NaN` にする。

## 候補区間

各 `window × R²閾値 × 最小持続長` について、

```text
slope > 0 かつ R² >= 閾値
```

を満たす最大連続区間をすべて保存する。`duration` は1 step刻みのf系列上の保存点数で、区間両端を含む。同じ区間が複数の候補パラメータを満たす場合も削除しない。

成長終了候補は各区間の次step以後を同じ回帰窓のslopeで探索する。条件A/B/Cと持続長の全組合せを保存し、見つからない組合せも `not_found` 行として残す。

## rank候補

既存列 `rank_q` が、厳密不等号

```text
count(q_j > 1e-8 × q1), j=1,...,4
```

と全行で一致することを検査する。比較用閾値 `1e-6`〜`1e-12` でも同じ後処理を行うが、既存定義を置き換えない。

q比はゼロ除算時に `NaN` とし、固定閾値によるイベント採用には使わない。

## 実行順序

```bash
PYTHONDONTWRITEBYTECODE=1 python3 paper7_N5_event_calibration/verify_inputs.py
PYTHONDONTWRITEBYTECODE=1 python3 paper7_N5_event_calibration/analyze_growth_candidates.py
PYTHONDONTWRITEBYTECODE=1 python3 paper7_N5_event_calibration/analyze_rank_candidates.py
PYTHONDONTWRITEBYTECODE=1 python3 paper7_N5_event_calibration/make_calibration_figures.py
PYTHONDONTWRITEBYTECODE=1 python3 paper7_N5_event_calibration/make_calibration_report.py
```

実際にはリポジトリルートから各スクリプトの完全な相対パスを指定する。各スクリプトは直前工程の成功記録を検査し、自身の出力が既にあれば上書きせず停止する。

## 最終状態

報告書は `reports/paper7_N5_event_calibration_report.md` に生成する。総合状態は次のいずれかだけである。

- `CALIBRATION_DATA_COMPLETE`
- `CALIBRATION_DATA_INCOMPLETE`
- `INPUT_MISMATCH`

Stage A1完了後は停止し、Stage A2へ進まない。
