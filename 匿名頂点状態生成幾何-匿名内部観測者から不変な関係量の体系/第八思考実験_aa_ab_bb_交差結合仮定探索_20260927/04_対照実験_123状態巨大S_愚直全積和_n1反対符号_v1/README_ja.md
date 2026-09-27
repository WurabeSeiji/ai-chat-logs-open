# 第八思考実験・対照実験 v1 — 123状態・巨大S・愚直全積和版

## 実験条件

今回残す条件は1つだけ。

- `n=1`
- `aa: (C,D)=(0.91,0)`
- `ab: (C,D)=(1.09,0.36)`
- `bb: (C,D)=(0.91,0)`

非対角41×41ブロックは全て構造的に0とする対照実験であり、交差物理は追加しない。

## 実装上の固定条件

- 状態は最初から最後まで `psi[123]` 1本。
- `aa/ab/bb` の41状態部分ベクトルを作らない。
- `reshape`, `concatenate`, 41状態 `transition()` を生成経路で使わない。
- 相互作用配列は `S[123,123]` 1枚を最初に確保し、最後まで同じ配列オブジェクトを使う。
- 41×41小行列を作らない。block copyしない。
- 毎microstep、`S` の全15,129セルを再評価して代入する。
- `q=0`、`psi[j]=0`、`S[i,j]=0` を理由に生成・積和をスキップしない。
- one-hot 11相は `argmax(q)` で1相だけ選ばず、11候補全てを `S` のq列として同時に構成する。
- 積和は全 `i=0..122`, `j=0..122` を訪問する二重ループで行う。
- 1 macrostep = 11 microsteps。
- microstep状態は保存しない。
- 初期状態と、完了した全macrostepの `psi[123]` は1行も間引かずHDF5へ保存する。
- 実験本体はPaper 6/7の参照rawを一切読まない。

## 停止条件

3チャネル全てについて `P <= 20` になったmacro境界で終了する。
早く `P<=20` に到達したチャネルも、巨大な1系の一部なので凍結せず同じ写像で更新し続ける。

## raw保存

各rowは

- `macro_step`: int64
- `psi123[123]`: float64

のみを正本として保存する。
初期状態を `macro_step=0` として保存し、その後全macrostepを保存する。

100,000 rowごとに `raw_macro_partNNN.h5` へpart分割する。part分割はディスク保存の都合だけであり、計算状態 `psi[123]` は分割しない。

## バックグラウンド進捗

本体は以下を行う。

1. `print(..., flush=True)` により標準出力へ定期表示。
2. `progress.json.tmp` を書いて `os.replace()` することで `progress.json` を原子的更新。
3. 30秒経過時にはmacro途中でも `micro_in_current_macro` を更新。
4. 100 macrostepごとにHDF5をflushし、macro境界進捗を表示。

`progress.json` には少なくとも以下を記録する。

- status
- macro_completed
- micro_in_current_macro
- aa/ab/bb のP
- C,D
- q sum / nonzero count
- raw_rows_saved
- current_raw_part
- elapsed_seconds
- macrosteps_per_second
- physical_progress_fraction
- logical_product_terms_visited

進捗表示専用の `show_progress_full123_v1.py` は `progress.json` を読むだけで、実験状態へ書き込まない。

## ファイル

- `run_paper8_control_full123_bruteforce_n1_opposite_v1.py` — 実験本体
- `show_progress_full123_v1.py` — 読み取り専用進捗表示
- `README_ja.md` — この監査仕様

**この段階では実行しない。コード監査後にのみ実行する。**
