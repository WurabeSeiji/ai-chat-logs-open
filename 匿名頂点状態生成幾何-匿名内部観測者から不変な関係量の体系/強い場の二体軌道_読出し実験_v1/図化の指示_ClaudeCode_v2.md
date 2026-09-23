# 図化の指示 v2（Claude Code 用）― 第三思考実験の論文の図

このファイルは、同じフォルダの「図化の指示_ClaudeCode.md」（読み出しの図だけを扱った旧版）を置き換える。

## 目的

論文 `thought_experiment_03_two_body_orbit_readout_ja_v1.md` が参照する図を、すべてそろえる。

## 場所

親フォルダ：Google Drive「匿名頂点状態生成幾何-匿名内部観測者から不変な関係量の体系」（folderId: 1auu7DnHjkrTkTfGqK-8R3kiwMyeZP9_p）

Windows で同期しているときの例：
`G:\マイドライブ\匿名頂点状態生成幾何-匿名内部観測者から不変な関係量の体系\`
ドライブ文字や「マイドライブ」の名前は環境によって違う。見つからなければフォルダ名で探す。

以下のパスは、すべてこの親フォルダからの相対パスである。

## 論文が参照する図と、その状態

| 図 | ファイル（親フォルダからの相対パス） | 状態 |
|---|---|---|
| 図1a | `円錐曲線軌道_形状点列データ生成_v1/figures/overview_input_ja.png` | 未作成（作業1） |
| 図1b | `円錐曲線軌道_形状点列データ生成_v1/figures/overview_truth_ja.png` | 未作成（作業1） |
| 図2 | `強い場の二体軌道_観測1観測2データ生成_v1/figures/obs_D1_ja.png`（obs_D2_ja.png、obs_D3_ja.png も） | 作成済み |
| 図3 | `強い場の二体軌道_観測1観測2データ生成_v1/figures/truth_D1_ja.png`（truth_D2_ja.png、truth_D3_ja.png も） | 作成済み |
| 図4 | `強い場の二体軌道_読出し実験_v1/figures/readout_summary_ja.png` | 未作成（作業2） |
| 図5 | `強い場の二体軌道_読出し実験_v1/figures/readout_details_ja.png` | 未作成（作業2） |

論文の中では、論文の置き場所（`強い場の二体軌道_読出し実験_v1/`）からの相対パスで図を差し込んである。

## 準備

numpy と matplotlib が無ければ入れる。

```
python -m pip install numpy matplotlib
```

## 作業1：数値実験01のデータと図を作り直す

場所：`円錐曲線軌道_形状点列データ生成_v1/`

```
python generate_conic_orbit_samples.py
python reproduce_from_meta.py
python plot_conic_orbit_samples.py
```

- 1本目で `data/input/`（orbit_D01〜D16.csv の16本）と `data/truth/`（正解、meta、manifest.json、validation_report.json）ができる。
- 2本目で `data/truth/reproduction_check.json` ができる。
- 3本目で `figures/overview_input_{ja,en}.{png,svg}` と `figures/overview_truth_{ja,en}.{png,svg}`（計8本）ができる。

乱数のシードは固定（20260922）なので、同じデータができる。プログラムは変えない。

確かめること：

- `validation_report.json` の自己検証が、すべて合格していること。
- `reproduction_check.json` で、再生成が一致していること。
- 論文の §3.1〜3.2 の数値と合っていること。論文には、ケース8通り、16データセット、各 318〜1274 点、最大残差 5×10⁻¹⁵ p 未満、再生成の相対誤差 3×10⁻¹⁵ と書いてある。合わない数値があれば、論文は直さずに報告する。

## 作業2：読み出しの図を作る

場所：`強い場の二体軌道_読出し実験_v1/`

```
python plot_readout_results.py
```

- `figures/readout_summary_{ja,en}.{svg,png}` と `figures/readout_details_{ja,en}.{svg,png}`（計8本）ができる。
- 入力は `results/summary.csv`、`results/features.csv`、`results/reconstruction.csv` の3つだけ。
- 日本語が文字化け（□）していないことを確かめる。「日本語フォントが見つからない」と出て英語版しかできなかったら、スクリプトの `JA_FONT_CANDIDATES` に、その環境にある日本語フォント名を加えて作り直す（Windows なら Yu Gothic か Meiryo）。

図の中身（確認用）：

- `readout_summary_*`：読み出した量の相対誤差の一覧。横軸は17項目で、D1〜D3 の3色の点、読み出しが見積もった 1σ（横棒）、1回の観測の分解能 1/1000 の破線を描く。白抜きの印は、D2（等質量）で正解が 0 の差を τ_A または d_A で割った誤差。
- `readout_details_*`：3行 × 4列。行は D1（a=120, b=112, q=4）、D2（a=30, b=27, q=1）、D3（a=60, b=50, q=2）。列は (a) 近点の方向の進み、(b) 動径周期の減少、(c) 相対距離 r(t) の相対誤差、(d) 固有時間の差 τ_A − τ_B。
- D2 の (d) で、重力波の読み出しで再構成した τ_A − τ_B が約 36 M まで増えるのは正しい結果で、不具合ではない（論文 §6.3）。

## 作業3：論文で図が表示されるかを確かめる

論文 `強い場の二体軌道_読出し実験_v1/thought_experiment_03_two_body_orbit_readout_ja_v1.md` を、Markdown を表示できるエディタ（VS Code のプレビューなど）で開く。図1a〜図5 がすべて表示されることを確かめる。

## してはいけないこと

- `results/` と、数値実験02の `data/` のファイルを書き換えない。
- `readout_binary_pn.py`（読み出し本体）と `generate_binary_pn_observations.py`（数値実験02の生成）は、頼まれない限り実行しない。実行すると、保存してある結果やデータが上書きされる。
- 生成プログラムと図化プログラムは、エラーやフォントの問題がない限り変えない。変えた場合は、どこをなぜ変えたかを報告する。
- この作業では、論文の本文を変えない。

## 報告すること

- 作った図のファイル名とサイズ（作業1と作業2）
- 数値実験01の自己検証と、再生成の照合の結果
- 日本語の表示に問題がないか
- 論文で図がすべて表示されたか
