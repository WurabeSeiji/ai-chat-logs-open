# 図化の指示（Claude Code 用）— 数値実験02 読出しの図

## 目的

読出し実験の結果（`results/` の CSV）から図を作り、このフォルダの `figures/` に保存する。
作る図は SVG と PNG の日本語版・英語版、計8本。Google Drive のデスクトップ版で同期していれば、`figures/` に書き出すだけで Drive に保存される。

## 場所

Google Drive:
`匿名頂点状態生成幾何-匿名内部観測者から不変な関係量の体系/強い場の二体軌道_読出し実験_v1/`
（folderId: 1N6RnwSGlbjgRDXdo3FvhhCOd7PELmMBH）

Windows で同期しているときの例:
`G:\マイドライブ\匿名頂点状態生成幾何-匿名内部観測者から不変な関係量の体系\強い場の二体軌道_読出し実験_v1\`
ドライブ文字や「マイドライブ」の名前は環境によって違う。見つからなければフォルダ名で探す。

使うファイル（すべてこのフォルダにある）:
- `plot_readout_results.py`（図化プログラム）
- `results/summary.csv`、`results/features.csv`、`results/reconstruction.csv`（入力はこの3つだけ）

## 手順

1. このフォルダに移動する。
2. numpy と matplotlib が無ければ入れる: `python -m pip install numpy matplotlib`
3. 図を作る: `python plot_readout_results.py`
4. `figures/` に次の8本ができたことを確かめる。
   - `readout_summary_ja.svg` / `.png`、`readout_summary_en.svg` / `.png`
   - `readout_details_ja.svg` / `.png`、`readout_details_en.svg` / `.png`
5. `readout_summary_ja.png` と `readout_details_ja.png` を開き、日本語が文字化け（□）していないことを確かめる。
   「日本語フォントが見つからない」と表示されて英語版しかできなかった場合は、スクリプトの `JA_FONT_CANDIDATES` に、その環境にある日本語フォント名を加えて作り直す。Windows なら Yu Gothic か Meiryo があるはず。
6. Drive の同期が終わったら、ブラウザで Drive の `figures/` に8本あることを確かめる。

同期していない場合は、上の4ファイルを同じフォルダ構成（`results/` の下に CSV）でダウンロードして手順3〜5を行い、できた `figures/` フォルダを Drive のこのフォルダにドラッグする。

## 図の中身（確認用）

- `readout_summary_*`: 読み出した量の相対誤差（正解との差）の一覧。
  - 横軸は17項目。q、ν、a、b、T_M、D、τ_A − τ_B、d_A − d_B を、観測1と重力波（観測2）の方法ごとに並べる。
  - データセット D1〜D3 の3色の点、読出しが見積もった 1σ（横棒）、1回の観測の分解能 1/1000 の破線を描く。
  - 白抜きの印は、D2（等質量）で正解が 0 の差を τ_A または d_A で割った誤差。
- `readout_details_*`: 3行 × 4列。行は D1（a=120, b=112, q=4）、D2（a=30, b=27, q=1）、D3（a=60, b=50, q=2）。
  - (a) 近点の方向の進み
  - (b) 動径周期の減少
  - (c) 相対距離 r(t) の相対誤差
  - (d) 固有時間の差 τ_A − τ_B
- D2 の (d) で、重力波の読出しで再構成した τ_A − τ_B が約 36 M まで増えるのは正しい結果で、不具合ではない。重力波から読んだ ν = 0.24979 が 0.25 よりわずかに小さいので、質量差 0.03 として再構成されるため（README の「注意」を参照）。

参考: こちらの環境（Noto Sans CJK）で作ったときのサイズ（bytes）。フォントなどの違いで変わるので、目安にとどめる。
- PNG: summary_ja 114421、summary_en 99490、details_ja 558373、details_en 536753
- SVG: summary_ja 42853、summary_en 48611、details_ja 260534、details_en 293203

## してはいけないこと

- `results/` の CSV と JSON を書き換えない。
- `readout_binary_pn.py`（読出し本体）は、頼まれない限り実行しない。実行すると約11分かかり、`results/` が上書きされる。
- 図化プログラムは、エラーやフォントの問題がない限り変えない。変えた場合は、どこをなぜ変えたかを報告する。
- `upload_to_drive.py` はこの手順では使わない（同期で保存されるため）。

## 報告すること

- 8本のファイル名とサイズ
- 日本語の表示に問題がないか
- Drive の `figures/` に保存されたか
