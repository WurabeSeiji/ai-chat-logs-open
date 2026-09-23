# 公開指図（Claude Code 用）― 第三思考実験の公開準備

## 0. 前提

- 論文の正本は `強い場の二体軌道_読出し実験_v1/thought_experiment_03_two_body_orbit_readout_ja_v1.5.md`（版 v1.5、公開前最終版）。同じフォルダにある v1・v1.1〜v1.4 は旧版で、使わない。
- 論文の本文は、この指図に書いた機械的な変更（図のパス、所在の節、DOI）以外は変えない。数値実験を再実行しない。`results/` と `data/` を変えない。
- 用語の統一：距離量は「時計率換算距離」。旧名称「レーダー距離」は README とプログラムのコメントから消す。ただし CSV・JSON のキー名（`radar d_A mean` など）は互換性のため残し、README に旧名称であると一行書く。
- 分からないことがあれば、推測で進めずに止まって報告する。

作業フォルダ（Google Drive をデスクトップで同期しているパスの例。ドライブ文字や「マイドライブ」の名前は環境で違う）：
- 数値実験01：`匿名頂点状態生成幾何-匿名内部観測者から不変な関係量の体系/円錐曲線軌道_形状点列データ生成_v1/`
- 数値実験02：`同/強い場の二体軌道_観測1観測2データ生成_v1/`
- 読み出し：`同/強い場の二体軌道_読出し実験_v1/`
- 公開先（GitHub の同期フォルダ）：`マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/`
- 先行論文（体裁と用語の手本）：`同/基礎物理モデル研究プロジェクト設計_20260920/`（論文0・1）、`同/第二思考実験_クーロン力の積読み出し_20260921/v1_2_実数配列版_20260921/`（論文2 v1.2）

`次元の生成構造/第三思考実験_三状態関係力学_20260921/` は別の未完の試みで、本稿とは無関係。触らない。

## 1. 公開フォルダを作る

`次元の生成構造/第三思考実験_二体軌道の読出し_20260923/` を作る。以下のファイルはすべてこの中に置く。

## 2. 論文本体（日本語）

1. v1.5 を `thought_experiment_03_two_body_orbit_readout_ja_v1.md` の名前でコピーする（版の表記は v1.5 のまま。論文2が「ファイル名は _v1、版は v1.2」としているのと同じ扱い）。
2. 図のパスを、§3 で公開フォルダに置く名前に書き換える。§4.4 と §4.5 の「D2、D3 は同じフォルダの obs_D2_ja.png、obs_D3_ja.png に示す」などの文も、新しいファイル名に直す。
3. 「データとプログラムの所在」の節を、Google Drive のフォルダ名ではなく、公開フォルダ内のファイル名（§4・§5 の一覧）で書き直す。
4. DOI の欄は「（公開時に付与）」のまま残す（§9 で入れる）。

## 3. 図

各実験の `figures/` から公開フォルダにコピーし、次の名前にする。PNG と SVG の両方、日本語版と英語版の両方。

| 論文の図 | 元のファイル（{lang} = ja, en） | 公開名 |
|---|---|---|
| 図1a | 円錐曲線…/figures/overview_input_{lang}.png/svg | fig01a_conic_shape_input_{lang}_v1 |
| 図1b | overview_truth_{lang} | fig01b_conic_shape_truth_{lang}_v1 |
| 図2 | 強い場…データ生成/figures/obs_D1_{lang}（D2、D3 も） | fig02_pn_observations_D1_{lang}_v1（D2、D3 も） |
| 図3 | truth_D1_{lang}（D2、D3 も） | fig03_pn_truth_D1_{lang}_v1（D2、D3 も） |
| 図4 | 読出し実験/figures/readout_summary_{lang} | fig04_readout_summary_{lang}_v1 |
| 図5 | readout_details_{lang} | fig05_readout_details_{lang}_v1 |

数値実験01の figures/ が無ければ、`円錐曲線軌道_形状点列データ生成_v1/` で `generate_conic_orbit_samples.py`、`reproduce_from_meta.py`、`plot_conic_orbit_samples.py` を順に実行して作る（乱数のシードは固定なので同じデータができる）。

## 4. README とプログラム

公開フォルダにコピーし、用語を直す。

- `readout_binary_pn.py`、`plot_readout_results.py`、読み出しの `README.md`：コメント・docstring・本文の「レーダー距離」を「時計率換算距離」に直す。プログラムの動作は変えない（コメントと文字列だけ）。直したあと `python -c "import ast; ast.parse(open('readout_binary_pn.py', encoding='utf-8').read())"` で構文を確かめる。
- 読み出しの `README.md` に一行加える：「results の CSV・JSON のキー `radar d_A mean` などは旧名称で、本文の時計率換算距離 d_A に当たる」。
- 読み出しの README の結果表と注意は、論文 v1.5 の §6・§7.4 と矛盾しないように直す（χ²/自由度を適合度として書いている箇所があれば、論文 §6.3 と同じ説明にする）。
- `generate_binary_pn_observations.py`、`plot_binary_pn_observations.py`、数値実験02の `README.md`、`generate_conic_orbit_samples.py`、`reproduce_from_meta.py`、`plot_conic_orbit_samples.py`、数値実験01の `README.md` はそのままコピー（用語が出てくれば同じ修正）。

## 5. 再現用ファイル

- `results/`：summary.csv、features.csv、reconstruction.csv、readout_results.json、stepsize_check_D2.csv
- `data/`：数値実験02の data/input（obs1_D1〜3.csv、obs2_D1〜3.csv）と data/truth（D1〜D3_meta.json、manifest.json、D1〜D3_timeseries.csv）。数値実験01の data/input と data/truth。大きければ `reproduction_data_v1.zip` に一つにまとめてもよい（中身の一覧を README に書く）。
- 公開フォルダの全ファイルの `SHA256SUMS_v1.txt` を作る（論文2の v1_2 フォルダと同じ形式）。DOI を入れて再ビルドしたあとに作り直す。

## 6. 英語版

`thought_experiment_03_two_body_orbit_readout_en_v1.md` を作る。

- 日本語版 v1.5 の忠実な翻訳。節の構成、式、表、数値、参考文献はそのまま。意味を足したり削ったりしない。
- 用語と札の英語は、論文1・2の英語版（`thought_experiment_01_…_en_v1.md`、`thought_experiment_02_…_en_v1.md`）の訳語に合わせる（読み出し、設計仮定、【原則】【設計仮定】【数値】【照合】、無名、ゼロ閉塞 など）。合わせた訳語の対応表を報告に付ける。
- 題名の案：*Third Thought Experiment: Reading Out the State of Two Bodies a, b from the Orbit of Two Celestial Bodies a, b — Mass ratio, relative distance, and relative time from orbit and gravitational-wave readouts, with two values of undetermined meaning and the knowledge G = c = M = 1 alone*。論文1・2の英語題名の体裁に合わせて整えてよい。
- 新しく作った用語の英訳：時計率換算距離 = clock-rate-converted distance、試験観測者 = test observer、視野角 = viewing angle、三体問題への入口 = the entry point to the three-body problem。
- 図は英語版のファイル（fig…_en_v1）を参照する。

## 7. TeX と PDF

- 日本語版・英語版それぞれの `.tex` と `.pdf` を作る（`…_ja_v1.tex/.pdf`、`…_en_v1.tex/.pdf`）。
- 体裁と変換手順は論文2の `.tex` に合わせる（プリアンブル、フォント、図の入れ方）。手順が分からなければ、論文2の tex を読んで再現できる範囲で行い、できない点を報告する。
- PDF の図は PNG を使う。数式・表・参考文献が md と一致することを、PDF を開いて確かめる。

## 8. note 記事と告知文

- `note_thought_experiment_03_two_body_orbit_readout_ja.md` と `…_en.md`：一般向け。数式と表は使わない。ハッシュタグを付ける。画像は埋め込めないので、挿入位置とファイル名を本文中にはっきり書く。曖昧な言い訳や但し書きは入れず、起きた事実を簡潔に書く。内容は論文 v1.5 の範囲を越えない（素粒子スケールへの展開などは書かない）。
- `fb_message_…_ja/en.md`、`x_message_…_ja/en.md`：論文2の同名ファイルと同じ体裁で短文を用意する。DOI は §9 のあとで入れる。

## 9. DOI と最終ビルド

1. ここまでを報告し、木原さんが Zenodo に登録して Version DOI と Concept DOI を受け取る。
2. 受け取った DOI を、日本語版・英語版の md の冒頭（Version DOI、Concept DOI）に入れる。日付を公開日に直す。
3. tex と pdf を作り直す。note 記事と告知文に DOI を入れる。
4. `SHA256SUMS_v1.txt` を作り直す。
5. `git add` と `git commit`（メッセージ例：`Add thought experiment 03: two-body orbit readout (v1.5)`）まで行い、push は木原さんの確認を得てから行う。

## 10. 報告すること

- 公開フォルダのファイル一覧とサイズ
- 論文 md で書き換えた箇所（図のパス、所在の節）の一覧
- 用語修正を行ったファイルと箇所
- 英語版の訳語の対応表と、判断に迷った箇所
- TeX/PDF の変換で使ったコマンドと、できなかったこと
- 未実行のもの（DOI 待ちの項目）
