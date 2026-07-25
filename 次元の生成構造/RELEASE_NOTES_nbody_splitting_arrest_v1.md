# Release Notes: 自発的分裂の停止と新しい直交回転平面の創発（次元の生成構造 第6論文）

Release date: 2026-07-25

## DOI

- Concept DOI: https://doi.org/10.5281/zenodo.21543070
- Version DOI (v1): https://doi.org/10.5281/zenodo.21543071
- Zenodo record: https://zenodo.org/record/21543071

（公開直後は Version DOI の doi.org 解決に DataCite 伝播の遅延がある。Zenodo レコード API は公開直後から open で解決する。）

## Zenn

- https://zenn.dev/noriaki_kihara/articles/nbody-splitting-arrest-orthogonal-plane-emergence

## note（一般向け）

- 日本語: https://note.com/kiharanoriaki/n/n28eeb5465c32
  （元原稿: `次元の生成構造/波はなぜ分裂を止めるのか_note完成稿.md`）
- 英語: https://note.com/kiharanoriaki/n/ne2967054b888
  （元原稿: `次元の生成構造/why_does_the_wave_stop_splitting_note_en.md`）

## Facebook

- 日本語投稿: https://www.facebook.com/kihara.noriaki/posts/pfbid027R1e2p5G8xi24iEmm7159EYn78H6mjDQQ9EHEns53p8jRCnwoT9quk6D2dRgP426l

## 公開ファイル（Zenodo, 7点）

- `nbody_splitting_arrest_and_orthogonal_plane_emergence_ja.md` / `_en.md`
- `nbody_splitting_arrest_and_orthogonal_plane_emergence_ja.tex` / `_en.tex`
- `nbody_splitting_arrest_and_orthogonal_plane_emergence_ja.pdf` / `_en.pdf`
- `nbody_splitting_arrest_reproduction_v1.zip`（原本2本＋計装7本＋図4枚＋結果md2件）

## 主張

単一波の自発的分裂は、新しい直交回転方向の創生とともに停止する。分裂量 f は支配回転平面から流出した保存ノルムの割合に等しく（恒等式 f = 1 − E_P1、偏差 ≤3.7e-14）、流出先はその他の直交回転部分空間と生成子核である。第2非零回転平面は N≤4 で不在・N≥5 で出現し、前論文の持続的分裂境界と一致する。準安定振幅は A_rel ≈ 1/√M（M=N(N−1)/2、38試行、べき指数 1.00820）。

## 方法の透明性

平面流入は厳密法（密行列 eig(K)、閾値非依存、N=5,40）と近似法（低ランク JG、明示閾値 σ_rel=1e-6、N=300）に分離。N=40 で両法が機械精度 2e-15 で一致することを確認。近似の閾値は CLI・出力JSON・図タイトルに明記。N≥300 の厳密計算は大メモリを要する今後の課題。

## 引用

- 自己引用: 生成子ランク三方向飽和（Concept DOI 21465898）、平面分解読出し（21468959）、自発的分裂の開始と帰結分類（21486233）、波の数は分解能（21486544）
- 二公理の定義: 基本公理系 v9.1（Concept DOI 21315735）
- 外部: Horn-Johnson（歪対称正準形）、Taghavi-Chabert（純スピノル/ツイスター）、Walker（回転からの時間創発）、Smith / Furey / Todorov / Manogue-Dray-Wilson（D4・八元数・E8 との比較対象）

## 生成手順

- 数式は GitHub 記法（$$ / $）、tex/PDF は /tmp/tex_compile で生成（日本語=lualatex、英語=pdflatex、各2回、LaTeX Error 0）。英語版はテーブルの Unicode（−, σ）を ASCII/数式へ、CJK データ名をローマ字へ前処理。
- CSV は .gitignore（*.csv）により追跡外。
