# Zenodo deposit（四元数の関係読み出しと、結晶化としてのインフレーション v1.1）

- 作成日：2026-10-21（API、新規 deposit。HTTP 201）
- Deposit ID：23523238
- Version DOI（prereserve）：10.5281/zenodo.23523238
- Concept DOI：10.5281/zenodo.23523237（conceptrecid 23523237）
- 状態：**published 2026-10-21**（files 11、version v1.1）
- bucket：https://zenodo.org/api/files/98029a9a-92d5-4df3-92dd-03a122579093
- 指示：「コミットプッシュはあとにして、さきに、公開プロセスを進めます。DOI コンセプトDOIを取得、日本語論文に記入、英訳、tex pdf 作成　アップロード、zennアーティクル作成、リリースノート更新、PAPERS更新　まで一連の公開プロセスを実行してください。」（2026-10-21）
- 版番号：v1 は公開前稿（Zenodo に出していない）。書誌の修正を入れた v1.1 を Zenodo の初版として公開した。

## メタデータ

- upload_type publication / publication_type article / language jpn / license cc-by-4.0 / version v1.1
- creators：Kihara, Noriaki（WF System Co., Ltd.、ORCID 0009-0004-6753-4020）
- 英題：Relational Readout of Quaternions and Crystallization as Inflation: Anonymity-Preserving Autonomous Structuring and Partial Projection
- keywords：quaternions, anonymity, partial projection, stereographic projection, Delaunay triangulation, rigidity, geometrical frustration, 600-cell, supercooling, inflation, split octonions
- related_identifiers（すべて references、Concept DOI）：23451494（複素数の集合から時空への写像）、23444451（目録）、23384928（統合仮説論文）、22729146（自己無撞着な関係波閉鎖系の力学）、21315735（基本公理系 v9）、19427779（中心投影による4次元空間の幾何学的定式化）
- 作成に使ったメタデータ JSON と応答はスクラッチ（`zenodo_qci/metadata_create.json`、`zenodo_qci/create_resp.json`）

## 進捗

1. DOI 取得（済）→ 日本語 md v1.1 のヘッダーと版の履歴に Version DOI・Concept DOI を記入（済）
2. 英訳 `quaternion_crystallization_inflation_en_v1.1.md`（済。日本語版が正本で英語版は翻訳と明記。§2.2 の割り算・原点・観測者の宣言の表と、§13 の主張の区分の表はそのまま訳した。図のリンクは `図表/en/*_en.svg`）
3. 日英 tex/PDF（済）。`build_tex.py` は前稿（四元数による時空への写像）のものを元に、(a) 図の直後の説明文の段落（**図 N**／**Figure N**）を図の説明文にまとめて `figs/*.png` を読む置換、(b) 冒頭の著者・日付・DOI の行を行ごとに改行、(c) 図の PNG を `図表/ja`・`図表/en` からコピー、を加えた。/tmp/tex_compile/qci で pandoc → lualatex ×2、日本語 20 頁・英語 21 頁、LaTeX エラー 0・欠字 0。英語版の PDF に残る日本語はフォルダ名とファイル名だけ
4. bucket へ PUT（済、11 ファイル、curl -4、サーバー側 MD5 は手元と一致）：日英 md・tex・PDF、build_tex.py、figures_and_programs.zip（`図表/`・図版生成プログラム 2 本・フォルダ直下の図 9 枚）、書誌照合の記録、追記（対話での気づき）、SHA256SUMS.txt
5. publish（済、2026-10-21、HTTP 202、state done。無認証 GET：doi 10.5281/zenodo.23523238、conceptdoi 10.5281/zenodo.23523237、publication_date 2026-10-21、files 11。PDF 日英の download URL は 200）
6. Zenn 記事 `articles/quaternion-crystallization-inflation.md`、RELEASE_NOTES.md、PAPERS_ja/en.md の 53 番・docs/index.html・docs/papers.json・docs/p/23523238.html の追記（済、全再生成はしていない）
7. コミット・プッシュ：木原指示により後で行う（未実施）

## 未着手（引き継ぎメモ §2-3）

- 引用の対応表 2 本（`自己引用の対応表_20261021.md`、`外部引用の対応表_20261021.md`）は作っていない。
