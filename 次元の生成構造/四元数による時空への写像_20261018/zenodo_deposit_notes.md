# Zenodo deposit（複素数の集合から時空 (1＋3) への写像 v1.0）

- 作成日：2026-10-18（API、新規 deposit。初版なので newversion ではない。HTTP 201）
- Deposit ID：23451495
- Version DOI（prereserve）：10.5281/zenodo.23451495
- Concept DOI：10.5281/zenodo.23451494（conceptrecid 23451494）
- 状態：**published 2026-10-18**（files 15、version v1.0）
- bucket：https://zenodo.org/api/files/fc4510e9-1401-4468-81ce-bdb83243ae4f
- 指示：「では、DOI、コンセプトDOIを取得してて日本語mdに埋め込み、英訳、tex pdf 作成　アップロード、zenn アーティクル作成、リリースノート作成、PAPERS 追記、コミットプッシュを実行してください。」（2026-10-18）

## メタデータ

- upload_type publication / publication_type article / language jpn / license cc-by-4.0 / version v1.0
- creators：Kihara, Noriaki（WF System Co., Ltd.、ORCID 0009-0004-6753-4020）
- 英題：A Map from a Set of Complex Numbers to Spacetime (1+3): Readout by Quaternions and Partial Projection
- related_identifiers（すべて references、Concept DOI）：23444451（対称性とマジックナンバーの目録）、23384928（統合仮説論文）、21902805（ゼロ閉塞は4次元だった）、20665661（双対幾何 論文14）、19427779（中心投影による4次元空間の幾何学的定式化）、21495305（零二乗和拘束とスケール不変性の幾何学的正体）、22112008（自己無撞着インフレーション機構）、23315557（第十四思考実験）
- 作成に使ったメタデータ JSON と応答はスクラッチ（`zenodo_qsm/metadata_create.json`、`zenodo_qsm/create_resp.json`）

## 進捗

1. DOI 取得（済）→ 日本語 md のヘッダーと版管理に Version DOI・Concept DOI を記入（済）
2. 英訳（済。日英で h2 20・h3 26・表の行 46・参考文献 S30・X26・図 6・本文の引用 59 が一致）
3. 日英 tex/PDF（済。build_tex.py は目録 v1.0 のコピーで、prep の引数（ファイル名・題名・副題・日付と DOI）を変更し、前処理に三つ追加：図のパス `](図表/figures/` → `](figs/`、`≅` → `$\cong$`（DejaVu Serif に字形がないため）、図の説明文の先頭の「図 N　」「Figure N. 」を外す（LaTeX が番号を付けるので二重になる）。/tmp/tex_compile/qsm で pandoc → lualatex ×2、日英とも 22 頁、LaTeX エラー 0・欠字 0、はみ出しは raw の URL の 0.2pt のみ）
4. bucket へ PUT（済、15 ファイル、サーバー側 MD5 は手元と一致）：md・tex・pdf 各 2、build_tex.py、自己引用・外部引用の対応表、SHA256SUMS.txt、figures_tables_v1.zip（`図表/` を zip にしたもの。リポジトリには置かず Zenodo のみ）、HTML 4 本（`../四元数を基本量とする思考実験_20261018/`）。最初の試行は Python の urllib が Zenodo への IPv6 接続の確立待ちで止まり、1 本目のあと進まなかった。curl -4 に切り替えて残り 14 本を送った
5. publish（済、2026-10-18、HTTP 202、state done・submitted True。無認証 GET：doi 10.5281/zenodo.23451495、conceptdoi 10.5281/zenodo.23451494、version v1.0、publication_date 2026-10-18、files 15、status published。PDF 日英と HTML の download URL は 200）
6. Zenn 記事 `articles/quaternion-spacetime-map.md`、RELEASE_NOTES.md、PAPERS_ja/en.md の 52 番・docs/index.html・docs/papers.json・docs/p/23451495.html の追記、コミット・プッシュ（同日）
