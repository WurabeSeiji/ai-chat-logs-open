# Zenodo deposit（対称性とマジックナンバーの目録 v1.0）

- 作成日：2026-10-18（API、新規 deposit。初版なので newversion ではない。HTTP 201）
- Deposit ID：23444452
- Version DOI（prereserve）：10.5281/zenodo.23444452
- Concept DOI：10.5281/zenodo.23444451（conceptrecid 23444451）
- 状態：**published 2026-10-18**（publish は HTTP 202、state done・submitted True。無認証 GET：doi 10.5281/zenodo.23444452、conceptdoi 10.5281/zenodo.23444451、version v1.0、publication_date 2026-10-18、files 9、status published。PDF の download URL は日英とも 200）
- bucket：https://zenodo.org/api/files/9565a752-186f-495f-a532-7e124648c6e2
- 指示：「公開ルールにしたがい、DOI コンセプトDOIを取得、英訳、tex/pdf作成、zennアーティクル作成、アップロード、リリースノート更新、PAPERS更新、コミットプッシュまでをおこなってください。」（2026-10-18）

## メタデータ

- upload_type publication / publication_type article / language jpn / license cc-by-4.0 / version v1.0
- creators：Kihara, Noriaki（WF System Co., Ltd.、ORCID 0009-0004-6753-4020）
- 英題：A Catalogue of Symmetries and Magic Numbers for Reading Out Spacetime from a State Space: The Answer-Check Standard of a Research Programme Taking Anonymity and Self-Consistency as Axioms
- related_identifiers（すべて references、Concept DOI）：21315735（基本公理系 v9）、21521899（無名部分投射と E₈）、22028072（零閉包・有限位数・自己無撞着幾何）、22112008（自己無撞着インフレーション機構）、22851944（研究プロジェクト設計）、23315557（第十四）、23384928（統合仮説論文）
- 作成に使ったメタデータ JSON と応答はスクラッチ（`zenodo_smc/metadata_create.json`、`zenodo_smc/create_resp.json`）

## 進捗

1. DOI 取得（済）→ 日本語 md に Version DOI・Concept DOI・日付・版管理を記入（済）
2. 自己引用の主張のうち [S12]（自己無撞着インフレーション機構 v2 §4「ゼロ閉塞は定理である」）と [S24]（N＝1, K＝4 の最小厳密閉系）を原文で確認（済）
3. 英訳（済。日英で h2 22・h3 32・表の行 169・参考文献 S32・X33・表示式 10 が一致）。付録 A に一度書いた架空のローマ字パスを削除し、実際のパスだけにした
4. 日英 tex/PDF（済。build_tex.py は統合仮説論文 v1.1 のコピーで、prep の引数（ファイル名・題名・副題・日付と DOI）を変更し、tex 化の前処理に ℓ → $\ell$ を一つ追加（[X16] の Cℓ(8) が DejaVu Serif で欠字になるため）。/tmp/tex_compile/smc で pandoc → lualatex ×2、ja 26 頁・en 28 頁、LaTeX エラー 0・欠字 0・はみ出し 0。英語版の undefined は和文斜体の字形の情報のみ。pdftotext で表の区切りの崩れがないこと、表の主要な行と DOI を確認）
5. bucket へ PUT（済、9 ファイル、すべて HTTP 201、サーバー側 MD5 は手元と一致）：md・tex・pdf 各 2、build_tex.py、自己引用の対応表_20261017.md、SHA256SUMS.txt（8 ファイル分）。最初の試行はシェルのファイル列の扱いの誤りで失敗（空の SHA256SUMS.txt が 400 で拒否され、bucket は空のまま）、Python でやり直した
6. publish（済、2026-10-18）
7. Zenn 記事、RELEASE_NOTES.md、PAPERS_ja/en.md・docs/index.html・docs/papers.json・docs/p/23444452.html の追記、コミット・プッシュ（同日）
