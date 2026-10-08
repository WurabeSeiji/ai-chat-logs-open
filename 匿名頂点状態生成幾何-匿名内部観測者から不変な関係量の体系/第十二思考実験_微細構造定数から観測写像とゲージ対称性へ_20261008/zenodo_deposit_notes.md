# Zenodo deposit（第十二思考実験 v1.0）

- 作成日：2026-10-08（API、新規 deposit。初版なので newversion ではない）
- Deposit ID：23226259
- Version DOI（prereserve）：10.5281/zenodo.23226259
- Concept DOI：10.5281/zenodo.23226258（conceptrecid 23226258）
- bucket：https://zenodo.org/api/files/c6ee7dea-4044-4e62-9007-7a903a56593f
- 状態：**published 2026-10-08**（HTTP 202。無認証 GET：doi 10.5281/zenodo.23226259、conceptdoi 10.5281/zenodo.23226258、files 15、version v1.0、publication_date 2026-10-08）

## メタデータ（登録済み）

- upload_type publication / publication_type article / language jpn / license cc-by-4.0 / version v1.0
- creators：Kihara, Noriaki（WF System Co., Ltd.、ORCID 0009-0004-6753-4020）
- related_identifiers（すべて references）：10.5281/zenodo.22851944（設計書）、23199299（第十一思考実験 Concept）、21422505（非自明二乗閉鎖 Concept）、19902677（xyztRQ の再検討 Concept）、20060728（CP-Comp Concept）
- 英題：The Twelfth Thought Experiment: From the Fine-Structure Constant to Observation Maps and Gauge Symmetry - Starting from the Fine-Structure Constant, Rereading Spacetime, Mass, Charge, Color Charge and Gauge Symmetry from the Light-Cone Equation and the Choice of Observable Axes

## 進捗

1. 日本語 md に DOI 埋め込み（済）→ 英訳 md（済）→ 日英 tex/PDF（済。`build_tex.py`（第十一のコピー、ファイル名・題・日付・DOI のみ変更）で /tmp/tex_compile/te12 に生成：pandoc → lualatex ×2、ja 24 頁・en 20 頁、LaTeX エラー 0、欠字 0、はみ出し 0。英語版の箱 2 か所を array で 2 行に分けた）
2. bucket へ PUT（済、全 15 ファイル HTTP 201）：md ×2、tex ×2、pdf ×2、fig01〜03 の png・svg、make_figures.py、build_tex.py、citation_ledger_20261008.md（引用文献台帳_20261008.md を改名）
3. publish（済、2026-10-08、木原氏の指示「DOI 取得 … アップロード … 再度コミットプッシュ」）→ 無認証 GET で確認（済）
4. Zenn 記事 articles/fine-structure-constant-light-cone-observation-maps.md、RELEASE_NOTES.md 追記、PAPERS_ja/en.md 47 番・docs/ 追記（全再生成はしない）、コミット・プッシュ（同日）
