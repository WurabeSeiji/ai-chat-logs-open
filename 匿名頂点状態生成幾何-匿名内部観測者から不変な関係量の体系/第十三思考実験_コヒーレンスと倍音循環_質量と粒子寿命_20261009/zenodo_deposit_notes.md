# Zenodo deposit（第十三思考実験 v1.0）

- 作成日：2026-10-10（API、新規 deposit。初版なので newversion ではない）
- Deposit ID：23266251
- Version DOI（prereserve）：10.5281/zenodo.23266251
- Concept DOI：10.5281/zenodo.23266250（conceptrecid 23266250）
- bucket：https://zenodo.org/api/files/e5a2cbb1-1767-4fce-87f3-8805bd19a955
- 状態：**published 2026-10-10**（HTTP 202。無認証 GET：doi 10.5281/zenodo.23266251、conceptdoi 10.5281/zenodo.23266250、files 10、version v1.0、publication_date 2026-10-10）

## メタデータ（登録済み）

- upload_type publication / publication_type article / language jpn / license cc-by-4.0 / version v1.0
- creators：Kihara, Noriaki（WF System Co., Ltd.、ORCID 0009-0004-6753-4020）
- related_identifiers（すべて references）：10.5281/zenodo.23226258（第十二思考実験 Concept）、21035831（局在波の二重スリット）、21333766（二波交換の先行実験 Concept）、21421366（有限位数共鳴 Concept）、21396760（交換重み Concept）、19902677（xyztRQ の再検討 Concept）
- 英題：The Thirteenth Thought Experiment: Why Does a Particle That Has Acquired Mass Not Dissipate and Vanish? - Coherence as a Relation, Harmonic Exchange of Localized Waves, Internal Clocks and Particle Lifetimes Considered from a Single Interaction

## 進捗

1. 日本語 md に DOI 埋め込み（済）→ 英訳 md（済）→ 日英 tex/PDF（済。`build_tex.py`（第十二のコピー。ファイル名・題・日付・DOI を変更し、`\texttt` 内の `_`・`/` の後に改行可能点、裸の URL を `\url{}` で囲む後処理と xurl を追加）で /tmp/tex_compile/te13 に生成：pandoc → lualatex ×2、ja・en 各 18 頁、LaTeX エラー 0、欠字 0、はみ出し 0。§15 の三行の箱と英語版 §14 の箱を array で改行）
2. bucket へ PUT（済、全 10 ファイル HTTP 201）：md ×2、tex ×2、pdf ×2、figures の fig01・fig02 の png、build_tex.py、verification_R070_exchange_20261009.zip（`検証_R070交換の正規化と振幅_20261009/` 一式）
3. publish（済、2026-10-10、木原氏の指示「DOI と コンセプトDOI を取得 … アップロード … コミットプッシュ」）→ 無認証 GET で確認（済）
4. Zenn 記事 articles/thought-experiment-13-harmonic-exchange-mass-lifetime.md、RELEASE_NOTES.md 追記、PAPERS_ja/en.md 48 番・docs/ 追記（全再生成はしない）、コミット・プッシュ（同日）
5. 補遺 `thought_experiment_13_appendix_…_ja_v0.1.md` と外部文献一覧 `external_references_thought_experiment_13_ja_v1.0.md` は本文に統合済みのため Zenodo には上げていない
