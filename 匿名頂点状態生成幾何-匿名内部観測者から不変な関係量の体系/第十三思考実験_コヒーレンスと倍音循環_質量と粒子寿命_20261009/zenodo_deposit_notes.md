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

# Zenodo deposit（第十三思考実験 v2.0）

- 作成日：2026-10-10（API、v1.0 レコード 23266251 の `newversion`。Concept DOI 維持）
- Deposit ID：23280070
- Version DOI（prereserve）：10.5281/zenodo.23280070
- Concept DOI：10.5281/zenodo.23266250（不変）
- bucket：https://zenodo.org/api/files/591d799e-5cef-4f79-b936-73e57c024c01
- 状態：**published 2026-10-10**（HTTP 202。無認証 GET：doi 10.5281/zenodo.23280070、conceptdoi 10.5281/zenodo.23266250、files 14、version v2.0、publication_date 2026-10-10）

## メタデータ

- v1.0 を継承し、version v2.0、説明文に v2.0 の段落（§6.2 の経緯と結果、周期は入力した R の帰結、2 回の実行でビット単位一致、変更履歴は論文末尾）を追加、検証 zip の記述を 2 本に、キーワードに finite-order root を追加

## 進捗

1. 日本語 md に Version DOI 23280070 を埋め込み → 英訳 md（en v1.0 に v2.0 の変更箇所の訳を反映。日英で見出し 36・図 6・表の行 53・参考文献 22・引用回数が一致）→ 日英 tex/PDF（build_tex.py のファイル名と DOI を v2.0 に変更、/tmp/tex_compile/te13v2 で pandoc → lualatex ×2、ja・en 各 25 頁、LaTeX エラー 0、欠字 0、はみ出し 0）
2. 継承した v1.0 の md・tex・pdf 6 本を DELETE（204）。bucket へ PUT（201）：md・tex・pdf 各 2（v2.0）、figures の fig03〜fig06、build_tex.py（v2.0 用）。verification_R12423_exact_root_20261010.zip（58 MB、npz 2 本を含む）は PUT が 2 回とも HTTP 502（上り約 13〜16 KB/s、45 分後と 6.5 分後）で失敗したため、木原氏の判断で Zenodo には上げず、[S8] の参照先を GitHub に改めた（本文・英訳・README・RELEASE_NOTES）。fig01・fig02・verification_R070_exchange_20261009.zip は継承のまま
3. 英語 PDF の PUT が Google ドライブ上のファイルの読み込みで止まったため中断し、scratchpad にコピーしたファイルから PUT し直した
4. `検証_R12423系の第12節追加実験_20261010/` は結果を採用しなかったため Zenodo に上げていない
