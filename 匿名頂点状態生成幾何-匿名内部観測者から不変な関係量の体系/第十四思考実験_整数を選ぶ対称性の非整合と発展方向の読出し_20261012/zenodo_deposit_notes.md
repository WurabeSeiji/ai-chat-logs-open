# Zenodo deposit（第十四思考実験 v1.0）

- 作成日：2026-10-12（API、新規 deposit。初版なので newversion ではない）。1 回目の作成要求は Zenodo 側の HTTP 504 で失敗し、再試行で作成
- Deposit ID：23315558
- Version DOI（prereserve）：10.5281/zenodo.23315558
- Concept DOI：10.5281/zenodo.23315557（conceptrecid 23315557）
- bucket：https://zenodo.org/api/files/34cf7785-aa64-4920-85e6-a9bdf0e33356
- 状態：**published 2026-10-12**（publish 要求は HTTP 504 だったが、deposit の state=done・submitted=True を確認。無認証 GET：doi 10.5281/zenodo.23315558、conceptdoi 10.5281/zenodo.23315557、files 12、version v1.0、publication_date 2026-10-12。Concept DOI は最新版へ 302）

## メタデータ

- upload_type publication / publication_type article / language jpn / license cc-by-4.0 / version v1.0
- creators：Kihara, Noriaki（WF System Co., Ltd.、ORCID 0009-0004-6753-4020）
- related_identifiers（すべて references）：10.5281/zenodo.23226258（第十二思考実験 Concept）、10.5281/zenodo.23266250（第十三思考実験 Concept）、10.5281/zenodo.21421366（有限位数共鳴 Concept）
- 英題：The Fourteenth Thought Experiment: Why Can Different Symmetries Never Be Fully Reconciled? - Non-commensurability, Pull-in and Metastability of Integer-Selecting Symmetries, and Readout along the Direction of Evolution

## 進捗

1. 日本語 md に DOI 埋め込み（済）→ 英訳 md（済。サブエージェントは書き出し前に停止したため Claude が直接訳出。日英で行数 743・見出し 46・表の行 51・数式ブロック 62・参考文献 30・式番号の並びが一致）→ 日英 tex/PDF（済。build_tex.py は第十三のコピーでファイル名・題・日付・DOI の 8 行のみ変更、/tmp/tex_compile/te14 で pandoc → lualatex ×2、ja 24 頁・en 25 頁、LaTeX エラー 0、欠字 0、はみ出し 0）
2. bucket へ PUT（済、12 ファイル）：md・tex・pdf 各 2、build_tex.py、verify_te14.py・verify_phase_lattice.py・verify_e8.py・verify_null_from_products.py、SHA256SUMS.txt（11 ファイル分）。en tex と ja pdf の 2 本は 1 回目が HTTP 504 で、再送で 201。サーバー側の MD5 が 12 本とも手元と一致することを確認。/tmp のコピーから送信
3. publish（済、2026-10-12、木原氏の指示「DOI コンセプトDOI を取得 … 英訳、tex/pdf 作成、アップロード、zenn アーティクル作成、リリースノート更新、PAPERS 更新、コミットプッシュ」）
4. Zenn 記事 articles/thought-experiment-14-competing-integer-symmetries.md、RELEASE_NOTES.md 追記、PAPERS_ja/en.md 49 番・docs/index.html・docs/papers.json（元の indent=1 を保持）・docs/p/23315558.html 追記（全再生成はしない）、コミット・プッシュ（同日）
