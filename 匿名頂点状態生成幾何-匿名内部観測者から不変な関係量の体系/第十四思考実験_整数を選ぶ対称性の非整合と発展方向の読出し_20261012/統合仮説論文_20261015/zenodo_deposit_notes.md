# Zenodo deposit（統合仮説論文 v1.0）

- 作成日：2026-10-15（API、新規 deposit。初版なので newversion ではない。HTTP 201、1 回で成功）
- Deposit ID：23384929
- Version DOI（prereserve）：10.5281/zenodo.23384929
- Concept DOI：10.5281/zenodo.23384928（conceptrecid 23384928。公開前の API 応答では conceptdoi は None、公開後に確定）
- bucket：https://zenodo.org/api/files/a5ff9fe5-653b-4231-bf93-67afaad4e73e
- 状態：**published 2026-10-15**（publish は HTTP 202、state done・submitted True。無認証 GET：doi 10.5281/zenodo.23384929、conceptdoi 10.5281/zenodo.23384928、version v1.0、publication_date 2026-10-15、files 15、status published。Concept DOI は最新版へ 302）

## メタデータ

- upload_type publication / publication_type article / language jpn / license cc-by-4.0 / version v1.0
- creators：Kihara, Noriaki（WF System Co., Ltd.、ORCID 0009-0004-6753-4020）
- related_identifiers（すべて references、Concept DOI）：23315557（第十四）、23266250（第十三）、23199299（第十一）、23226258（第十二）、21421366（有限位数共鳴）、22112008（インフレーション機構）、22763329（時間頂点化 K₄）、22851944（設計書）
- 英題：Reading Out Time, Space, Mass, Charge, Translational Velocity and Spin from the Relations and Interference of Complex Waves without a Background Spacetime: An Integrated Hypothesis in the Order of the Thought Experiment, from Residents of a Closed Circle to Residents of a Four-Dimensional Spherical Shell, and on to α, the Four-Dimensional Lattice and E8
- 作成に使ったメタデータ JSON と応答はスクラッチ（`zenodo/metadata_create.json`、`zenodo/create_resp.json`）

## 進捗

1. DOI 取得（済）→ 日本語 md v1.0（DOI・HTML 直接ダウンロード URL・版管理、済）→ 英訳（済。日英で h2 34・h3 92・引用 161・図 10・数式 20・表 7 が一致）→ 日英 tex/PDF（済。build_tex.py は第十四のコピーで prep の引数のみ変更、/tmp/tex_compile/uhp で pandoc → lualatex ×2、ja 65 頁・en 68 頁、LaTeX エラー 0、欠字 0（≅・𝒞・∮ を $\cong$・$\mathcal{C}$・$\oint$ に置換、二本並びの表示式 3 か所を分割）、はみ出しは付録 A の表のファイル名で最大 4.26pt、図の alt を空にして見出しの重複を解消）
2. bucket へ PUT（済、15 ファイル、すべて HTTP 201、1 回で成功）：md・tex・pdf 各 2、対話型 HTML 7 本、build_tex.py、SHA256SUMS.txt（14 ファイル分）。サーバー側 MD5 が 15 本とも手元と一致
3. publish（済、2026-10-15、木原氏の指示「DOI コンセプト DOI を取って、日本語 md を修正、HTML の直接ダウンロード URL を適切に修正して、英訳、tex/pdf 作成、アップロード、zenn アーティクル作成、リリースノート更新、PAPERS 更新して、再度アップロード」）
4. Zenn 記事 articles/unified-hypothesis-closed-sphere-residents.md、RELEASE_NOTES.md 追記、PAPERS_ja/en.md 50 番・docs/index.html・docs/papers.json（indent=1 を保持）・docs/p/23384929.html 追記（全再生成はしない）、コミット・プッシュ（同日）
