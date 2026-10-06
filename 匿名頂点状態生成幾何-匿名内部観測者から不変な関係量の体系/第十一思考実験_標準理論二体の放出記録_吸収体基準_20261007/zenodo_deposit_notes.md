# Zenodo deposit（第十一思考実験 v1.0）

- 作成日：2026-10-07（API、新規 deposit。初版なので newversion ではない）
- Deposit ID：23199300
- Version DOI（prereserve）：10.5281/zenodo.23199300
- Concept DOI：10.5281/zenodo.23199299（conceptrecid 23199299）
- bucket：https://zenodo.org/api/files/8f7d9cb0-79c8-4ca2-b19c-e2936e41811a
- 状態：**published 2026-10-07**（HTTP 202、state done。無認証 GET：doi 10.5281/zenodo.23199300、conceptdoi 10.5281/zenodo.23199299、files 28、version v1.0、publication_date 2026-10-07）

## メタデータ（登録済み）

- upload_type publication / publication_type article / language jpn / license cc-by-4.0 / version v1.0
- creators：Kihara, Noriaki（WF System Co., Ltd.、ORCID 0009-0004-6753-4020）
- related_identifiers（すべて references）：10.5281/zenodo.22851944（設計書）、22974631（第六思考実験 Concept）、22985299（補遺 Concept）、20060728（CP-Comp Concept）、20462569（球面投影 Concept）
- 英題：The Eleventh Thought Experiment: Emission Records of Two-Body Systems without a Background Spacetime - Computing Standard Theory (Dirac-Coulomb Bound States, 1PN Two-Body Gravity, Einstein A Coefficients) Referenced Only to the Rest Frame of an Absorber at Infinity, and Reading Out Emission-Absorption Relations for the Hydrogen Atom, Electron-Electron, Proton-Proton and the Hydrogen Molecule

## 進捗

1. 日本語 md に DOI 埋め込み（済）→ 英訳 md（済）→ 日英 tex/PDF（済。`build_tex.py` で /tmp/tex_compile に生成：pandoc 3.9 → lualatex ×2、日英とも DejaVu Serif/Sans/Sans Mono、`\ltjsetparameter{jacharrange={-2,-3}}` で記号・ギリシャ文字を Latin フォントへ。ja 32 頁、en 31 頁、LaTeX エラー 0、欠字 0、はみ出し ≤ 9.4 pt）
2. bucket へ PUT（済、2026-10-07、全 28 ファイル HTTP 201、計 9.0 MB）：md ×2、tex ×2、pdf ×2、fig01〜fig16.png、make_fig16_lightcone.py、exchange_rel_20261005.zip、exchange_general_20261006.zip、exchange_engine_20261006.zip、figures_paper_20261006.zip（__pycache__・*.csv・.DS_Store を除外）、citation_ledger_20261006.md
3. publish（済、2026-10-07、木原氏の指示「では、実行してください」）→ 無認証 GET で確認（済）
4. RELEASE_NOTES.md 追記、PAPERS_ja/en.md・docs/ へ追記（全再生成はしない）、コミット・プッシュ（同日）
