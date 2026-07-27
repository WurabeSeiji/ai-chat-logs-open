# Release Notes: 三方向生成の時間構造の因果分離——二段階seed除去（次元の生成構造 第8論文）

Release date: 2026-07-27

## DOI

- Concept DOI: https://doi.org/10.5281/zenodo.21614402
- Version DOI (v1): https://doi.org/10.5281/zenodo.21614403
- Zenodo record: https://zenodo.org/record/21614403

## note（一般向け）

- 日本語: https://note.com/kiharanoriaki/n/n48a02cd70f47 （2026-07-27 公開）
  （元原稿: `次元の生成構造/タネがなくてもインフレーションは起きた_note完成稿.md`）
- 英語: https://note.com/kiharanoriaki/n/nb584455b0aa5 （2026-07-27 公開）
  （元原稿: `次元の生成構造/inflation_without_a_seed_note_en.md`）

## Facebook

- 日本語投稿: https://www.facebook.com/kihara.noriaki/posts/pfbid0362xhGjsbjGjYoJ2RffBr5YcnCQDUbDuyDyDM3doHYr9qL82vhKthFxBkE7pHsqgMl （2026-07-27、日本語note URLはコメント欄）
- 英語投稿: https://www.facebook.com/kihara.noriaki/posts/pfbid02uTCZukGCJ19oxfnNga3poC9YHjPnMtXFHQfRFaUUGNN3QsjRr1mVytFiMywcmdyHl （2026-07-27、英語note URLはコメント欄）

## X

- 日本語ポスト: https://x.com/NoriakiKihara/status/2081580195061104704 （2026-07-27）
- 英語ポスト: https://x.com/NoriakiKihara/status/2081580633030365211 （2026-07-27）

## 公開ファイル（Zenodo, 7点）

- `nbody_two_stage_seed_removal_causal_separation_ja.md` / `_en.md`
- `nbody_two_stage_seed_removal_causal_separation_ja.tex` / `_en.tex`
- `nbody_two_stage_seed_removal_causal_separation_ja.pdf` / `_en.pdf`
- `nbody_two_stage_seed_removal_reproduction_v1.zip`（原本6本＋計装32本＋図9枚＋報告書11件＋表11件＋実行メタ8件、MANIFEST付き）

## 主張

第7論文が発見した三方向生成（長い低変化領域→幾何級数的急拡大→三方向的準安定構造の自発形成）の時間構造を、二段階の明示的seedを独立除去して因果分離した。seedは分裂の発生・時刻・最初の生成方向・準安定振動のいずれの原因でもない（crossing差 ≤5 step、同一step方向overlap中央値 0.999999987）。潜伏領域は毎step検査で負差分ゼロの単一連続増大であり下位の底を持たない（λ=0.04937, R²=0.9999993）。rank_q=4 の早期出現は数値床応答であり方向成立と同一視できない。方向部分空間は急拡大中に回転・混合して再編成される（早期対後期overlap 0.14）。準安定状態は t=110000 まで第二の急拡大を生まない。初期二方向状態は零二乗閉鎖の帰結（全M関係波非零のまま実rank 2）であり、実験者の置いた二軸ではない。ボゾン的／フェルミオン的区別の自発形成の予兆は観測量に現れない（残余回転占有は全Nで数値零へ減衰）。

## 方法の透明性

- 第7論文原本を SHA-256 固定・read-only import で再利用（原本無変更）。第7論文軌道のビット一致再現（Stage A0）、無seed軌道の独立2実行ビット一致（Stage A2a）、N=300 五色分解の既存条件Aとの共通27列全行一致を確認。
- rank_q は Q=[B₀|B_dom]（M×4）の rank で構造上4が上限のため追加方向の検出器にならないことを明記し、追加方向不在の判定は五色分解の残余回転占有の直接測定による。
- 生CSVは決定論的に再生成可能なため版管理外・zip非同梱。各CSVの SHA-256 は実行メタ JSON に記録。

## 引用

- 自己引用: 基本公理系（Concept DOI 21315735）、第7論文（21543070）、零二乗和解説ノート（21495305）、局在性交換（21333766）、白猫黒猫灰色猫（21353208）
- 外部: Benettin ほか（FPUT準安定）、Mirollo–Strogatz（Kuramotoロック状態スペクトル）、Daumont–Dauxois–Peyrard（変調不安定性）、Mori ほか（前熱化）、Kim–Nishimura–Tsuchiya（IKKT行列模型 (3+1)次元創発）

## 生成手順

- 数式は GitHub 記法（$$ / $）。tex/PDF は /tmp/tex_compile で生成（日本語=ltjsarticle+lualatex 2回、英語=article+pdflatex 2回、TeX/LaTeX Error 0・Missing character 0）。
- 前処理: 日本語版は Unicode 上付き・下付き（×10⁻ⁿ, qₙ 等）を数式へ変換（Latin Modern にグリフが無く欠落するため）。英語版は ≤ → ` $\le$ `（後続数字と密着させない。密着すると pandoc が数式と認識しない）、→・π・− を数式/ASCIIへ。
- CSV は .gitignore（*.csv）により追跡外。
