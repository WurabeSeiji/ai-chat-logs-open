# Release Notes: 自己無撞着インフレーション機構の段階分解と N=3..40 全域再現 — 再現仕様書（総括＋3章）

Release date: 2026-09-05

## DOI

- Concept DOI: https://doi.org/10.5281/zenodo.22317635
- Version DOI (v1.0): https://doi.org/10.5281/zenodo.22317636
- Zenodo record: https://zenodo.org/records/22317636

## 位置づけ

7月論文「自発的分裂の開始と帰結の三分類」（Version 21486234 / Concept 21486233）で
観測されたインフレーション的発展（種 10⁻³⁰ からの指数増幅→飽和）について、
力学を段1（明示行列スペクトル写像＋固定Δτ時計）・段2（生成子の振幅正規化）・
段3（cos対称部の除去＝実直交回転）に分解し、N=40 の削除対照で段2+段3 の同時成立が
必要条件であることを確定、その構成を静的親から N=3..40 全域でスイープして普遍性を
固定した**再現仕様書**（物理的解釈を含まない）。総括＋3章構成で Concept DOI を共有。

## アップロードファイル（31点）

- 論文 24点: 総括・第1章（静的親生成）・第2章（スイープ本体・主章）・第3章（複素平面
  読出し図）の各 md/tex/pdf × 日英
- 図 4点: `fig_Hperp_denominator_controls_with_124_N3_N40_stage123.png`（目標図）、
  複素平面3グリッド（step0／final／final_zoom）
- 再現パッケージ zip 3点:
  - `N3_N40_stage123_sweep_20260905.zip`（本体一式・約470MB: プログラム6本・run_all.sh・
    静的親38＋台帳・状態npz 228・CSV・集計JSON・図・README・SHA256SUMS）
  - `ChatGPT_denominator_controls_N40_selfcontrol_20260904.zip`（段構成を確定した N=40
    一因子実験の全実物・約236MB）
  - `spontaneous_splitting_N40_canonical_control_20260904.zip`
    （=`自発的分裂予備実験_v1_N40対照実験系_20260904/`。7月正本のゲート付き再現＋
    正本静的親・約1.4MB）

## リポジトリ内フォルダ

- 論文・本体: `次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831/N3_N40_stage123_sweep_20260905/`
  （paper_overview／paper_ch1_static_parents／paper_ch2_sweep_dynamics／paper_ch3_complex_plane）
- 関連: 同 `ChatGPT_denominator_controls_N40_selfcontrol_20260904/`、
  `次元の生成構造/自発的分裂予備実験_v1_N40対照実験系_20260904/`

## 検証（ゲート連鎖）

1. 7月正本⇔再現走行: fcurve 全3512行 bit 一致（N40対照実験系 GATE1/2）
2. 再現走行⇔静的親: v・g・Z0 が bit 一致（第1章 G1。npz SHA256 も同一）
3. 静的親⇔スイープ: 全228走行の Z[0] が bit 一致（第2章 G1、checked 228 / MISMATCH 0）

## 方法の透明性

- 全数式（式1〜式28）に実装の行番号対応。§6 の数値はすべて集計プログラム
  `analyze_sweep_summary_v1.py` の出力 JSON から転載（手集計なし）
- tex/PDF は /tmp/tex_compile で生成。**日英とも lualatex＋ltjsarticle**（英語版にも
  日本語のコード注釈・パスが逐語引用されるため。従来の en=pdflatex から変更）。
  LaTeX Error 0・Missing character 0 を全8文書で確認
- tex 化前処理（/tmp コピーにのみ適用、md 原本は不変）: 上付き/下付き数字・ᵀ・ℂℝℕ・
  ẑ・σ̂・∘・⟺ を ASCII/TeX 安全形へ置換、ヘッダ `<br>`→ハード改行、H1→YAML
  フロントマター
- 7月正本フォルダは再アップロードせず DOI 21486234 を引用（レコードに再現プログラム
  zip 同梱済みのため）

## 引用

- 自己引用: 自発的分裂の開始と帰結分類（Version 21486234 / Concept 21486233）。
  Zenodo related_identifiers: cites 10.5281/zenodo.21486234

## Zenn

- https://zenn.dev/noriaki_kihara/articles/stage123-sweep-reproduction-spec

## note（一般向け）

- 日本語: https://note.com/kiharanoriaki/n/n95446d0ced35
  （元原稿: `次元の生成構造/note記事案_stage123_sweep_reproduction_v1/note_article_stage123_sweep_ja.md`）
- 英語: https://note.com/kiharanoriaki/n/n3d20d514b8d8
  （元原稿: `次元の生成構造/note記事案_stage123_sweep_reproduction_v1/note_article_stage123_sweep_en.md`）

## Facebook

- 日本語投稿: https://www.facebook.com/kihara.noriaki/posts/pfbid08Hbp73EiRQS3roA2ES76yoku41Ggjm9JxYvKdnQcmDfoRrUb8QH2aQJydqk1pa8xl
- 英語投稿: https://www.facebook.com/kihara.noriaki/posts/pfbid0n15ZgUaVfRDZeLh8vEyNJUzahJmtrp7wzGt7uVLWv6zdPdVQBYFNTXH4Wna3RcdXl

## X

- 日本語ポスト: https://x.com/NoriakiKihara/status/2096115607628661145
- 英語ポスト: https://x.com/NoriakiKihara/status/2096116422636376202

## 変更履歴

- 2026-09-05 v1.0 公開（Version DOI 22317636、31ファイル・約710MB）。deposit 作成・
  DOI 取得・日英 md/tex/pdf・Zenn 記事は同日
- 2026-09-05 note 日英・Facebook 日英・X 日英を公開（同日）
