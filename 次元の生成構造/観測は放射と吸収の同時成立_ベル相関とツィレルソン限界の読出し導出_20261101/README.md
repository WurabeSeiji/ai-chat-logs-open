# 観測は放射と吸収の同時成立である（論文フォルダ）

**作成：** 2026-11-01　**状態：** v1.4 公開済み（Zenodo Version DOI [10.5281/zenodo.23857885](https://doi.org/10.5281/zenodo.23857885)、Concept DOI [10.5281/zenodo.23857884](https://doi.org/10.5281/zenodo.23857884)）
**題名：** 観測は放射と吸収の同時成立である ― 視野角と大きさの読出し規則からボルン則・ベル相関・ツィレルソン限界を導く ―（英題 Observation as Simultaneous Emission and Absorption: Deriving the Born Rule, Bell Correlations and the Tsirelson Bound from a Viewing-Angle Readout Rule）
**位置づけ：** 前稿「複素数を正しく理解していたか」（Concept DOI 10.5281/zenodo.23585436）の次。読出しだけで閉じた論文で、作用の形と数値実験（E0〜E5）に依存しない。初稿は 2026-11-01 の新横浜〜東京の車中で Claude.ai と作成。

## 主旨
観測に混在する三つの相互作用（離散状態の切替・位相だけのずれ・位相関係の切断）を分け、観測の本体を二番目に置く。前提は P1（距離＝視野角と大きさの比、比しか読めない）・P2（放射と吸収は一回の関係、光速は換算係数）・P3（関係の位相は倍音次数 n）の三つだけ。これからボルン則の形、E＝∓cos nθ、CHSH の上限 2√2、信号不能性が出る。確率に読む正規化は大きさの比を落とす割り算＝η＝0 の断面の選択（前稿 §4.4 と整合）。倍音があれば相関はパルス形の自己相関で、上限は不変、等号は単音のみ、矩形パルス（M→∞）は局所隠れ変数で S≤2。

## このフォルダの中身

| ファイル | 内容 |
|---|---|
| `paper_observation_handshake_bell_tsirelson_ja_v1.4.md` | 日本語原稿 v1.4（正本）。要旨 → §1 三つの相互作用 → §2 前提 P1〜P3（図 1〜3）→ §3 導出（3.1 正規化、3.2 ボルン則、3.3 二者相関、3.4 ツィレルソン限界（図 4）、3.5 信号不能性）→ §4 帰結と照合（図 5）→ §5 先行研究 → §6 議論 → §7 区分 → 参考文献 [S1]〜[S5]・[X1]〜[X24] → 付録 A 図の再現 → 版の履歴 |
| `paper_observation_handshake_bell_tsirelson_en_v1.4.md` | 英語版（翻訳。日本語版が正本。図は `図表/figures/fig*_en.svg`、自己引用に英題を併記） |
| `paper_observation_handshake_bell_tsirelson_{ja,en}_v1.4.tex` / `.pdf` | `build_tex.py` で生成（pandoc → lualatex ×2、ローカルディスクで組版）。日本語 11 頁・英語 11 頁、LaTeX エラー 0・欠字 0 |
| `build_tex.py` | tex 生成スクリプト（前稿のものの忠実コピー。入力ファイル名・題名・著者行・日付と、英語版の図説見出し「Fig. N」の正規表現だけを変更） |
| `図表/make_figures.py`、`図表/figures/` | 図 1〜5（日英、PNG・SVG・PDF）。点はすべて式から計算、乱数なし。図 3 は有限和 N＝13 |
| `figures_and_programs.zip` / `SHA256SUMS.txt` | Zenodo に上げた `図表/` 一式の zip と、Zenodo 掲載ファイルの SHA-256 |
| `zenodo_deposit_notes.md` | Zenodo deposit の作業記録（メタデータ・ファイル・進捗） |
| `paper_observation_handshake_bell_tsirelson_ja_v1.md` / `_v1.1.md` / `_v1.2.md` / `_v1.3.md` | 草稿の履歴（v1 初稿、v1.1 引用番号、v1.2 修正点適用、v1.3 図の追加）。版の履歴は本文末尾 |
| `修正点_自己分析_v1.1_20261101.md` | v1.1 → v1.2 の修正点（箇所・現行・置換・理由、A 必須 4・B 必須 7・C 追加 3・D 1・E 1・F 1） |
| `申し送り_ClaudeCode_引用作成とベル数値実験E0-E5_20261101.md` | 課題 A（引用作成、完了）と課題 B（ベル数値実験 E0〜E5、**今後の課題**。3 パルス系の拡張、ボルン則を手で入れない読出し規則） |

## 再現
図：`cd 図表 && python3 make_figures.py && python3 make_figures.py en`。tex：ローカルの作業フォルダで `python3 build_tex.py <このフォルダ>` → `lualatex ja.tex` ×2、`lualatex en.tex` ×2。
