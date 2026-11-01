# EPR パラドックスはパラドックスではない（補遺論文フォルダ）

**作成：** 2026-11-01　**状態：** v1 公開済み（Zenodo Version DOI [10.5281/zenodo.23860606](https://doi.org/10.5281/zenodo.23860606)、Concept DOI [10.5281/zenodo.23860605](https://doi.org/10.5281/zenodo.23860605)、本丸 23857884 の isSupplementTo）
**題名：** EPR パラドックスはパラドックスではない ― 「観測は放射と吸収の同時成立である」の補遺：位置と運動量の対を閉包の恒等式として読む ―（英題 The EPR Paradox Is Not a Paradox: A Supplement to "Observation as Simultaneous Emission and Absorption": Reading the Position–Momentum Pair as a Closure Identity）
**位置づけ：** 本丸は公開済みの「観測は放射と吸収の同時成立である」（Concept DOI 10.5281/zenodo.23857884）。本稿はその系で、新しい前提も新しい数も置かず、EPR の三段論法を閉包の恒等式と既知の破壊的読出しで読み替える。

## 主旨
対は源で一回成立した閉包 q_A＋q_B＝0。これから x_B＝−x_A と p_B＝−p_A が同時に恒等式として出る。片端の値は存在せず（片側は一様）、A の読出しで B に出来事はなく、相関は順序にも距離にも依存しない。EPR の (i) 完全相関は恒等式、(ii) 局所性は正しい、誤っているのは (iii) の帰属（関係量を片端に帰属）。破壊的読出しの存在は既知の事実（P4）として置き、機構は別の課題。「瞬時の選択」は書かない。先行（Peres、Mermin、Smerlak–Rovelli、Griffiths、Jarrett、Spekkens、Wharton–Argaman、Cramer）を正面に置き、本稿にしかないものは四つに限る。

## このフォルダの中身
| ファイル | 内容 |
|---|---|
| `paper_epr_not_a_paradox_supplement_ja_v1.md` | 日本語原稿（正本）。要旨 → §1 EPR の論法 → §2 前提（前稿の P1〜P3、閉包、位置と運動量の読み、P4 既知の事実）→ §3 恒等式（図 1・2）→ §4 読み替え → §5 先行研究 → §6 落ちる形と別の課題 → §7 区分 → 参考文献 [S1]〜[S3]・[X1]〜[X18] → 付録 A → 版の履歴 |
| `paper_epr_not_a_paradox_supplement_en_v1.md` | 英語版（翻訳。日本語版が正本。図は fig*_en） |
| `図表/make_figures.py`、`図表/figures/` | 図 1（閉包の対：対蹠、時間発展、A の読出しは B に出来事を起こさない）、図 2（順序非依存：ローレンツ変換で二つの枠）。日英、PNG・SVG・PDF |
| `図表/verify_identities.py`、`図表/results/checks.txt` | 恒等式 C1〜C4 の数値確認（固定格子、乱数なし） |
| `paper_epr_not_a_paradox_supplement_{ja,en}_v1.tex` / `.pdf` | `build_tex.py` で生成（pandoc → lualatex ×2、ローカルディスクで組版）。日本語 8 頁・英語 7 頁、LaTeX エラー 0・欠字 0 |
| `build_tex.py` | 本丸の build_tex.py の忠実コピー（入力名・題名だけ変更） |
| `figures_and_programs.zip` / `SHA256SUMS.txt` | Zenodo に上げた `図表/` 一式の zip と、Zenodo 掲載ファイルの SHA-256 |
| `zenodo_deposit_notes.md` | Zenodo deposit の作業記録 |

外部文献 16 本は Crossref で実在確認済み（2026-11-01）。
