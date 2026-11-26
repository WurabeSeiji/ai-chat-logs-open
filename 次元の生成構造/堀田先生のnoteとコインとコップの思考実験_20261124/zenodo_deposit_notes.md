# Zenodo deposit の作業記録（2026-11-26）

- **deposit id：** 24562868　**Version DOI：** 10.5281/zenodo.24562868　**Concept DOI：** 10.5281/zenodo.24562867（conceptrecid 24562867）
- **bucket：** https://zenodo.org/api/files/2b67236d-89bd-4b26-859a-a4bd1094d28f
- **メタデータ：** title＝英題、upload_type＝publication、publication_type＝article、creators＝Kihara, Noriaki（WF System Co., Ltd.、ORCID 0009-0004-6753-4020）、language＝jpn、license＝cc-by-4.0、version＝v1.0、publication_date＝2026-11-26、keywords＝英語版のキーワード 15 個、related_identifiers（references）＝19731593・24097611・20680269・20785539、description＝英語要旨
- **手順：** 新規 deposit 作成（API）→ prereserve DOI を日本語 md のヘッダーに記入 → 英語版の図（`make_figures_en.py`）→ 英訳 → `build_tex.py` で tex/PDF → ファイルを bucket へ PUT → サーバー側 MD5 と手元を照合 → publish

## アップロードと公開

| ファイル | MD5（手元＝サーバー） | PUT |
|---|---|---|
| コインとコップの思考実験.md | 4d8ea1b222b31cd6860ab18f961bab06 | 201 |
| coin_and_cup_thought_experiment_en.md | 141a6916e3f3888987597c08c418283c | 201 |
| コインとコップの思考実験.tex | 5c6a639b3a0103aa9e506b5a2b75f793 | 201 |
| coin_and_cup_thought_experiment_en.tex | eb07cdd0ffd4ea84ea2ebde8af6535b3 | 201 |
| コインとコップの思考実験.pdf | e2044d07f02c92cf52ad1f67b48b66f7 | 201 |
| coin_and_cup_thought_experiment_en.pdf | d8b167501d8692ca5bde45582a197be7 | 201 |
| build_tex.py | e080248d290c6f0ef5c59c30d7c2a1f4 | 201 |
| figures_and_programs.zip | 488bbf175cb9a271e1b885f561c1918b | 201 |
| SHA256SUMS | d40876370dfb8dd599e3b3d254ab07a1 | 201 |

publish HTTP 202（state done）。無認証 GET：doi 10.5281/zenodo.24562868、conceptdoi 10.5281/zenodo.24562867、version v1.0、publication_date 2026-11-26、files 9、access open。英語 PDF の download URL は 200、Concept は 302 で最新版へ転送。

`figures_and_programs.zip` は、`figures/`・`figures_en/`・`results/`・`fig_cup_coin_table.py`・`fig_cup_well_cross.py`・`fig_cup_well_cross_derived_R.py`・`make_figures_en.py` をまとめたものである。中身はすべてこのフォルダにあるので、zip そのものはリポジトリに入れていない。

## tex/PDF

pandoc → lualatex を 2 回ずつ（/tmp で実行）。日本語 41 頁・英語 45 頁、LaTeX エラー 0・欠字 0。図は SVG を rsvg-convert で PDF に変換して入れた。英語 PDF の図 17 枚が英語であることを、頁の画像で目視確認した。
