# Zenodo deposit の作業記録（2026-11-10）

- **deposit id：** 24097612　**Version DOI：** 10.5281/zenodo.24097612　**Concept DOI：** 10.5281/zenodo.24097611（conceptrecid 24097611）
- **bucket：** https://zenodo.org/api/files/9d41574b-8100-489b-9d9a-ee3c24e8d75f
- **メタデータ：** title＝英題、upload_type＝publication、publication_type＝article、creators＝Kihara, Noriaki（WF System Co., Ltd.、ORCID 0009-0004-6753-4020）、language＝jpn、license＝cc-by-4.0、version＝v1.0、publication_date＝2026-11-10、keywords＝observer resolution, gauge, quaternion, harmonic order, winding number, phase, Pontryagin duality, finite Weyl relation, Peter-Weyl theorem, finite subgroups of SU(2), Niven's theorem, cyclotomic field, McKay correspondence, extra dimension, viewing angle, thought experiment、related_identifiers（references）＝23451494（四元数の写像）・23444451（対称性とマジックナンバーの目録）・23857884（観測は放射と吸収の同時成立である）・23585436（複素数論文）、description＝英語要旨
- **手順：** 新規 deposit 作成（API）→ prereserve DOI を日本語 md のヘッダーと版の履歴に埋め込み → 英訳 → build_tex.py で tex/PDF → ファイルを bucket へ PUT → サーバー側 MD5 と手元を照合 → publish
- **ファイル（9）：** ja/en の md・tex・pdf、build_tex.py、figures_and_programs.zip（図表/・概念図_次数を逆に読む/・run_all.sh）、SHA256SUMS

## アップロードと公開（2026-11-10）
| ファイル | MD5（手元＝サーバー） | PUT |
|---|---|---|
| paper_resolution_gauge_reverse_reading_ja_v1.0.md | 34e3fc37bd093d2a1a5fcebb2a542c67 | 201 |
| paper_resolution_gauge_reverse_reading_en_v1.0.md | 304e02dde8093c6b43c4608ba24fbbe5 | 201 |
| paper_resolution_gauge_reverse_reading_ja_v1.0.tex | ea1bbed870d9843c3221658e3a859fad | 201 |
| paper_resolution_gauge_reverse_reading_en_v1.0.tex | 3d52073ccf41e719d58aeeb466b98507 | 201 |
| paper_resolution_gauge_reverse_reading_ja_v1.0.pdf | 1d7b17aa6090cc714a8ddd171f44ad8a | 201 |
| paper_resolution_gauge_reverse_reading_en_v1.0.pdf | 2362662151017bdcb8e75f453304ae20 | 201 |
| build_tex.py | a22f547b49e9d955ae7124655396d82b | 201 |
| figures_and_programs.zip | 2f393afbc73997de51a1cd230c82def3 | 201 |
| SHA256SUMS | 6fafeef1f8aeb537e591ddbbd67d6e4a | 201 |

curl -4 で PUT、サーバー側 MD5 は 9 本とも手元と一致。publish HTTP 202（state done）。無認証 GET：doi 10.5281/zenodo.24097612、conceptdoi 10.5281/zenodo.24097611、version v1.0、publication_date 2026-11-10、files 9、access open。日英 PDF の download URL は 200、Concept は 302 で最新版へ転送。

## tex/PDF

pandoc → lualatex を 2 回ずつ（/tmp で実行）。日本語 15 頁・英語 15 頁、LaTeX エラー 0・欠字 0・はみ出し 0。§11 の交差項の式は 1 行だと版面を 76pt はみ出したので、日英の md で aligned の 2 行に分けた。表の列幅は md の区切り行のダッシュの比で指定した。英語 PDF の図 3 が英語であることを、頁の画像で目視確認した。図 1・2 は文字のない切り抜きで、日英共通。
