# Zenodo deposit の作業記録（2026-11-01、補遺）

- **deposit id：** 23860606　**Version DOI：** 10.5281/zenodo.23860606　**Concept DOI：** 10.5281/zenodo.23860605（conceptrecid 23860605）
- **bucket：** https://zenodo.org/api/files/77328b7b-6a9d-45ca-8b94-700027ab7c0a
- **メタデータ：** title＝英題、upload_type＝publication、publication_type＝article、creators＝Kihara, Noriaki（WF System Co., Ltd.、ORCID 0009-0004-6753-4020）、language＝jpn、license＝cc-by-4.0、version＝v1、keywords＝EPR paradox, entanglement, closure, relational readout, position and momentum, locality, no-signaling, elements of reality, emission and absorption, anonymity、related_identifiers＝isSupplementTo 23857884（本丸）、references 23585436（複素数論文）・21315735（基本公理系 v9）、description＝英語要旨
- **手順：** 新規 deposit 作成（API）→ prereserve DOI を日英 md のヘッダーと版の履歴に埋め込み → build_tex.py で tex/PDF（日本語 8 頁・英語 7 頁、エラー 0）→ bucket へ PUT → MD5 照合 → publish
- **ファイル（9）：** ja/en の md・tex・pdf、build_tex.py、figures_and_programs.zip（図表/ 一式：図 1・2 日英、make_figures.py、verify_identities.py、results/checks.txt）、SHA256SUMS.txt

## アップロードと公開（2026-11-01）
| ファイル | MD5（手元＝サーバー） | PUT |
|---|---|---|
| SHA256SUMS.txt | d7cb60147e1099ea49946354cf7e652f | 201 |
| build_tex.py | b5c0e3c0b59eab3fbdc8fde7f04baefe | 201 |
| figures_and_programs.zip | 5cf80b3695ab496b1c9c0cc630b8737b | 201 |
| paper_epr_not_a_paradox_supplement_en_v1.md | 2839d8c99e4fe9805098d14afae936eb | 201 |
| paper_epr_not_a_paradox_supplement_en_v1.pdf | 2ef1c2564237bc7fc5312103718ad331 | 201 |
| paper_epr_not_a_paradox_supplement_en_v1.tex | f70c376620e1cccbd499f2806a082e47 | 201 |
| paper_epr_not_a_paradox_supplement_ja_v1.md | 003e1c2486180127b942b60a389187d0 | 201 |
| paper_epr_not_a_paradox_supplement_ja_v1.pdf | 07d950ab43615ab6356b0313b3dd8c16 | 201 |
| paper_epr_not_a_paradox_supplement_ja_v1.tex | 6cfeed0dfa65065df7f21a78b3a795df | 201 |

curl -4 で PUT、サーバー側 MD5 は 9 本とも手元と一致。publish HTTP 202（state done）。無認証 GET：doi 10.5281/zenodo.23860606、conceptdoi 10.5281/zenodo.23860605、version v1、publication_date 2026-11-01、files 9、access open、related_identifiers＝isSupplementTo 23857884・references 23585436・21315735。日英 PDF の download URL は 200、Concept は 302 で最新版へ転送。
