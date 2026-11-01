# Zenodo deposit の作業記録（2026-11-01）

- **deposit id：** 23857885　**Version DOI：** 10.5281/zenodo.23857885　**Concept DOI：** 10.5281/zenodo.23857884（conceptrecid 23857884）
- **bucket：** https://zenodo.org/api/files/01cdffa5-7fb8-4fe3-ac15-03fd37767537
- **メタデータ：** title＝英題、upload_type＝publication、publication_type＝article、creators＝Kihara, Noriaki（WF System Co., Ltd.、ORCID 0009-0004-6753-4020）、language＝jpn、license＝cc-by-4.0、version＝v1.4、keywords＝observation, emission and absorption, handshake, readout rule, viewing angle, Born rule, Bell inequality, CHSH, Tsirelson bound, no-signaling, absorber theory, transactional interpretation, harmonics, anonymity、related_identifiers（references）＝23585436（複素数論文）・22909660（第三思考実験）・22851944（研究プロジェクト設計）・21315735（基本公理系 v9）・23199299（第十一思考実験）、description＝英語要旨
- **手順：** 新規 deposit 作成（API）→ prereserve DOI を日英 md のヘッダーと版の履歴に埋め込み → build_tex.py で tex/PDF → ファイルを bucket へ PUT → サーバー側 MD5 と手元を照合 → publish
- **ファイル（9）：** ja/en の md・tex・pdf、build_tex.py、figures_and_programs.zip（図表/ 一式）、SHA256SUMS.txt

## アップロードと公開（2026-11-01）
| ファイル | MD5（手元＝サーバー） | PUT |
|---|---|---|
| paper_observation_handshake_bell_tsirelson_ja_v1.4.md | 3212127130e704ffb9a7421bbe40c306 | 201 |
| paper_observation_handshake_bell_tsirelson_en_v1.4.md | ee6306f330b8867b4c3570b8bd727bcf | 201 |
| paper_observation_handshake_bell_tsirelson_ja_v1.4.tex | edcc2cbbe1366b7079d9ecf1a902e37d | 201 |
| paper_observation_handshake_bell_tsirelson_en_v1.4.tex | 45894d5bf24bb7ce92adc6dca1d24ff4 | 201 |
| paper_observation_handshake_bell_tsirelson_ja_v1.4.pdf | 639c13aa0dcce279d679918a0444a1ab | 201 |
| paper_observation_handshake_bell_tsirelson_en_v1.4.pdf | 5907ee9843bbfc8f04eef68a84a1c32d | 201 |
| build_tex.py | d271de15cca365e013e2b5e52a611b06 | 201 |
| figures_and_programs.zip | 06dee057b68e62f1e3c7882a4aa0deda | 201 |
| SHA256SUMS.txt | 84cd956dbbdf70ef4e258b8e8f7871a1 | 201 |

curl -4 で PUT、サーバー側 MD5 は 9 本とも手元と一致。publish HTTP 202（state done）。無認証 GET：doi 10.5281/zenodo.23857885、conceptdoi 10.5281/zenodo.23857884、version v1.4、publication_date 2026-11-01、files 9、access open。日英 PDF の download URL は 200、Concept は 302 で最新版へ転送。
