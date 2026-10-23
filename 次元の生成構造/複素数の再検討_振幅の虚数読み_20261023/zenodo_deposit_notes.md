# Zenodo deposit（複素数を正しく理解していたか v1）

- 作成日：2026-10-23（API、新規 deposit。HTTP 201）
- Deposit ID：23585437
- Version DOI（prereserve）：10.5281/zenodo.23585437
- Concept DOI：10.5281/zenodo.23585436（conceptrecid 23585436）
- 状態：**published 2026-10-23**（HTTP 202、files 12、version v1）
- bucket：https://zenodo.org/api/files/d9e82b41-9920-4447-9557-de0c47988fb2
- 指示：「DOI コンセプト doi 取得　日本語 md への doi 登録　英訳　Tex pdf 作成　アップロード　zen アーティクル作成　リリースノート　Papers 更新　コミットプッシュを実行してください」（2026-10-23）
- 版番号：v1（日本語原稿 v1 を acfe4975 で公開前稿としてプッシュしたあと、同じ v1 に DOI を記入して Zenodo の初版とした）

## メタデータ

- upload_type publication / publication_type article / publication_date 2026-10-23 / language jpn / license cc-by-4.0 / version v1
- creators：Kihara, Noriaki（WF System Co., Ltd.、ORCID 0009-0004-6753-4020）
- 英題：Did We Understand Complex Numbers Correctly? Where the Amplitude Lives, the Imaginary Reading, and Attribution across ℂ → ℍ → 𝕆
- keywords：complex numbers, quaternions, octonions, amplitude and phase, duality, imaginary reading, rapidity, Born rule, Lorentz factor, modified Rodrigues parameters, associator, anonymity, partial projection, attribution
- related_identifiers（すべて references、Concept DOI）：23523237（結晶化とインフレーション、前稿）、23451494（複素数の集合から時空への写像）、23444451（目録）、23384928（統合仮説論文）、21315735（基本公理系 v9）、20521566（位相空間の面積としての保存量と Wick 回転）、20967081（観測者–系ビートによるボルン統計）、22282689（複素定常波真空からの凝縮体の読出し）
- 作成に使ったメタデータ JSON と応答はスクラッチ（`zenodo_cri/metadata_create.json`、`zenodo_cri/create_resp.json`）

## 進捗

1. DOI 取得（済）→ 日本語 md v1 のヘッダーと版の履歴に Version DOI・Concept DOI を記入（済）
2. 英訳 `complex_revisited_amplitude_imaginary_en_v1.md`（済。日本語版が正本で英語版は翻訳と明記。図は `図表/figures/fig*_en.svg`。自己引用の題名は Zenodo の英題（docs/papers.json）、外部文献は日本語版と同じ）
3. 日英 tex/PDF（済。`build_tex.py` は前稿のものを元に、図の参照先 `図表/figures/*.svg`・DOI のヘッダー読み・黒板太字 ℂ ℍ 𝕆 の \mathbb 置換・参考文献の段落分けを加えた。/tmp/tex_compile/cri で pandoc → lualatex ×2、日本語 23 頁・英語 24 頁、LaTeX エラー 0・欠字 0。1 頁目・図のある頁・参考文献の頁を PNG にして目視確認）
4. bucket へ PUT（済、12 ファイル、curl -4、サーバー側 MD5 は 12 本とも手元と一致）：日英 md・tex・PDF、build_tex.py、figures_and_programs.zip（`図表/` 一式）、自己引用の対応表、外部引用の対応表、課題の整理、SHA256SUMS.txt
5. publish（済、2026-10-23、HTTP 202、state done。無認証 GET：doi 10.5281/zenodo.23585437、conceptdoi 10.5281/zenodo.23585436、version v1、publication_date 2026-10-23、files 12。PDF 日英の download URL は 200、Concept の API は 302 で最新版へ）
6. Zenn 記事 `articles/complex-numbers-revisited-amplitude.md`、RELEASE_NOTES.md、PAPERS_ja/en.md の 54 番・docs/index.html・docs/papers.json・docs/p/23585437.html の追記（済、全再生成はしていない）
7. コミット・プッシュ（公開直後に実施）
