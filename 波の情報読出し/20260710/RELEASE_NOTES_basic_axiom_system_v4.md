# RELEASE NOTES: 無名等振幅複合波モデル基本公理系 v4

**日付:** 2026-07-12
**対象:** 無名等振幅複合波モデル基本公理系 v4
**Version DOI:** 10.5281/zenodo.21316620
**Concept DOI:** 10.5281/zenodo.21315735
**Zenn:** https://zenn.dev/noriaki_kihara/articles/basic-axiom-system-v4

---

## v4 公開内容

本リリースでは、基本公理系 v4 を次の二層で公開した。

| 区分 | 内容 |
|---|---|
| 純化定義論文 | 公理、帰結、作業公理、定義のみを含む正本 |
| 解釈メモ | 現時点での解釈、説明、作業仮説、読出し例を含む補助文書 |

---

## 主要変更

1. 公理3「投影軸縮退と虚方向読出し」を追加した。
2. 投影法における投影軸方向は、内在座標系の住民には方向として区別できないと明記した。
3. 投影軸方向の値は消去されるのではなく、内在座標系では虚方向値として扱うことを定義した。
4. 補償項 `z` を `iz` として読み、二乗により `(iz)^2=-z^2` として観測可能方向の読出しに現れることを公理層で整理した。
5. 旧 v3 の公理3以降を、公理4以降へ繰り下げた。
6. 純化定義論文と解釈メモの英訳、TeX、PDFを作成した。

---

## 公開ファイル

| 種類 | ファイル |
|---|---|
| 純化定義論文 JA | `無名等振幅複合波モデル基本公理系v4_純化定義論文.md` |
| 純化定義論文 EN | `anonymous_equal_amplitude_composite_wave_model_basic_axiom_system_v4_pure_definition_en.md` |
| 解釈メモ JA | `基本公理系 v4.md` |
| 解釈メモ EN | `basic_axiom_system_v4_en.md` |
| 純化定義論文 JA TeX | `basic_axiom_system_v4_pure_definition_ja.tex` |
| 純化定義論文 EN TeX | `basic_axiom_system_v4_pure_definition_en.tex` |
| 解釈メモ JA TeX | `basic_axiom_system_v4_interpretation_note_ja.tex` |
| 解釈メモ EN TeX | `basic_axiom_system_v4_interpretation_note_en.tex` |
| 純化定義論文 JA PDF | `basic_axiom_system_v4_pure_definition_ja.pdf` |
| 純化定義論文 EN PDF | `basic_axiom_system_v4_pure_definition_en.pdf` |
| 解釈メモ JA PDF | `basic_axiom_system_v4_interpretation_note_ja.pdf` |
| 解釈メモ EN PDF | `basic_axiom_system_v4_interpretation_note_en.pdf` |

---

## 備考

純化定義論文は、解釈や応用説明を入れず、公理系そのものの定義に限定した。

解釈メモは、今後拡張または変更される可能性がある。公理そのものに誤り、欠落、または変更の必要が発見された場合は、解釈メモで読み替えるのではなく、正本の公理系を明示的に改訂する。
