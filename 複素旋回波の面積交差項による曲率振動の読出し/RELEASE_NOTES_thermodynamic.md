---
title: 熱力学読出し論文 v1.0 リリースノート
type: release_note
modified: 2026-09-02
---

## 文書情報

**日本語タイトル**
時空を仮定しない自己無撞着複素関係系からの熱力学的読出し

**英語タイトル**
Thermodynamic Readout from a Self-Consistent Complex Relational System without Presupposing Spacetime

**著者**
Noriaki Kihara (木原 範昭)

**所属**
Independent Researcher

**ORCID**
0009-0004-6753-4020

**公開日**
2026年9月2日

**DOI（Version）**
10.5281/zenodo.22240034

**Concept DOI**
10.5281/zenodo.22240033

**ライセンス**
CC BY 4.0

**ドキュメントタイプ**
Research Note / 仮説・考察論文

---

## ファイル一覧

### 日本語版

| ファイル名 | 形式 | ページ数 | ファイルサイズ | 更新日 |
|---|---|---|---|---|
| `pre_spatiotemporal_thermodynamic_readout_from_zero_closure_ja.md` | Markdown | - | 28.8 KB | 2026-09-02 |
| `pre_spatiotemporal_thermodynamic_readout_from_zero_closure_ja.tex` | TeX (lualatex) | 39 | 31.7 KB | 2026-09-02 |
| `pre_spatiotemporal_thermodynamic_readout_from_zero_closure_ja.pdf` | PDF | 39 | 375 KB | 2026-09-02 |

### 英語版

| ファイル名 | 形式 | ページ数 | ファイルサイズ | 更新日 |
|---|---|---|---|---|
| `pre_spatiotemporal_thermodynamic_readout_from_zero_closure_en.md` | Markdown | - | 28.1 KB | 2026-09-02 |
| `pre_spatiotemporal_thermodynamic_readout_from_zero_closure_en.tex` | TeX (pdflatex) | 32 | 31.2 KB | 2026-09-02 |
| `pre_spatiotemporal_thermodynamic_readout_from_zero_closure_en.pdf` | PDF | 32 | 182 KB | 2026-09-02 |

---

## 公開プラットフォーム

- **Zenodo（正本）:** https://zenodo.org/records/22240034
  - Record ID: 22240034
  - Version DOI: 10.5281/zenodo.22240034
  - Concept DOI: 10.5281/zenodo.22240033
  - 全ファイル（MD/TeX/PDF×2言語）

- **note（日本語詳解版）:** 近日公開予定

- **Zenn（概要・図解記事）:** 
  - 日本語版：`zenn_thermodynamic_readout_ja.md` 公開済

- **GitHub（ソースコード・メタデータ）:**
  - リポジトリ: https://github.com/wurabeseiji/ai-chat-logs-open/
  - ブランチ: main
  - フォルダ: `複素旋回波の面積交差項による曲率振動の読出し/`

---

## Version History

### v1.0（初版）- 2026年9月2日

**リリース内容**

本論文は、時空・粒子・質量・エネルギー・温度を基礎概念として仮定せず、複素数で表された自己無撞着な関係系から、熱力学的量に相当する読出しを構成できるかを考察する初版である。

**中核的主張**

1. **基礎系の構成**
   - 有限個の複素状態からなる関係系
   - 自己無撞着条件としてのゼロ閉塞 $\sum_n x_n^2 = 0$
   - 有限位相周期 $U^N = I$

2. **状態数とエントロピー**
   - ゼロ閉塞の多様な分解から許容状態数 $\Omega$ を定義
   - Boltzmann 型エントロピー $S = \ln\Omega$ の導出

3. **エネルギーと温度**
   - 内部位相尺度からのエネルギー型計量の構成
   - 温度型読出し $1/T_{\mathrm{read}} = \Delta S / \Delta E_{\mathrm{read}}$ の定義

4. **観測読出し階層**
   - 読出し値は観測対象の選択に依存
   - 曲率読出し（重力波論文）と共通の観測階層に属する

5. **宇宙論的検証**
   - 離散度 $N \sim 10^{60}$ での二体関係数 $M \sim 10^{120}$
   - バリオン対光子比 $\sim 10^{-9}$ との矛盾検査（決定的検算）

**主要な発見**

- 時空・粒子なしの関係系から「エントロピー」と「温度」に同型の構造が必然的に出現
- 観測選択によって異なる読出し値は、相対性原理の示唆
- 巨大な関係空間での物質希薄性は観測宇宙の特性と矛盾しない

**前出論文との関係**

- **重力波論文（2026-09-01）** との統一性：共通の観測読出し階層
- **既往の双対幾何シリーズ** との継続：複素関係系の新たな読出し方法
- **自己無撞着インフレーション論文** との補完：マクロ系へのスケール拡張

---

## 修正履歴（git commits）

| Commit Hash | 日付 | メッセージ | 影響ファイル |
|---|---|---|---|
| （初版 Zenodo 登録） | 2026-09-02 | 熱力学読出し論文 v1.0 公開 | 全6ファイル（MD/TeX/PDF） |

---

## 技術仕様

**マークダウン標準**
- エンコーディング: UTF-8
- 改行: LF（Unix 形式）
- 数式: LaTeX フェンス（`$...$` inline、`$$...$$` display）

**TeX/PDF 生成**
- 日本語: lualatex（luatexja + ltjsarticle）、2回コンパイル
- 英語: pdflatex（lmodern）、2回コンパイル
- テンプレート: /tmp/tex_compile/ 配置
- コンパイル環境: /tmp/tex_compile/（Google Drive 外で実行）

**メタデータ**
- Zenodo フォーマット: JSON（metadata 準拠）
- Creator: ORCID 付き著者情報
- Keywords: 9項（複素関係系、ゼロ閉塞、自己無撞着、状態数など）
- License: CC BY 4.0

---

## 参考文献

**先行研究（木原）**

1. 複素旋回波の面積交差項による曲率振動の読出し（2026-09-01）
   - 重力波への最小構成
   - DOI 10.5281/zenodo.22240022 (Concept) / 22240021 (Version)

2. 自己無撞着インフレーション論文シリーズ
   - Version 2 (2026-08-30)
   - DOI 10.5281/zenodo.22176949 (Version) / 22112008 (Concept)

3. 双対幾何シリーズ（論文1-8）
   - 複素関係系の幾何的基礎

**標準参考**

- Boltzmann, L. (1877). Über die Beziehung zwischen dem zweiten Hauptsatze der Wärmetheorie und der Wahrscheinlichkeitsrechnung.
- Gibbs, J. W. (1902). Elementary Principles in Statistical Mechanics.

---

## 注記

**限定事項**

本稿は既存熱力学の完全な導出ではない。以下を示す限定的な考察である：

- 時空・粒子を仮定しない複素関係系からの読出し構造
- 自己無撞着性の必然性
- 観測階層の相対性

未解決課題：動力学機構、相互作用、第二法則、完全統一フレーム。

**読者層**

- 理論物理学者（既成枠組みの批判的検討に関心のある研究者）
- 数学者（複素関係系の組合せ論的側面に関心のある者）
- 哲学者（基礎的仮定の必要性を問う者）

**引用形式**

BibTeX:
```
@article{Kihara2026Thermodynamic,
  author={Kihara, Noriaki},
  title={Thermodynamic Readout from a Self-Consistent Complex Relational System without Presupposing Spacetime},
  journal={Zenodo},
  year={2026},
  month={9},
  day={2},
  doi={10.5281/zenodo.22240034}
}
```

---

## ライセンス

CC BY 4.0 — 帰属を明記すれば、自由に複製・改変・配布可。

---

**最終確認日**: 2026年9月2日 UTC
**確認者**: Noriaki Kihara, Independent Researcher
