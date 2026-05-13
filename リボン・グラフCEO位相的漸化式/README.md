# リボン・グラフ／CEO位相的漸化式 学習資料

## 用途

2026-05-30（土）14:00–18:00、大阪公立大学 森之宮キャンパス 中教室404にて開催される **南部研究所大阪城セミナー 第2回**「ファットグラフに基づくCEO位相的漸化式入門」（藤博之・神戸大学CMDS）の予習用文献。

聴講参加を念頭に、約2週間で基本概念を習得することを目標に文献を選定。

---

## 第1部：概念マップ — 「なぜこの分野が重要か」

### 1. ファットグラフ（リボン・グラフ）とは何か

通常のグラフ（点と線）に対し、**ファットグラフ／リボン・グラフ**は各辺を「リボン（帯）」として描き、各頂点で辺の **巡回順序** を指定したグラフです。これが何を生むかというと：

- リボンを境界とする曲面（**有向リーマン面**）が一意に決まる
- グラフの種数 g（穴の数）、面の数 F、頂点数 V、辺数 E の間に Euler 関係 **V − E + F = 2 − 2g** が成立
- リボン・グラフ自体が **コンパクト・リーマン面の cell 分割（三角化）** を与える

つまりファットグラフは **「組合せ論的データから位相的曲面を生成する装置」** であり、純粋に組合せ論的・離散的な対象を扱いながら、連続的なリーマン面のモジュライ空間に接続します。

### 2. 't Hooft（1974）— 大Nゲージ理論からファットグラフへ

't Hooft は U(N) ゲージ理論の Feynman 図を、グルーオン場 A^a_b を「二重線」（色添字 a と b に対応する2本の線）として描くと、図の **トポロジー（種数）が 1/N 展開の次数を決定** することを発見しました。具体的には：

- 平面図（種数0）が支配的：N² × N^(V−E) × ...
- 種数 g の図は N^(−2g) で抑制される

これがファットグラフ＝リボングラフの物理的起源です。**色付きQCDの大N極限が、リーマン面の幾何に変換される**という驚異的な結果。弦理論／行列模型の出発点。

### 3. 行列模型 — 組合せ論の生成関数

エルミート行列 M（N×N）上の積分

$$Z_N = \int dM \exp\left(-\frac{N}{g}\text{Tr}\, V(M)\right)$$

を Feynman 図展開すると、**ファットグラフの数え上げ生成関数** になります（BIPZ 1978）。具体的には、4次相互作用 V(M) = M⁴ なら、自由エネルギーは 4-valent ファットグラフを種数別に数え上げる。

これにより：
- 行列積分の漸近展開 = リーマン面のモジュライ空間体積の計算
- 行列模型は **離散的曲面（三角形分割）の重み付き和** として解釈できる

### 4. Kontsevich（1992）— Witten 予想の証明

Witten は「モジュライ空間 M̄_{g,n} 上の psi-class の交点数の生成関数が KdV 階層のτ関数になる」と予想しました。Kontsevich は、これを **特殊な行列模型（Kontsevich Airy 模型）** の評価により証明：

$$Z_{Kontsevich}(\Lambda) = \int dM \exp\left(\frac{i}{6}\text{Tr}\, M^3 + \frac{i}{2}\text{Tr}\, \Lambda M^2\right)$$

この積分は 3-valent ファットグラフの和に展開され、リーマン面の **trivalent ribbon graph による cell 分割**（Penner 1987 の decorated Teichmüller 理論の基礎）を介して、ψ-class の交点数を計算します。

### 5. Chekhov-Eynard-Orantin（CEO）位相的漸化式

CEO 漸化式は、**スペクトラル曲線**（複素代数曲線 + 微分形式 + 分岐構造のデータ）から、**全種数の不変量 ω_{g,n}**（n変数の有理形微分形式）を再帰的に計算する装置です。

```
ω_{g,n}(z₁, ..., zₙ) = Σ Res 〔 K(z, z₁) [ ω_{g-1,n+1}(z, σ(z), z₂, ..., zₙ)
                                    + Σ ω_{g₁,|I|+1} ω_{g₂,|J|+1} ] 〕
```

- 入力：スペクトラル曲線（複素曲線 Σ、二つの関数 x, y、二重点 B）
- 出力：すべての (g, n) に対する多重微分 ω_{g,n}、および自由エネルギー F_g

**驚異的な普遍性**：
- 行列模型 → ω_{g,n} は相関関数
- リーマン面のモジュライ空間 → ω_{g,n} は ψ-class 交点数
- トーリック Calabi-Yau のミラー曲線 → ω_{g,n} は Gromov-Witten 不変量（BKMP予想）
- Lambert 曲線 → ω_{g,n} は Hurwitz 数
- 双曲リーマン面 → ω_{g,n} は Mirzakhani の Weil-Petersson 体積
- 結び目補空間のA-多項式 → ω_{g,n} は colored Jones 多項式の漸近（量子曲線）

**一つの再帰関係式が、組合せ論・代数幾何・位相的弦・結び目理論・双曲幾何を統一的に生成する**。これが現代数理物理の中核ツールです。

---

## 第2部：W8（中心射影幾何）との接続点

木原様の W8 論文（Schläfli双対 + B₄ 同変性 の組合せ論的・位相的構成）と、ファットグラフ／CEO 漸化式が共有する構造的特徴：

1. **離散組合せ構造から連続幾何を生成**：
   - W8：Schläfli シンボル {p, q, r, ...} と B₄ 群作用から多胞体／曲率データを構成
   - CEO：リボングラフの組合せから ω_{g,n} を構成、リーマン面モジュライへ橋渡し

2. **双対性の組合せ的表現**：
   - W8：Schläfli 双対（多胞体の双対性）が B₄ 群構造に内蔵
   - ファットグラフ：planar dual graph により行列模型の M⁴ 模型と M³ 模型が双対関係（quadrangulation ↔ triangulation）

3. **群作用と分岐**：
   - W8：B₄（4次元 Coxeter 群）の同変構造
   - CEO：スペクトラル曲線の分岐点における対称性が漸化式の核 K(z, z₁) に直接表れる

4. **「種数展開」的階層構造**：
   - W8：次元別の射影構造（n次元 → n+1次元の埋め込み）の階層性
   - CEO：種数 g 別の階層的漸化

5. **共通の数学言語**：Riemann–Hurwitz の公式、cohomology of moduli spaces、Mirror Symmetry の B-model — どちらも触れる可能性。

**討論で持ち出せる具体的問い（メモ用）**：

- 中心射影幾何の組合せ的構成は、何らかのスペクトラル曲線として表現できるか？
- W8 の B₄ 同変性は、Eynard-Orantin の global symplectic invariant と何か関係するか？
- α 恒等式（W7） (π²/2)α² + 137α − 1 = 0 は、何らかの spectral curve の特殊化として現れるか？

---

## 第3部：収録文献リスト（カテゴリ別）

すべて PDF をローカル保存済み（一部 URL リンクのみ）。

### A. 講演者・藤博之氏の関連論文（直接の予習対象）

| # | 文献 | arXiv | ファイル |
|---|------|-------|---------|
| A1 | Andersen, Fuji, Manabe, Penner, Sulkowski (2016) "Partial chord diagrams and matrix models" | [1612.05840](https://arxiv.org/abs/1612.05840) | `A_藤博之氏論文/A1_*.pdf` |
| A2 | Andersen, Fuji, Manabe, Penner, Sulkowski (2016) "Enumeration of chord diagrams via topological recursion and quantum curve techniques" | [1612.05839](https://arxiv.org/abs/1612.05839) | `A_藤博之氏論文/A2_*.pdf` |
| A3 | Andersen, Chekhov, Penner, Reidys, Sulkowski (2012) "Topological recursion for chord diagrams, RNA complexes, and cells in moduli spaces" | [1205.0658](https://arxiv.org/abs/1205.0658) | `A_藤博之氏論文/A3_*.pdf` |
| A4 | Fuji, Iwaki, Manabe, Satake (2017→2019) "Reconstructing GKZ via topological recursion" | [1708.09365](https://arxiv.org/abs/1708.09365) | `A_藤博之氏論文/A4_*.pdf` |
| A5 | Fuji, Manabe (2023→2024) "Mirzakhani's Recursion and Masur-Veech Volumes via Topological Recursions" | [2303.14154](https://arxiv.org/abs/2303.14154) | `A_藤博之氏論文/A5_*.pdf` |
| A6 | Fuji, Gukov, Sulkowski (2012→2013) "Super-A-polynomial for knots and BPS states" | [1205.1515](https://arxiv.org/abs/1205.1515) | `A_藤博之氏論文/A6_*.pdf` |
| A7 | Fuji, Gukov, Sulkowski (2012) "Volume Conjecture: Refined and Categorified" | [1203.2182](https://arxiv.org/abs/1203.2182) | `A_藤博之氏論文/A7_*.pdf` |

### B. CEO 位相的漸化式 原典・基本文献

| # | 文献 | arXiv | ファイル |
|---|------|-------|---------|
| B1 | Eynard, Orantin (2007) "Invariants of algebraic curves and topological expansion" ★原典 | [math-ph/0702045](https://arxiv.org/abs/math-ph/0702045) | `B_CEO原典/B1_*.pdf` |
| B2 | Chekhov, Eynard (2006) "Hermitian matrix model free energy" ★原典 | [hep-th/0504116](https://arxiv.org/abs/hep-th/0504116) | `B_CEO原典/B2_*.pdf` |
| B3 | Eynard, Orantin (2009) "Algebraic methods in random matrices and enumerative geometry" ★レビュー | [0811.3531](https://arxiv.org/abs/0811.3531) | `B_CEO原典/B3_*.pdf` |
| B4 | Eynard (2014) "A short overview of the Topological recursion" ★ICM講演 | [1412.3286](https://arxiv.org/abs/1412.3286) | `B_CEO原典/B4_*.pdf` |
| B5 | Borot (2017) "Lecture notes on topological recursion and geometry" ★現代的入門 | [1705.09986](https://arxiv.org/abs/1705.09986) | `B_CEO原典/B5_*.pdf` |
| B6 | Eynard "Counting Surfaces" Birkhäuser (2016) ★教科書（書籍） | ISBN 978-3-7643-8796-9 | 未収録（書籍購入要） |

### C. 歴史的原典

| # | 文献 | アクセス | ファイル |
|---|------|---------|---------|
| C1 | 't Hooft (1974) "A planar diagram theory for strong interactions" Nucl. Phys. B72 | [著者ホスト PDF](https://webspace.science.uu.nl/~hooft101/gthpub/planar_diagram_theory.pdf) | **未収録**（サーバー応答せず／URL有効） |
| C2 | Brézin, Itzykson, Parisi, Zuber (1978) "Planar diagrams" Commun. Math. Phys. 59 | [自由公開 PDF](https://filippo-colomo.github.io/random_matrices/Brezin-Itzykson-Parisi-Zuber_78.pdf) | `C_歴史的原典/C2_*.pdf` |
| C3 | Kontsevich (1992) "Intersection theory on the moduli space of curves and the matrix Airy function" Commun. Math. Phys. 147 | [IHES 公開 PDF](https://www.ihes.fr/~maxim/TEXTS/Intersection%20theory%20and%20Airy%20function.pdf) | `C_歴史的原典/C3_*.pdf` |
| C4 | Penner (1987) "The decorated Teichmüller space of punctured surfaces" Commun. Math. Phys. 113 | [Project Euclid](https://projecteuclid.org/euclid.cmp/1104160216) | 未収録（要購読／open access あり） |

### D. モジュライ空間・ミラー対称性接続

| # | 文献 | arXiv | ファイル |
|---|------|-------|---------|
| D1 | Bouchard, Klemm, Mariño, Pasquetti (2007→2009) "Remodeling the B-model" ★BKMP予想 | [0709.1453](https://arxiv.org/abs/0709.1453) | `D_モジュライ・ミラー対称性/D1_*.pdf` |
| D3 | Eynard, Orantin (2007) "Weil-Petersson volume of moduli spaces, Mirzakhani's recursion and matrix models" | [0705.3600](https://arxiv.org/abs/0705.3600) | `D_モジュライ・ミラー対称性/D3_*.pdf` |
| D2 | Mirzakhani (2007) "Simple geodesics and Weil-Petersson volumes" Invent. Math. 167 | [AMS JAMS 関連](https://www.ams.org/journals/jams/2007-20-01/S0894-0347-06-00526-1/S0894-0347-06-00526-1.pdf) | 未収録（リンクのみ） |

### E. Hurwitz 数・結び目不変量への応用

| # | 文献 | arXiv | ファイル |
|---|------|-------|---------|
| E1 | Bouchard, Mariño (2008) "Hurwitz numbers, matrix models and enumerative geometry" | [0709.1458](https://arxiv.org/abs/0709.1458) | `E_Hurwitz・結び目応用/E1_*.pdf` |

結び目応用は A6 (super-A polynomial)、A7 (Volume Conjecture) も同時に該当。

---

## 第4部：ゼロから2週間で入門する推奨学習パス

| Day | 文献 | 所要 | 焦点 |
|-----|------|------|------|
| 1 | C1 't Hooft 1974（URL 経由で取得）または B4 § 1 のみ | 半日 | 1/N 展開とファットグラフの起源を15分で掴む |
| 2 | C2 BIPZ 1978 §1–§3 | 1日 | エルミート行列積分 = ファットグラフの和（手計算） |
| 3 | B4 Eynard 2014 全文 | 1日 | CEO 漸化式の全体像を最小コストで俯瞰 |
| 4–7 | B1 Eynard-Orantin 2007 §1–§4 | 4日 | スペクトラル曲線、二重微分 B、ω_{g,n}、再帰公式（原典で確認） |
| 8–10 | B5 Borot lecture notes 2017 | 3日 | Airy structure、Cohomological Field Theory、現代的視点 |
| 11 | C3 Kontsevich 1992 §1–§2 | 1日 | ファットグラフによる moduli space の cell 分割 |
| 12–13 | A1 + A2 Andersen-Fuji-Manabe-Penner-Sulkowski 2016 | 2日 | **講演者本人の論文。** Partial chord diagrams が行列模型 + CEO で解ける構造 |
| 14 | A5 Fuji-Manabe 2023 | 1日 | 藤氏の最新方向性（Masur-Veech, Mirzakhani 一般化、JT 重力） |

**最終余裕日（あれば）**：D1 BKMP、D3 Mirzakhani×CEO、E1 Hurwitz、A4 GKZ で応用面の幅を確認。

---

## 第5部：セミナー当日の対話に向けた準備

### 講演の予測される構成（概要欄より）

**基礎編**：
1. ファットグラフとその自己同型群
2. 行列模型と CEO 位相的漸化式

**応用編**（時間によって）：
3. 弦の場の理論の Hamiltonian に基づく定式化への応用
4. ファットグラフのデータサイエンスへの応用

### 「応用編」の意義

藤氏は応用編を **「完成した理論の説明のみとせず、今後の新たな研究の方向を探る機会」** と明記。これは、聴講者の独自視点・別分野からの問題提起を歓迎するサインです。

**4. ファットグラフのデータサイエンス応用**は、藤氏グループの A3 論文（RNA 複合体への応用）の延長線上にあり、組合せ的構造から創発する「構造的データ」をどう扱うかという視点。**離散組合せから連続構造を構成する** という共通言語が、W8 と直接的に結びつきうる領域です。

### 持ち込めるかもしれない問い（再掲・整理）

1. 「中心射影幾何（W8）では Schläfli 双対と B₄ 同変性を組合せ論的・位相的に構成しています。CEO の枠組みで、Schläfli 双対のような **次元横断的双対** をスペクトラル曲線の変形として表現できますか？」
2. 「α 恒等式 (π²/2)α² + 137α − 1 = 0 のような **自己整合方程式** が、ある spectral curve の特殊化として現れる例はありますか？（W7 の文脈）」
3. 「ファットグラフのデータサイエンス応用について、組合せ的構造から **創発する物理量** を見いだす方法論として、何が現在のボトルネックですか？」

---

## 第6部：補足リソース

- 講演者 researchmap：https://researchmap.jp/fuji360679/?lang=english
- 大阪城セミナー Indico ページ：https://indico.nitep.omu.ac.jp/event/221/
- NITEP 公式：https://www.nitep.osaka-cu.ac.jp/
- 連絡先：osaka_castle_seminar@gsuite.kobe-u.ac.jp

---

## 改訂履歴

- 2026-05-13：初版作成。19文献中16本のPDFをローカル保存（C1, C4, D2, B6 はURL/書籍のみ）。
