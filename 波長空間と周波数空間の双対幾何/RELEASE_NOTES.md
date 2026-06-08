# RELEASE NOTES — 波長空間と周波数空間の双対幾何（Dual Geometry of Wavelength Space and Frequency Space）

本シリーズは、波長空間と周波数空間の間に逆数双対条件を置いたときに現れる幾何学的・位相的構造を観察するモデルである。物理理論の導出・修正を主張するものではなく、$\nu_n$ を時間周波数・エネルギー・運動量等の物理量と同一視しない。

---

## 論文1：波長空間と周波数空間の双対幾何（v0.3）

Dual Geometry of Wavelength Space and Frequency Space: A Geometric and Topological Observational Model of Reciprocal Conditions, Logarithmic Representation, and Uncertainty-Weighted Counting

- **Concept DOI**: 10.5281/zenodo.20588036
- **Version DOI (v0.3)**: 10.5281/zenodo.20588037
- **Zenodo deposit**: 20588037
- **公開日**: 2026-06-08
- **ライセンス**: CC BY 4.0
- **関連**: 論文2（10.5281/zenodo.20588039）を `isSupplementedBy` で関連付け

### 収録ファイル（10点）
- `paper1_wavelength_frequency_dual_geometry_v0_3.md`（日本語）
- `paper1_wavelength_frequency_dual_geometry_v0_3_en.md`（英語）
- `paper1_wavelength_frequency_dual_geometry_v0_3.tex` / `.pdf`（日本語）
- `paper1_wavelength_frequency_dual_geometry_v0_3_en.tex` / `.pdf`（英語）
- 図（日本語）: `λ_ν_双対条件の模式図.png` / `λ_ν_双対の対数表現.png`
- 図（英語）: `figure1_lambda_nu_dual_constraint_EN.png` / `figure2_log_representation_lambda_nu_duality_EN.png`

### 内容
- 逆数双対条件 $\lambda_n=1/\nu_n$ ＋ 各空間の二乗和一定条件。
- 1次元では解がほぼ一点に固定、5成分1制約では4自由度が残る。
- 対数変換 $q=\log\lambda$, $p=\log\nu$ で逆数双対は符号反転対称 $p=-q$ に単純化。
- 周波数／波長空間を4次元格子に選ぶと二乗和条件は単位セル数え上げに変換。$\delta=0$ の完全内接数 $N_0(R)$、$\delta>0$ の重み付き数 $N_\delta(R)$ を定義（$\delta$・重み関数の値は導出せず後続課題）。
- §9 関連研究と位置づけ／参考文献（Gabor 1946、Shannon 1948、Aliev–Henk 2023、Hirschhorn 1987）を数理的背景として最小限に付与。

---

## 論文2：4次元格子における単位セル完全内接数の半径スイープ（v0.1）

Radius Sweep of Fully-Inscribed Unit-Cell Counts on a 4-Dimensional Lattice: An Enumeration Table from R = 0.5 to 10.0 with a Reproducible Formulation

- **Concept DOI**: 10.5281/zenodo.20588038
- **Version DOI (v0.1)**: 10.5281/zenodo.20588039
- **Zenodo deposit**: 20588039
- **公開日**: 2026-06-08
- **ライセンス**: CC BY 4.0
- **関連**: 論文1（10.5281/zenodo.20588037）を `isSupplementTo` で関連付け

### 収録ファイル（7点）
- `paper2_4d_lattice_cell_count_radius_sweep_v0_1.md`（日本語）
- `paper2_4d_lattice_cell_count_radius_sweep_v0_1_en.md`（英語）
- `paper2_4d_lattice_cell_count_radius_sweep_v0_1.tex` / `.pdf`（日本語）
- `paper2_4d_lattice_cell_count_radius_sweep_v0_1_en.tex` / `.pdf`（英語）
- `paper2_radius_sweep_R_0_5_to_10_step_0_5.csv`（数え上げ表データ）

### 内容
- 4次元整数格子上の一辺1の単位セルが半径 $R$ の4次元超球体に完全内接する個数 $N_0(R)$ を $R=0.5$〜$10.0$（0.5刻み）で数え上げ。
- 完全内接条件 $\sum_{i=1}^{4}(|k_i|+\tfrac12)^2\le R^2$。再現用擬似コードとCSVを同梱。
- $N_0(1)=1$, $N_0(2)=9$, $N_0(3)=137$。$R=3$ の137は殻分解 $1+8+24+40+64$ により得られる純粋な数え上げ結果（物理定数との対応は主張しない）。

---

## 論文3：閉じた4自由度構造と4次元格子数え上げの対応（v0.2、論文1の追補）

Closed Four-Degree-of-Freedom Structure and Its Correspondence with 4-Dimensional Lattice Counting: A Geometric Organization from the 5-Component Sum-of-Squares Constraint to the Unit-Cell Counting Region

- **Concept DOI**: 10.5281/zenodo.20589261
- **Version DOI (v0.2)**: 10.5281/zenodo.20589262
- **Zenodo deposit**: 20589262
- **公開日**: 2026-06-08
- **ライセンス**: CC BY 4.0
- **関連**: 論文1（10.5281/zenodo.20588037）を `isSupplementTo` で関連付け

### 収録ファイル（6点）
- `paper3_closed_4d_structure_and_lattice_counting_supplement_v0_2.md`（日本語）
- `paper3_closed_4d_structure_and_lattice_counting_supplement_v0_2_en.md`（英語）
- `paper3_closed_4d_structure_and_lattice_counting_supplement_v0_2.tex` / `.pdf`（日本語）
- `paper3_closed_4d_structure_and_lattice_counting_supplement_v0_2_en.tex` / `.pdf`（英語）

### 内容
- 論文1の5成分二乗和制約 $\sum_{n=1}^{5}x_n^2=R^2$ は5次元中の4次元超球面 $S^4_R$（4自由度）を定義し、論文2の数え上げは4次元球体 $B^4_R$ の内部で行う、という対応を物理解釈なしに整理。
- 球面投影 $\Pi_R(y)=R\,y/\|y\|$ は、制約を満たす点（$\|\lambda\|=\Lambda$）に対して恒等写像 $\lambda'=\lambda$, $\nu'=\nu$。投影は値の変換ではなく、制約点を半径一定の閉じた4自由度構造上の点として読む幾何学的記述。
- 用語は「球面投影」に統一（v0.2）。$S^4_R$ 上の測地的分割・物理対応は対象外（今後の課題）。

---

## 公開先リンク

- **Zenn（技術解説・数式表示あり）**: https://zenn.dev/noriaki_kihara/articles/wavelength-frequency-dual-geometry
- **note（日本語版、図1・2付き）**: https://note.com/kiharanoriaki/n/n08aeb3c4e8ae ＜2026-06-08 公開＞
- **note（英語版、図1・2付き）**: https://note.com/kiharanoriaki/n/nf2b3e4392ea1 ＜2026-06-08 公開＞
- **Facebook（日本語版、カバー図＋コメントに note URL 誘導）**: ＜2026-06-08 公開＞（投稿原稿 `facebook_article_ja.md`）
- **Facebook（英語版、カバー図＋コメントに note URL 誘導）**: ＜2026-06-08 公開＞（投稿原稿 `facebook_article_en.md`）
- **X / Twitter（@NoriakiKihara、日本語版・英語版、note URL＋ハッシュタグ）**: ＜2026-06-08 公開＞（投稿原稿 `x_post.md`、note OGPカバー図が自動表示）

---

## 変更履歴

### 2026-06-08
- 新シリーズ「波長空間と周波数空間の双対幾何」を開始。
- 論文1 v0.3 を初版公開（Concept DOI 10.5281/zenodo.20588036、Version DOI 10.5281/zenodo.20588037）。
- 論文2 v0.1 を初版公開（Concept DOI 10.5281/zenodo.20588038、Version DOI 10.5281/zenodo.20588039）。
- 数式区切りを `$`/`$$` に統一、標準ヘッダ（著者・ORCID・版・DOI・ライセンス）を付与。
- 論文1に参考文献4件（Gabor / Shannon / Aliev–Henk / Hirschhorn）を追加。
- 日本語・英語の md / tex / pdf と図を各 Zenodo レコードへアップロード。
- Zenn 紹介記事を公開（`articles/wavelength-frequency-dual-geometry.md`）。
- note 用 日英記事（図1・2付き・ハッシュタグ）を作成。**note 日本語版**（https://note.com/kiharanoriaki/n/n08aeb3c4e8ae ）・**note 英語版**（https://note.com/kiharanoriaki/n/nf2b3e4392ea1 ）を公開。
- Facebook 用 日英記事（タイトル→コメント誘導→概要→ハッシュタグ構成、コメントに note URL）を作成し、**日本語版・英語版とも公開**（カバー図付き）。
- X / Twitter 用 日英ポスト（note URL＋ハッシュタグ、280字制限内）を作成し、**@NoriakiKihara で日本語版・英語版とも公開**（note OGPカバー図が自動表示）。
- **論文3 v0.2（論文1の追補）を初版公開**（Concept DOI 10.5281/zenodo.20589261、Version DOI 10.5281/zenodo.20589262、論文1へ `isSupplementTo`）。数式を `$`/`$$` に統一、用語を「球面投影」に統一、日英 md/tex/pdf をアップロード。Zenn 記事に追補節を追加（FB・X・note は当面据え置き、note は後日追補予定）。
