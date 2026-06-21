# RELEASE NOTES — 平方数を基本量とした場合の検討

本フォルダは、正の平方写像

$$
X_i=x_i^2
$$

を「平方量読出し」として扱い、平方和制約・線形単体化・初等運動形式・曲率補正候補の出現次元を整理する論文および研究記録を収録する。

物理法則の証明、観測事実の主張、既存物理の置換は行わない。中心的な位置づけは、自己参照論文0の測地単位セル歪みの結果と、平方量読出しを並置して、長さだけで閉じる式と二次元面積セルを検討すべき式を分ける整理原理である。

---

## 論文1：平方量読出しによる線形単体化と曲率補正候補の次元別整理（v0.1-r1）

Linear Simplexification by Square-Quantity Readout and a Dimension-wise Organization of Curvature-Correction Candidates

- **Concept DOI**: 10.5281/zenodo.20785539（外部参照用・最新版へ自動転送）
- **Version DOI (v0.1-r1, 最新)**: 10.5281/zenodo.20785540
- **Zenodo deposit**: 20785540 / record https://zenodo.org/records/20785540
- **公開日**: 2026-06-21
- **ライセンス**: CC BY 4.0
- **自己参照**: 論文0「正曲率定曲率空間における測地的単位セルの歪み」（Concept DOI 10.5281/zenodo.20680269 / Version DOI 10.5281/zenodo.20684135）
- **位置づけ**: 観察・整理論文。平方写像そのものの新規性ではなく、平方量読出しと論文0の次元別歪み構造を並置し、曲率補正候補がどの幾何次元で問題になるかを整理する。

### 収録ファイル（公開予定）

- `paper_square_quantity_readout_ja_v0_1.md`（日本語）
- `paper_square_quantity_readout_en_v0_1.md`（英語）
- `paper_square_quantity_readout_ja_v0_1.tex` / `.pdf`（日本語）
- `paper_square_quantity_readout_en_v0_1.tex` / `.pdf`（英語）
- 図5点（PNG/SVG、英語ラベル）
  - `figures/fig01_2d_square_map.{png,svg}`
  - `figures/fig02_3d_square_map.{png,svg}`
  - `figures/fig03_motion_readouts.{png,svg}`
  - `figures/fig04_quadratic_readings.{png,svg}`
  - `figures/fig05_area_coefficient_ks.{png,svg}`
- 作図スクリプト `figures/make_square_quantity_figures.py`
- 図仕様 `figures/figure_manifest.md`

### 内容

- 正の平方写像 $X_i=x_i^2$ により、平方和制約 $\sum_i x_i^2=E$ が平方量側で線形単体制約 $\sum_i X_i=E$ として読めることを定式化。
- ただし、これは計量の等長平坦化ではなく、制約式の代数的単純化であることを明記。
- 平方量側の一次式 $X=AT$ と二次式 $X=\frac12 BT^2$ が、正の平方根読出しで等速運動形・等加速度運動形と同型になることを示す。
- 一次元完全弾性衝突の式では、質量・運動量・エネルギーを外から置いても、一次元速度だけで閉じる限り面積補正が出ないことを確認。
- 遠心型モデルを、接線方向と半径方向が同時に現れるため二次元面積セルを検討する自然な入口として位置づける。
- 自己参照論文0の一辺1の測地正方形面積係数

  $$
  k_s(R)
  =
  R^2
  \left[
  4\arccos\!\left(-\tan^2\frac{1}{2R}\right)-2\pi
  \right]
  $$

  を用いるが、これを遠心型力学式の一意な導出とは主張しない。
- 平方写像の非等長性と、論文0の測地正方形面積超過を同一視しないことを明記。

### 改訂履歴

- **v0.1-r1 (2026-06-21)**: ChatGPT/Claude Code の査読コメントと著者助言を反映。主張範囲を「発見論文」から「観察・整理論文」へ整理し、平方写像の非等長性から $k_s(R)$ を因果的に導く読みを明示的に否定。遠心型モデルへの $k_s(R)$ 挿入を「最小試験モデル」として限定。英語図5点を作成。

### Zenn 記事

- [square-quantity-readout-simplexification](https://zenn.dev/noriaki_kihara/articles/square-quantity-readout-simplexification)

---

## 研究記録

初期の思考実験・分割ログ・査読応答は、研究過程のタイムスタンプとして同フォルダに保存している。正式論文の主張範囲は `paper_square_quantity_readout_ja_v0_1.md` / `paper_square_quantity_readout_en_v0_1.md` を正とする。
