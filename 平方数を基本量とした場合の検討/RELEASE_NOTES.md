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

### note 記事

- 日本語: [平方数で読むと、曲がった平方和が直線になる ── 曲率歪みと平方量読出しの整理](https://note.com/kiharanoriaki/n/n879eac8f6cb2)
- English: [Reading by Squares: When a Curved Sum of Squares Becomes a Straight Line](https://note.com/kiharanoriaki/n/n1f41301fdf49)

### SNS 投稿

- Facebook 日本語投稿: 個人タイムライン（2026-06-22、本文2行目「詳細はコメント参照」、コメントに日本語 note URL、図添付）
- Facebook English post: personal timeline (2026-06-22, second line "See the comment for details", English note URL in comment, figure attached)
- X 日本語投稿: [@NoriakiKihara](https://x.com/NoriakiKihara)（2026-06-22、日本語 note URL、図添付）
- X English post: [@NoriakiKihara](https://x.com/NoriakiKihara) (2026-06-22, English note URL, figure attached)

---

## 論文2：半波長位相区間における一定振幅奇数倍音和の孤立ピーク波とその局在性に関する観察（v0.4）

An Observation on the Isolated Peak Wave of a Constant-Amplitude Odd-Harmonic Sum on a Half-Wavelength Phase Interval and Its Localization

- **Concept DOI**: 10.5281/zenodo.20833096（外部参照用・最新版へ自動転送）
- **Version DOI (v0.4, 最新)**: 10.5281/zenodo.20834424
- **Zenodo deposit**: 20834424 / record https://zenodo.org/records/20834424（旧 v0.3: 20833097）
- **公開日**: 2026-06-25（v0.3 公開 → 同日 v0.4 改訂公開）
- **ライセンス**: CC BY 4.0
- **位置づけ**: 観察・整理論文。半波長区間 $[-\pi/2,\pi/2]$ 上で一定振幅の奇数倍音を余弦で重ね合わせると、中央に主ピークをもち両端で零となる「孤立ピーク波」が形成されること、その局在幅が $1/(N+1)$ で縮むこと、指定した局在幅から必要な最高倍音次数を逆算する式を、初等的なフーリエ和の性質として整理する。物理的解釈は与えない。

### 収録ファイル

- `paper_odd_harmonic_localization_ja_v0_1.md`（日本語）
- `paper_odd_harmonic_localization_en_v0_1.md`（英語）
- `paper_odd_harmonic_localization_ja_v0_1.tex` / `.pdf`（日本語）
- `paper_odd_harmonic_localization_en_v0_1.tex` / `.pdf`（英語）
- 図2点（PNG/SVG、英語ラベル）
  - `figures/fig01_odd_harmonic_localization.{png,svg}`（孤立ピーク波）
  - `figures/fig02_odd_harmonic_scaling.{png,svg}`（局在幅の $1/(N+1)$ スケーリング）
- 作図スクリプト `figures/make_odd_harmonic_figure.py` / `figures/make_odd_harmonic_scaling_figure.py`

### 内容

- 一定振幅奇数倍音和 $S_N(\varphi)=\sum_{m=0}^{(N-1)/2}\cos((2m+1)\varphi)$ を定義し、等差数列の余弦和の公式から閉形式 $S_N(\varphi)=\sin((N+1)\varphi)/(2\sin\varphi)$ を独立に導出（ディリクレ核を前提としない）。
- 拡大変数 $u=(N+1)\varphi$ を固定して中央主ピーク近傍を見ると、ピーク値で規格化した形が $N$ によらない普遍関数 $\sin u/u$ に近づくことを段階的に導出（分母 $2\sin\varphi$ が $2\sin(u/(N+1))$ へ書き換わり、小角展開とピーク値規格化で $N$ 依存が相殺する筋を明示）。
- 規格化二乗振幅 $\widehat{I}_N\approx(\sin u/u)^2$、中央主ピークの横幅が $1/(N+1)$ で縮小。
- $k$-局在半幅 $\Delta_k$ を「中心から離れて二度と $k$ を超えなくなる**最後の交点**（サイドローブ包絡線が $k$ に落ちきる外縁）」と定義し、特性値を $(\sin u/u)^2=k$ の**最大根** $u_k^{\mathrm{out}}$（$k=0.01$ で $8.4232$、$k=0.001$ で $30.151$）として必要な最高奇数倍音次数を逆算する数値解代入法と、包絡線上界 $u_k^{\mathrm{out}}\lesssim1/\sqrt{k}$ による安全側の近似閉形式 $N\approx1/(\pi\sqrt{k}\,\Delta_k)-1$ を提示。
- 付録Aはスケール比の算術例にすぎず、物理的実在・物理過程を主張しない（ボーア半径は任意の小スケールの一例）。

### 改訂履歴

- **v0.4 (2026-06-25)**: §2.4 の局在半幅 $\Delta_k$ の定義を、主ローブが最初に $k$ を切る点（誤）から、**サイドローブ込みで二度と $k$ を超えなくなる最後の交点**（正）へ修正。厳密逆算の特性値を主ローブの最初の根（$u\in(0,\pi)$）から $(\sin u/u)^2=k$ の最大根 $u_k^{\mathrm{out}}$ へ、近似式を $u_k\approx\pi/(1+\sqrt{k})$ から包絡線上界 $u_k^{\mathrm{out}}\lesssim1/\sqrt{k}$（$N\approx1/(\pi\sqrt{k}\,\Delta_k)-1$、安全側）へ変更。表と付録Aの数値を最後の交点基準に全面差し替え（$k=0.01$ で $\Delta$ が約2.95倍、$N=99$ で $0.908\%\to2.682\%$）。日英 md/tex/pdf 再生成。
- **v0.3 (2026-06-24 整理 → 2026-06-25 公開)**: §2.3 を新設し閉形式から $\sin u/u$ 近似への導出を $u=(N+1)\varphi$ 固定の極限として補強。用語を「中央主ピーク」に統一、孤立ピーク波の定義文を追加、局在幅を $1/(N+1)$ に整合、式番号を (2.1)–(2.16) の正規連番に振り直し。図1（孤立ピーク波）の本文埋め込みを追加。

### Zenn 記事

- [odd-harmonic-localization](https://zenn.dev/noriaki_kihara/articles/odd-harmonic-localization)

---

## 研究記録

初期の思考実験・分割ログ・査読応答は、研究過程のタイムスタンプとして同フォルダに保存している。正式論文の主張範囲は `paper_square_quantity_readout_ja_v0_1.md` / `paper_square_quantity_readout_en_v0_1.md` を正とする。
