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
- **Version DOI (v0.5, 最新)**: 10.5281/zenodo.20981890
- **Zenodo deposit**: 20981890 / record https://zenodo.org/records/20981890（旧 v0.4: 20834424、v0.3: 20833097）
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

- **v0.5 (2026-06-28)**: §2.5「正規化とスケール不変性（注記）」と結論 (6) を追加。$\nu=1,\lambda=1$ は相対正規化（絶対スケールを想定しない）であり、奇数倍音は波形（相対的な局在の鋭さ $1/(N+1)$）を変えるのみで基本振動数 $\nu=1$・基本波長 $\lambda=1$ を変えないことを明示（$\gcd(1,3,\dots,N)=1$＝結論(2) の周期構造の言い換え、証明不要の観察）。これにより付録Aの $N\sim10^{38}$ は「絶対化したときの比」と相対化される。主張・数式・既存結果は不変。日英 md/tex/pdf 再生成。
- **v0.4 (2026-06-25)**: §2.4 の局在半幅 $\Delta_k$ の定義を、主ローブが最初に $k$ を切る点（誤）から、**サイドローブ込みで二度と $k$ を超えなくなる最後の交点**（正）へ修正。厳密逆算の特性値を主ローブの最初の根（$u\in(0,\pi)$）から $(\sin u/u)^2=k$ の最大根 $u_k^{\mathrm{out}}$ へ、近似式を $u_k\approx\pi/(1+\sqrt{k})$ から包絡線上界 $u_k^{\mathrm{out}}\lesssim1/\sqrt{k}$（$N\approx1/(\pi\sqrt{k}\,\Delta_k)-1$、安全側）へ変更。表と付録Aの数値を最後の交点基準に全面差し替え（$k=0.01$ で $\Delta$ が約2.95倍、$N=99$ で $0.908\%\to2.682\%$）。日英 md/tex/pdf 再生成。
- **v0.3 (2026-06-24 整理 → 2026-06-25 公開)**: §2.3 を新設し閉形式から $\sin u/u$ 近似への導出を $u=(N+1)\varphi$ 固定の極限として補強。用語を「中央主ピーク」に統一、孤立ピーク波の定義文を追加、局在幅を $1/(N+1)$ に整合、式番号を (2.1)–(2.16) の正規連番に振り直し。図1（孤立ピーク波）の本文埋め込みを追加。

### Zenn 記事

- [odd-harmonic-localization](https://zenn.dev/noriaki_kihara/articles/odd-harmonic-localization)

---

## 論文3：半波長奇数倍音孤立ピーク波の二コピー共通相対位相重ね合わせにおける波形不変性とコントラスト則の観察（v0.1）

Waveform Invariance and a Contrast Law for the Two-Copy Common-Relative-Phase Superposition of a Half-Wavelength Odd-Harmonic Isolated Peak Wave

- **Concept DOI**: 10.5281/zenodo.20923461（外部参照用・最新版へ自動転送）
- **Version DOI (v0.1, 最新)**: 10.5281/zenodo.20923462
- **Zenodo deposit**: 20923462 / record https://zenodo.org/records/20923462
- **公開日**: 2026-06-26
- **ライセンス**: CC BY 4.0
- **位置づけ**: 観察・整理論文。論文2の孤立ピーク波 $S_N$ の二つのコピーに、全倍音へ共通の相対位相 $\pm\alpha$ を与えて重ね合わせると、波積の三角恒等式から $\psi_\alpha=2\cos\alpha\cdot S_N$、$I_\alpha=4\cos^2\alpha\cdot I_N$ が厳密に成り立つ。これより (i) 規格化波形は相対位相に不変（単一コピーと一致）、(ii) 相対位相は二乗振幅を $\cos^2\alpha$ 倍するコントラスト因子としてのみ現れる、の二点を初等的なフーリエ和と三角恒等式の観察として記録する。物理的解釈は与えない。

### 収録ファイル

- `paper_relative_phase_contrast_ja_v0_1.md`（日本語）
- `paper_relative_phase_contrast_en_v0_1.md`（英語）
- `paper_relative_phase_contrast_ja_v0_1.tex` / `.pdf`（日本語）
- `paper_relative_phase_contrast_en_v0_1.tex` / `.pdf`（英語）
- 図2点（PNG、英語ラベル）
  - `coherent_self_interference_odd_modes.png`（奇数倍音と孤立ピーク波への合成）
  - `two_source_coherent_interference_corrected_v4_hires.png`（二コピー重ね合わせ、規格化波形が単一コピーと一致）

### 内容

- 二コピー重ね合わせ $\psi_\alpha=\sum[\cos((2m+1)\varphi-\alpha)+\cos((2m+1)\varphi+\alpha)]$ に和積の恒等式 $\cos(n\varphi-\alpha)+\cos(n\varphi+\alpha)=2\cos(n\varphi)\cos\alpha$ を適用し、$\psi_\alpha=2\cos\alpha\cdot S_N$ を厳密に導出。
- 波形不変性：規格化波形 $\widehat{I}_\alpha=\widehat{I}_N$ は相対位相に依存せず、形・両端 $\pm\pi/2$ の零・主ピーク近傍の局在幅が不変。
- コントラスト則：ピーク二乗振幅 $I_\alpha(0)=(N+1)^2\cos^2\alpha$、同相基準で比 $\cos^2\alpha$。$\alpha=0$ で最大、$\alpha=\pi/2$ で零。
- 具体例（$N=9$、$\alpha=15^\circ$）：$I_0(0)=100$、$I_\alpha(0)\approx 93.30$、$\cos^2 15^\circ\approx 0.9330$。
- 自己引用なし。外部引用は高木『解析概論』1件のみ。

### Zenn 記事

- [relative-phase-waveform-invariance](https://zenn.dev/noriaki_kihara/articles/relative-phase-waveform-invariance)

---

## 論文4：局在奇数倍音波の再生核性によるボルン分布の形の導出（v0.1）

Deriving the Form of the Born Distribution from the Reproducing-Kernel Property of a Localized Odd-Harmonic Wave: Reducing the Remaining Postulates to the Squaring Rule and Randomness

- **Concept DOI**: 10.5281/zenodo.20965526（外部参照用・最新版へ自動転送）
- **Version DOI (v0.1, 最新)**: 10.5281/zenodo.20965527
- **Zenodo deposit**: 20965527 / record https://zenodo.org/records/20965527
- **公開日**: 2026-06-27
- **ライセンス**: CC BY 4.0
- **自己参照**: 論文2「半波長位相区間における一定振幅奇数倍音和の孤立ピーク波とその局在性に関する観察」（Concept DOI 10.5281/zenodo.20833096）を `isSupplementTo` で参照
- **位置づけ**: 観察・整理論文。論文2の孤立ピーク波 $S_N$ を「観測の局在核」とみなし、(i) 観測＝局在核による位相差畳み込み＋二乗、(ii) 有限 $N$ 打ち切り、の二前提を置くと、$S_N$ の切り詰め再生核性により観測分布の**形**がベース波の $|\psi_{\rm base}|^2$（ボルン分布）に厳密一致することを示す。導出するのは「形」のみで、二乗則（振幅→確率）・確率解釈・ランダム性の起源は公準として残す。「ボルン則の導出」「測定問題の解決」は主張しない。

### 収録ファイル（claude/ サブフォルダ）

- `paper_born_from_localization_ja_v0_1.md` / `.tex` / `.pdf`（日本語、11 ページ）
- `paper_born_from_localization_en_v0_1.md` / `.tex` / `.pdf`（英語、13 ページ）
- 検証コード `born_from_localization.py`（sympy 記号＋数値＋複素拡張＋否定対照、機械精度再現）
- 図3点（PNG/SVG、英語ラベル）＋作図スクリプト
  - 図1 核の性質 `born_localization_kernel.{png,svg}`（`born_fig2.py`）
  - 図2 機構 `born_mechanism.{png,svg}`（`born_mechanism.py`）
  - 図3 再生の検証（4パネル）`born_from_localization.{png,svg}`（`born_fig.py`）
- 査読対話ノート `instruction_paper_born_claudecode.md` / `born_from_localization_review_iris.md` / `questions_to_claudeai_born_paper.md` / `answers_q1_q6_iris.md`

### 内容

- 区間 $[-\pi/2,\pi/2]$ 上で奇数倍音余弦が直交することから、$S_N=\sum_{m}\cos((2m+1)\varphi)$ は係数 1 の**切り詰め再生核**（ディリクレ核／RKHS の古典的事実）。新規性は恒等式ではなく〈論文2の物理的局在波＝再生核＝観測モデル〉という**写像**にある。
- 核の恒等式 $\int_{-\pi/2}^{\pi/2}\cos((2m+1)(\varphi_0-\varphi))\cos\varphi\,d\varphi=\frac{\pi}{2}\cos\varphi_0\,\delta_{m,0}$ から、$(S_N*\cos)(\varphi_0)=\frac{\pi}{2}\cos\varphi_0$、二乗で $\frac{\pi^2}{4}\cos^2\varphi_0$ が全奇 $N$ で厳密（端効果なし・境界含め点ごとに厳密）。
- 複素拡張：シフト核は複素奇数倍音基底の再生核で $(S_N*e^{ik\cdot})(\varphi_0)=\frac{\pi}{2}e^{ik\varphi_0}$。複素ベースで $|Z|^2=\mathrm{Re}^2+\mathrm{Im}^2$ の**真のモジュラス**（$\ne Z^2$）に到達。
- 否定対照：包絡読みは $1/\sin^2$ 発散・非可積分で確率になりえない。射影＋二乗が一意の操作化。
- 循環性への回答（§3.1）：歪まないのは等振幅のとき**だけ**（一般核は $c_m\to a_mc_m$ と歪む）。恒等写像は仮定でなく論文2の物理的局在波が再生核に一致するという事実の帰結。残す公準は二乗のべきのみ。
- Gleason／Zurek／決定理論的導出（二乗則そのものを導く）との対比で射程を明確化。

### 改訂履歴

- **v0.1 (2026-06-27)**: 初版公開。Claude Code 査読（ディリクレ核/RKHS の位置づけ明示、循環性への正面回答 §3.1、「フラクタル」→「局在核」、抄録の射程下方修正＋バンド条件明記）と機構図（図2）追加を反映。日英 md/tex/pdf＋検証コード＋図3点を Zenodo に収録（計13ファイル）。

### Zenn 記事

- [born-form-from-localized-kernel](https://zenn.dev/noriaki_kihara/articles/born-form-from-localized-kernel)

---

## 論文5：観測者–系ビートの等分布によるボルン統計の創発（予想論文、v0.2）

Emergence of Born Statistics from the Equidistribution of the Observer–System Beat: A Conjecture on the Localized-Kernel Model

- **Concept DOI**: 10.5281/zenodo.20967081（外部参照用・最新版へ自動転送）
- **Version DOI (v0.3, 最新)**: 10.5281/zenodo.20981910
- **Zenodo deposit**: 20981910 / record https://zenodo.org/records/20981910（旧 v0.2: 20967082）
- **公開日**: 2026-06-27
- **ライセンス**: CC BY 4.0
- **自己参照**: 論文4「局在奇数倍音波の再生核性によるボルン分布の形の導出」（Concept DOI 10.5281/zenodo.20965526）を `isSupplementTo` で参照
- **性格**: **予想（予言）論文**（証明論文でも観察論文でもない）。論文4が公準として残した「ランダム性の起源」に対する機構の候補を提示。二者査読（アイリス＋Claude Code）を経て v0.2。
- **位置づけ**: 走査変数 $\varphi_0$ は観測器の位相であり、観測器が有限振動数 $\psi$ をもつ物理対象ゆえ絶対参照が存在しない（主観空間の前提）。測定列は相対位相（ビート）$(\nu-\psi)t$ を走査し、無理比のとき等分布する（Weyl）。**(I) 条件付き証明済みの核**＝再生核性（論文4）＋ Weyl 等分布で位置分解強度プロファイルの形が $|\psi_{\rm base}|^2$ に一致（等分布の仮定 O1 の下で）。**(II) 予想の橋**＝強度→単発クリック確率（閾値検出に委ね導出しない）。「ボルン則の導出」「測定問題の解決」は主張しない。

### 収録ファイル（claude/ サブフォルダ）

- `paper_born_beat_conjecture_ja_v0_2.md` / `.tex` / `.pdf`（日本語、10 ページ）
- `paper_born_beat_conjecture_en_v0_2.md` / `.tex` / `.pdf`（英語、12 ページ）
- `paper_born_beat_conjecture_ja_v0_1.md`（初稿、研究過程の記録）
- 検証コード `born_beat_conjecture.py`（Weyl 等分布／再生核の実畳み込み両ケース／役割分離(D)／反証バッテリー(E)、機械精度）

### 内容

- **役割分担の三分離**：等分布（Weyl）＝一様測度の供給（偏りなし）／$|\psi|^2$ の形＝再生核の決定論的読み／強度→クリック＝橋（予想）。「Weyl だけで二乗の偏りが出る」という混同を排除（§5(D) で数値裏取り：等分布した $\varphi_0$ のヒストグラムは平坦 $1.0\times10^{-5}$）。
- **両ケース厳密導出**：局在入力 → 尖った出力／広がり入力 $\cos\varphi$ → 広がった出力 $\cos^2$ が同一再生核から出る（出力差＝入力差のみ）。実畳み込みで機械精度確認（$\le1.8\times10^{-15}$）、バンド条件（cos 全N／多モード N≥5）も明示。
- **反証バッテリー（§5 E）**：広がり入力が尖らない・局在入力が尖る・Weyl 単独は平坦・両ケースは同一恒等式、の四条件をすべてクリア（整合性の積極的証拠）。
- **先行研究の交点**：Khrennikov PCSFT ＋ 閾値検出（系譜A）／任意関数の方法 Poincaré–Feintzeig（系譜B）。固有の寄与＝等分布する高速変数のヘテロダイン同定＋再生核による形の厳密性。
- **反証可能な予言**：有限 $N$（離散格子）での $O(1/N)\sim\ell/R$ 補正＝ボルン則からの測定可能な逸脱。
- §5(D) の橋 MC は橋を**仮定した整合性デモ**であって導出・検証ではないことを明示し、核（証明済み）vs 橋（予想）の地位分離を維持。

### 改訂履歴

- **v0.3 (2026-06-28)**: (1) §7 の「勾配エネルギー $\sum\nu^2$ が有限化」というエネルギー発散を示唆する句を削除（反証可能予言の本体＝$O(1/N)$ 再生補正は不変）。(2) §2 末尾に「振動数の一致は自明」を明示（局在波は奇数倍音で波形を彫られても基本振動数 $\nu=1$ は不変＝スケール不変）し、点1（自明な振動数一致＋非自明な波形再現＝再生核）と点2（非自明な揺らぎ＝ビート）の役割分担を鮮明化。(3) 「ビートは一様サンプリングのみを与え、$|\psi|^2$ の偏りは再生核＋橋から来る」役割分担を §9 結論にも明示。主張・結論・数式は不変。日英 md/tex/pdf 再生成。
- **v0.2 (2026-06-27)**: 二者査読反映で初版公開。両ケース導出・§3(I) 役割分担の書き直し・反証バッテリー・バンド条件・「条件付き証明済み」統一・fd 自己引用削除・文献検証（arXiv:2409.16457「Lynnx」は単名表記と確認）・§5(D) 橋MCの地位明示。日英 md/tex/pdf＋検証コードを Zenodo に収録（計7ファイル）。

### Zenn 記事

- [born-beat-conjecture](https://zenn.dev/noriaki_kihara/articles/born-beat-conjecture)

---

## 研究記録

初期の思考実験・分割ログ・査読応答は、研究過程のタイムスタンプとして同フォルダに保存している。正式論文の主張範囲は `paper_square_quantity_readout_ja_v0_1.md` / `paper_square_quantity_readout_en_v0_1.md` を正とする。
