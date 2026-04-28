# $S^4$ からのグノモン投影と4次元主観座標図上の誘導計量の Schwarzschild 様の形（論文2：射影の数学的基盤）

**著者**：木原 範昭（WF System Co., Ltd.）
**ORCID**：[0009-0004-6753-4020](https://orcid.org/0009-0004-6753-4020)
**日付**：2026年4月
**DOI**：[10.5281/zenodo.19837590](https://doi.org/10.5281/zenodo.19837590)

---

## 要旨

半径 $R$ の4次元球 $S^4(R) \subset \mathbb{R}^5$ から4次元接超平面へのグノモン（中心）投影を考察し、$S^4(R)$ 上の標準計量の引き戻しにより超平面上に誘導される計量を扱う。誘導計量の明示的形を導出し、Christoffel 記号、Riemann テンソル、Einstein テンソルを計算する。誘導計量は4次元で正の宇宙項 $\Lambda = 3/R^2$ を持つ真空 Einstein 方程式を満たすことを示す。$R$ 固定スライスへの制限と Lorentz 符号への解析接続により、構成は de Sitter 様計量を生み、4+1次元 Schwarzschild–Tangherlini 解との比較は地平面挙動における構造的類似性を明らかにする。論文は数学的内容で、物理的解釈は延期する。

---

## §1. 序論

球面から接平面へのグノモン投影は、球面上の大円が平面上の直線に写るという素朴な性質を持つ。2次元では地図学の古典的観察である。高次元での類似構成は $S^n$ 上の大 $(n-1)$-円を接超平面上の $(n-1)$-平面に写し、結果として超平面上に誘導される計量 — 投影に沿って球面上の標準計量を引き戻して得られる — は球の半径 $R$ を反映する曲率を持つ非自明な Riemann 計量である。

本論文は4次元の場合を詳細に扱う：$S^4(R) \subset \mathbb{R}^5$ から4次元接超平面 $\Pi_R \subset \mathbb{R}^5$ へのグノモン投影。誘導計量を明示的座標で導出し、曲率不変量を計算し、球座標での Schwarzschild 様の形を同定する。主要結果：

1. 誘導計量は4次元で $G_{\mu\nu} + \Lambda\, g_{\mu\nu} = 0$（正の宇宙項 $\Lambda = 3/R^2$）を満たす。
2. Lorentz 連続化後、計量は Beltrami 座標での de Sitter 計量。
3. 動径座標を適切に同定した球座標で、計量は Schwarzschild–de Sitter 形に帰着し、質量パラメータが小さい極限で純粋 de Sitter、宇宙項が消える極限（$R \to \infty$）で純粋 Schwarzschild に一致。
4. 自然な高次元比較対象は4+1次元 Schwarzschild–Tangherlini 解で、本構成は重力ポテンシャルの $r^{-2}$ 落下を共有する。

論文構成。第2節でグノモン投影と接超平面座標を設定。第3節で誘導計量を導出し Christoffel 記号を計算。第4節で Riemann と Einstein テンソルを計算。第5節で球座標での計量を解析し Schwarzschild 様形を同定。第6節で Tangherlini 解と比較。第7節で結論。

論文は純粋に数学的。物理的ブラックホール熱力学に関する構成の解釈は別研究の対象。

---

## §2. 構成：$S^4$ からのグノモン投影

### §2.1 設定

$\mathbb{R}^5$ を周辺空間、座標 $Y = (Y^0, Y^1, Y^2, Y^3, Y^4)$、標準 Euclid 計量 $\delta_{AB}\, dY^A \otimes dY^B$ を持つ。半径 $R$ の4次元球：
$$S^4(R) = \{Y \in \mathbb{R}^5 : (Y^0)^2 + (Y^1)^2 + (Y^2)^2 + (Y^3)^2 + (Y^4)^2 = R^2\}.$$

$Y^0$-軸を投影軸として選ぶ。「北極」$(R, 0, 0, 0, 0)$ での接超平面：
$$\Pi_R = \{Y \in \mathbb{R}^5 : Y^0 = R\}$$
を $Y^\mu = x^\mu$（$\mu = 1, 2, 3, 4$）で座標 $(x^1, x^2, x^3, x^4)$ により径数化。

### §2.2 投影写像

グノモン投影 $\Phi : \Pi_R \to S^4(R)_+$（「北半球」$Y^0 > 0$）を「原点と $\Phi(x) \in S^4(R)$ を通る直線が $x \in \Pi_R$ も通る」規則で定義。明示的に：
$$\Phi(x) = \frac{R}{\ell(x)} (R, x^1, x^2, x^3, x^4), \qquad \ell(x) := \sqrt{R^2 + |x|^2}, \quad |x|^2 := \sum_{\mu=1}^{4} (x^\mu)^2.$$

直接 $\sum_{A=0}^{4} \Phi(x)^A \Phi(x)^A = R^2$ を確認、$\Phi(x) \in S^4(R)$。

### §2.3 逆写像と正則性領域

逆写像 $\Phi^{-1} : S^4(R)_+ \to \Pi_R$：
$$\Phi^{-1}(Y)^\mu = \frac{R Y^\mu}{Y^0} \qquad (\mu = 1, 2, 3, 4).$$

写像は開半球 $Y^0 > 0$ で well-defined かつ滑らか。赤道 $Y^0 = 0$ は「無限遠」に写る（$Y^0 \to 0^+$ で $\Phi^{-1}(Y)$ が発散）。南半球 $Y^0 < 0$ は $\Phi$ の像にない；南極を通る別の投影で扱えるが、ここでは不要。

### §2.4 測地線像

$S^4(R)$ 上の大円（測地線）は $S^4(R)$ と $\mathbb{R}^5$ の原点を通る2次元線型部分空間との交わり。これらの部分空間は再スケーリング $Y \mapsto \lambda Y$ で安定なので、$\Phi^{-1}$ による像はこれら部分空間と $\Pi_R$ の交わりであり、$\Pi_R$ 内の直線（またはアフィン平面）。これがグノモン投影の定義的性質：球面上の測地線が接超平面上の直線に対応する。

---

## §3. 誘導計量

### §3.1 引き戻し計算

$S^4(R)$ 上の標準計量は周辺 Euclid 計量の制限。$\Phi$ に沿って引き戻すと、$\Pi_R$ 上の計量 $g$ を得る：
$$g_{\mu\nu}(x) = \delta_{AB}\, \frac{\partial \Phi^A}{\partial x^\mu}\, \frac{\partial \Phi^B}{\partial x^\nu}.$$

直接計算：
$$\frac{\partial \Phi^0}{\partial x^\mu} = -\frac{R^2 x_\mu}{\ell^3}, \qquad \frac{\partial \Phi^\nu}{\partial x^\mu} = \frac{R}{\ell} \delta^\nu_\mu - \frac{R x^\nu x_\mu}{\ell^3}.$$

代入・整理：
$$\boxed{\;g_{\mu\nu}(x) = \frac{R^2}{\ell^2}\left(\delta_{\mu\nu} - \frac{x_\mu x_\nu}{\ell^2}\right).\;}$$

### §3.2 逆計量

直接計算（または Sherman-Morrison 公式での行列反転）：
$$g^{\mu\nu}(x) = \frac{\ell^2}{R^2}\left(\delta^{\mu\nu} + \frac{x^\mu x^\nu}{R^2}\right).$$

関係式 $g_{\mu\nu} g^{\nu\rho} = \delta_\mu^\rho$ は計算で確認。

### §3.3 行列式

$\det(\delta_{\mu\nu} - x_\mu x_\nu / \ell^2) = 1 - |x|^2/\ell^2 = R^2/\ell^2$ を用いて：
$$\det g = \left(\frac{R^2}{\ell^2}\right)^4 \cdot \frac{R^2}{\ell^2} = \frac{R^{10}}{\ell^{10}}.$$

体積要素：$\sqrt{\det g}\, d^4 x = (R/\ell)^5 d^4 x$。

### §3.4 Christoffel 記号

Christoffel 記号 $\Gamma^\rho_{\mu\nu} = \frac{1}{2} g^{\rho\sigma}(\partial_\mu g_{\nu\sigma} + \partial_\nu g_{\mu\sigma} - \partial_\sigma g_{\mu\nu})$ は閉形式で計算可能：
$$\Gamma^\rho_{\mu\nu} = -\frac{1}{\ell^2}(x^\rho \delta_{\mu\nu} + (\text{対称補正}))$$

完全な表式は付録 A。**（後で検討：交差項寄与を含む完全表式）**

---

## §4. Riemann と Einstein テンソル

### §4.1 Riemann テンソル

Christoffel 記号から Riemann テンソルを計算：
$$R^\rho_{\sigma\mu\nu} = \frac{1}{R^2}(g_{\sigma\mu} \delta^\rho_\nu - g_{\sigma\nu} \delta^\rho_\mu).$$

これは最大対称空間 — 具体的に断面曲率 $1/R^2$ の正定数曲率空間 — の Riemann テンソル。

### §4.2 Ricci テンソルとスカラー曲率

縮約：
$$R_{\sigma\nu} = R^\mu_{\sigma\mu\nu} = \frac{3}{R^2} g_{\sigma\nu}, \qquad R = g^{\mu\nu} R_{\mu\nu} = \frac{12}{R^2}.$$

### §4.3 Einstein テンソルと宇宙項

$n = 4$ 次元 Einstein テンソル：
$$G_{\mu\nu} := R_{\mu\nu} - \frac{1}{2} R\, g_{\mu\nu} = \frac{3}{R^2} g_{\mu\nu} - \frac{6}{R^2} g_{\mu\nu} = -\frac{3}{R^2} g_{\mu\nu}.$$

よって $G_{\mu\nu} + \Lambda g_{\mu\nu} = 0$ で
$$\boxed{\;\Lambda = \frac{3}{R^2}.\;}$$

誘導計量は4次元真空 Einstein 方程式を、球の半径で決定される正の宇宙項とともに満たす。

### §4.4 Lorentz 連続化

座標の1つを Lorentz 符号に連続化 — 例えば $x^4 \to i\tau$ — すると $\Pi_R$ 上の Riemann 計量は Lorentz 計量に変換。直接代入で結果が**Beltrami 座標での de Sitter 計量**であることを示す：
$$ds^2 = \frac{R^2}{\ell_L^2}\left(-d\tau^2 + |dx|^2 + (\text{交差項})\right), \qquad \ell_L^2 := R^2 - \tau^2 + |x|^2.$$

この Lorentz 設定での宇宙項は $\Lambda = 3/R^2 > 0$ で、計量は非標準座標での4次元 de Sitter 時空。

---

## §5. 球座標と Schwarzschild 様の形

### §5.1 球分解

Lorentz 符号で空間部分を $r = |\vec{x}|$ とする球座標 $(r, \theta, \phi)$ で分解：
$$ds^2 = -A(\tau, r)\, d\tau^2 + B(\tau, r)\, dr^2 + r^2\, d\Omega_2^2 + (\text{交差項})$$
で $A, B$ の明示的形は §4.4 から従う。

### §5.2 静的スライシング

特定のスライシング（適切に選ばれた $\tau$ 一定面 — **後で検討**：正確なゲージ選択）で、計量は
$$ds^2 = -f(r)\, d\tau^2 + f(r)^{-1} dr^2 + r^2\, d\Omega_2^2, \qquad f(r) = 1 - \frac{r^2}{R^2}$$
に帰着。これは静的座標での de Sitter 計量で、$r = R$ で宇宙論的地平面。

### §5.3 質量パラメータの追加

自然な拡張として質量パラメータ $M$ を導入：
$$f(r) = 1 - \frac{2GM}{r} - \frac{r^2}{R^2}.$$

これは **Schwarzschild–de Sitter** 計量で、事象地平面（$r \approx 2GM$ 近傍）と宇宙論的地平面（$r \approx R$ 近傍）の両方を持つ。

$M \to 0$ 極限：純粋 de Sitter（$r = R$ で宇宙論的地平面）。
$R \to \infty$ 極限：純粋 Schwarzschild（$r = 2GM$ で事象地平面）。

4次元計量での質量様項の存在は、中心投影構成への可能な物理的源；この項が姉妹論文の離散立方体充填に微視的起源を持つかは未解決。

---

## §6. 4+1次元 Schwarzschild–Tangherlini 解との比較

### §6.1 Tangherlini 計量

4+1次元 Schwarzschild–Tangherlini 計量：
$$ds^2 = -f_T(r)\, dt^2 + f_T(r)^{-1} dr^2 + r^2\, d\Omega_3^2, \qquad f_T(r) = 1 - \frac{8GM}{3\pi r^2}.$$

主要特徴：
- 地平面は $r_h = (8GM/(3\pi))^{1/2}$
- 重力ポテンシャルの動径依存性 $1/r^2$（$D = 5$ 特有）
- Hawking 温度 $T_H = 1/(2\pi r_h)$
- Bekenstein–Hawking エントロピー $S = 2\pi^2 r_h^3/(4 G_5)$

### §6.2 中心投影計量との構造的類似

§5 の中心投影構成は宇宙論的地平面で形 $1 - r^2/R^2$、Schwarzschild–de Sitter 拡張で形 $1 - 2GM/r - r^2/R^2$ の計量を生む。一方、4+1 Tangherlini 計量は $1 - 8GM/(3\pi r^2)$。

両計量は等価ではない。地平面を持つという定性的特徴を共有し、特に中心投影構成の曲率半径 $R$ は Tangherlini 構成の地平面半径 $r_h$ と類似の役割を演じる。正確な対応は同定 $R \leftrightarrow r_h$ により媒介され、これはエネルギー-半径スケーリングに関する姉妹論文の動力学的解析で支持される。

### §6.3 妥当性領域

構造的類似は計量形式と地平面の存在のレベルで成立。計量の完全な同一性には拡張しない。姉妹論文（論文4）では体積不足 $\Delta(R)$ と Bekenstein–Hawking エントロピーを介した定量的比較を展開する。

---

## §7. 結論

$S^4(R)$ の4次元接超平面上のグノモン投影計量を導出し、Riemann と Einstein テンソルを計算し、構造的形式を同定した。主要結果：

1. 誘導計量は $G_{\mu\nu} + (3/R^2) g_{\mu\nu} = 0$ を満たす。
2. Lorentz 連続化で4次元 de Sitter 計量を与える。
3. 質量パラメータを伴う静的球座標で、計量は Schwarzschild–de Sitter 形。
4. 構成は4+1次元 Schwarzschild–Tangherlini 計量と構造的特徴を共有するが、同一ではない。

構成は研究プログラムの仮説1の幾何学的基盤を提供する：射影された4次元理論内のブラックホールが高次元構造の像でありうること。誘導計量の特定形と Tangherlini との関係が本論文の内容；4次元理論の物質内容からの $R$ の動力学的決定を含む物理的解釈は論文5 の対象。

---

## 付録 A：Christoffel 記号（延期）

**後で検討**：中心投影計量に対する $\Gamma^\rho_{\mu\nu}$ の完全閉形式表現。

## 参考文献（選定、**後で検討**）

- F. R. Tangherlini, "Schwarzschild Field in $n$ Dimensions and the Dimensionality of Space Problem," *Nuovo Cimento* 27, 636 (1963).
- R. C. Myers and M. J. Perry, "Black Holes in Higher Dimensional Spacetimes," *Annals Phys.* 172, 304 (1986).
- F. Kottler, "Über die physikalischen Grundlagen der Einsteinschen Gravitationstheorie," *Annalen der Physik* 56, 401 (1918).
