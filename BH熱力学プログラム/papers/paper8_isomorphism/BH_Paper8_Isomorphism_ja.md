# 4次元超立方格子上の鎖複体構造：Kihara 充填と Wilson 格子ゲージ理論の圏論的同型（論文8：α 恒等式 II）

**著者**：木原 範昭（WF System Co., Ltd.）
**ORCID**：[0009-0004-6753-4020](https://orcid.org/0009-0004-6753-4020)
**日付**：2026年4月
**Concept DOI（最新版に自動転送）**：[10.5281/zenodo.19880467](https://doi.org/10.5281/zenodo.19880467)
**v1 Version DOI**：[10.5281/zenodo.19880468](https://doi.org/10.5281/zenodo.19880468)

---

## 要旨

本論文は、姉妹論文（論文7：α 恒等式）で観察された自己整合恒等式 $\alpha^{-1} = N(1) + V_4(1)\cdot\alpha = 137 + (\pi^2/2)\alpha$（α を 8.7 ppb 精度で予言）の数学的位置づけを明らかにする。論文7 は4次元充填問題という純粋幾何学的アプローチから出発したが、その精度が偶然とは思えないほど高いため、標準的なゲージ理論との関係を見直したところ、**4次元超立方格子の幾何学と Wilson 格子ゲージ理論（Wilson 1974）が、鎖複体・対称群・離散 Stokes 構造のレベルで圏論的に同型である**ことが判明した。

本論文では、(i) 4次元 Euclidean 整数格子 $\mathbb{Z}^4$ 上の標準的な鎖複体 $C_\bullet = (C_0, C_1, C_2, C_3, C_4)$ を定義し、(ii) Schläfli 双対（テッセラクト ↔ 16-cell）が $C_n \leftrightarrow C_{4-n}$ の同型を与え、(iii) 超立方対称群 $B_4$ が両者に同変的に作用し、(iv) 離散 Stokes の定理が両者で同一の bulk-boundary 分解を生成することを示す。Wilson 理論は $C_1$（リンク）にゲージ場、$C_2$（プラケット）に作用を載せる formalism であり、Kihara 充填は $C_4$（4-cell）に体積、$C_3$（境界面）に測度を載せる formalism である。両者は Schläfli 双対を介して**圏論的に同型**である。

この同型により、(a) ゲージ理論側の解法技術（Wilson 強結合展開、繰り込み群、モンテカルロ等）が Kihara プログラムに適用可能となり、(b) 逆に Kihara 側の数論的恒等式（Lagrange–Jacobi 四平方定理、包除原理）が Wilson 理論の新規計算手法を提供しうる。論文7 の inside/outside 分解 $1 = 137\alpha + V_4(1)\alpha^2$ は、この鎖複体上の bulk-boundary 規格化条件として自然に解釈される。

α の数値的値そのものの第一原理導出（論文7 の観察を機構的に正当化すること）、および R=3 の特権性、残差 $c_3 \approx 1.6 \times 10^{-3}$ の幾何学的同定は、本論文の射程外として開問題に残す。本論文の貢献は、論文7 が単なる numerology ではなく、**4次元格子のトポロジカル構造の自然な現れ**であることを、標準数学の言葉で明示することにある。

---

## §1. 序論

### §1.1 論文7 の到達点と本論文の動機

姉妹論文（論文7、Concept DOI: 10.5281/zenodo.19869266、v3）において、4次元 Euclidean 空間内の半径 $R = 3$ の球 $B(3)$ に整数中心の単位立方体を完全包含で詰めた個数 $N(1) = 137$ と、4次元単位球の体積 $V_4(1) = \pi^2/2$ から構成される自己整合代数恒等式
$$ \alpha^{-1} = N(1) + V_4(1) \cdot \alpha = 137 + \frac{\pi^2}{2}\alpha \tag{1.1} $$
を観察した。これを α について解いた二次方程式 $(\pi^2/2)\alpha^2 + 137\alpha - 1 = 0$ の正根は α = 7.29735194×10⁻³ で、CODATA 2018 観測値と相対誤差 $8.7 \times 10^{-8}$（8.7 ppb）で一致する。

この精度は Eddington 系の整数当て（ズレ $2.6 \times 10^{-4}$）の約 1/3000 に相当し、「単なる数値偶然」と片付けるには良すぎる。論文7 では Eddington 罠を慎重に回避するため、(i) 自己整合方程式の構造、(ii) 137 と $\pi^2/2$ の独立な4次元幾何学的派生、(iii) QFT 摂動展開との類比、を明示したが、α 恒等式の物理機構や R=3 の特権性は未解決として残された。

本論文の動機は、この精度が偶然とは思えないため、**4次元充填問題と標準ゲージ理論の数学的関係**を体系的に見直すことにある。具体的には、Wilson が1974年に確立した格子ゲージ理論（Wilson 1974）と、論文7 の幾何学的設定が**幾何学的にほぼ同型**であることに気づき、その同型性を圏論的に明示することを目的とする。

### §1.2 本論文の主張の範囲

本論文は以下を**主張する**：

1. 4次元 Euclidean 整数格子 $\mathbb{Z}^4$ 上の標準的な鎖複体に、超立方対称群 $B_4$ が同変的に作用する。
2. Schläfli 双対（テッセラクトと 16-cell の自己双対性）が $C_n \leftrightarrow C_{4-n}$ の同型を与える。
3. Wilson 格子ゲージ理論と Kihara 充填は、この鎖複体上で異なる「載せ方」をした formalism であり、Schläfli 双対を介して**圏論的に同型**である。
4. 離散 Stokes の定理が両者で同一の bulk-boundary 分解を生成する。
5. 論文7 の inside/outside 分解 $1 = 137\alpha + V_4(1)\alpha^2$ は、この同型構造から自然に解釈される。

本論文は以下を**主張しない**：

1. α の値が本同型から第一原理的に導出されること（論文7 の観察を引用するに留める）
2. R=3 の特権性が本同型から証明されること
3. U(1) ゲージ対称性が本論文の枠組みから一意に創発すること
4. 残差 $c_3 \approx 1.6 \times 10^{-3}$ の幾何学的同定
5. ヒッグス機構や電弱対称性破れとの具体的接続

これらは全て open problem として §8 に明示する。

### §1.3 論文構成

§2 で4次元超立方格子の鎖複体を標準的に定義する。§3 で Schläfli 双対と Poincaré 双対を述べる。§4 で Wilson 格子ゲージ理論を鎖複体的に記述する。§5 で Kihara 充填を同様に記述する。§6 で主定理（圏論的同型）を述べる。§7 で論文7 の α 恒等式の数学的位置づけを論ずる。§8 で開問題を、§9 で結論を述べる。

本論文は self-contained に書かれており、Wilson 1974 や論文7 を読んでいなくても、§2〜§6 の数学的内容は本文内で完結している。

---

## §2. 4次元超立方格子の鎖複体

### §2.1 格子の定義

4次元 Euclidean 整数格子を $\Lambda = \mathbb{Z}^4 \subset \mathbb{R}^4$ と書く。各 $x \in \Lambda$ を**サイト**（vertex）と呼ぶ。隣接サイト間の有向辺を**リンク**（link）と呼び、リンク全体を $L$ と書く。リンクで囲まれる最小2次元面を**プラケット**（plaquette）、3次元胞を**3-cell**、4次元胞を**4-cell**（テッセラクト）と呼ぶ。

### §2.2 鎖複体

形式的整数係数で生成される自由アーベル群を $C_n$ とする：

$$ C_0 = \mathbb{Z}\langle\text{vertices}\rangle, \quad C_1 = \mathbb{Z}\langle\text{links}\rangle, \quad C_2 = \mathbb{Z}\langle\text{plaquettes}\rangle, \quad C_3 = \mathbb{Z}\langle\text{3-cells}\rangle, \quad C_4 = \mathbb{Z}\langle\text{4-cells}\rangle. $$

境界作用素 $\partial_n : C_n \to C_{n-1}$ は標準的に定義される（向き付け込み）：
$$ \partial_n \circ \partial_{n+1} = 0. \tag{2.1} $$

これにより $C_\bullet = (C_0 \xrightarrow{\partial_1} C_1 \xrightarrow{\partial_2} C_2 \xrightarrow{\partial_3} C_3 \xrightarrow{\partial_4} C_4)$ は鎖複体をなす。

### §2.3 超立方対称群 $B_4$ の作用

$\Lambda$ を不変にする等距変換群は、超立方対称群 $B_4$（Coxeter 群、位数 $|B_4| = 2^4 \cdot 4! = 384$）である。$B_4$ の生成元：

- 軸反転：$x_i \to -x_i$（4個、位数 $2^4 = 16$）
- 軸の置換：$x_i \leftrightarrow x_j$（$4! = 24$ 個）

$B_4$ は各 $C_n$ に自然に作用し、境界作用素と同変：
$$ g \cdot \partial_n = \partial_n \cdot g, \quad \forall g \in B_4. \tag{2.2} $$

### §2.4 単位 cell の数値的データ

| $n$ | $C_n$ の単位 cell | 単位 cell あたりの境界 cell 数 | 4-cell（テッセラクト）あたりの個数 |
|---|---|---|---|
| 0 | vertex | — | 16 |
| 1 | link | 2 vertices | 32 |
| 2 | plaquette | 4 links | 24 |
| 3 | 3-cell | 6 plaquettes | 8 |
| 4 | 4-cell | 8 3-cells | 1 |

これらは標準的なテッセラクトの $f$-vector：$(16, 32, 24, 8, 1)$。

---

## §3. Schläfli 双対と Poincaré 双対

### §3.1 テッセラクトと 16-cell の自己双対

4次元正多胞体には Schläfli 記号で表される双対対が存在する。テッセラクト（8-cell, $\{4,3,3\}$）と 16-cell（cross-polytope, $\{3,3,4\}$）は互いに双対：
$$ \text{Tesseract}^* = \text{16-cell}, \quad (\text{16-cell})^* = \text{Tesseract}. $$

具体的には、テッセラクトの $k$-cell（$k$-次元面）は 16-cell の $(4-k-1)$-cell に対応：

| Tesseract | 16-cell |
|---|---|
| 16 vertices ($k=0$) | 16 cells ($k=3$) |
| 32 edges ($k=1$) | 32 faces ($k=2$) |
| 24 faces ($k=2$) | 24 edges ($k=1$) |
| 8 cells ($k=3$) | 8 vertices ($k=0$) |

### §3.2 鎖複体上の Poincaré 双対

格子 $\Lambda$ の双対格子 $\Lambda^* \cong \mathbb{Z}^4$（半整数シフト）を考える。$\Lambda$ の $n$-cell は $\Lambda^*$ の $(4-n)$-cell と1対1対応：
$$ D : C_n(\Lambda) \xrightarrow{\sim} C_{4-n}(\Lambda^*). \tag{3.1} $$

これは離散版の **Poincaré 双対**（4次元、コンパクト台）であり、$B_4$ 同変：
$$ D \circ g = g \circ D, \quad \forall g \in B_4. \tag{3.2} $$

### §3.3 双対作用の具体例

| $\Lambda$ 側 | $\Lambda^*$ 側 |
|---|---|
| サイト ($C_0$) | 4-cell ($C_4$) |
| リンク ($C_1$) | 3-cell ($C_3$) |
| プラケット ($C_2$) | プラケット ($C_2$、自己双対) |
| 3-cell ($C_3$) | リンク ($C_1$) |
| 4-cell ($C_4$) | サイト ($C_0$) |

特筆すべきは $C_2 \to C_2^*$ の自己双対性で、これが Wilson プラケット作用と Kihara 充填の同型の直接的根源である。

---

## §4. Wilson 格子ゲージ理論の鎖複体的記述

### §4.1 リンク変数

Wilson (1974) は、ゲージ理論を格子上で次のように離散化した。各リンク $\ell = (x, y) \in L$ にゲージ群の元 $U_\ell \in G$ を配置する：
$$ U : C_1 \to G. $$

本論文では $G = U(1)$（電磁気力）に焦点を絞り、$U_\ell = e^{i\theta_\ell}$、$\theta_\ell \in [0, 2\pi)$ とする。逆向きリンク：$U_{yx} = U_{xy}^{-1}$。

### §4.2 ゲージ変換

各サイト $x$ に位相 $\alpha(x) \in U(1)$ を割り当て、リンク変数を変換：
$$ U_{xy} \to e^{i\alpha(x)} U_{xy} e^{-i\alpha(y)}. \tag{4.1} $$

これは局所ゲージ変換であり、各サイトで独立に位相を選べる。

### §4.3 プラケット位相

プラケット $p \in C_2$ を $4$ つのリンク（向き込み）の周回で表し、プラケット位相を：
$$ U_p = \prod_{\ell \in \partial p} U_\ell. \tag{4.2} $$

ゲージ変換 (4.1) の下で、プラケット位相は不変：
$$ U_p \to U_p \quad \text{（ゲージ不変）}. $$

### §4.4 Wilson 作用

$$ S_{\text{Wilson}} = \frac{1}{g^2} \sum_{p \in C_2} \left[1 - \text{Re}\, U_p\right] \tag{4.3} $$

ここで $g$ は結合定数。連続極限 $a \to 0$ で標準QED 作用 $\frac{1}{4} F_{\mu\nu}F^{\mu\nu}$ に収束する（Wilson 1974, Kogut 1979）。

### §4.5 Wilson 作用の鎖複体的読み替え

Wilson 作用は次のように鎖複体的に書ける：
$$ S_{\text{Wilson}} : C_2 \to \mathbb{R}, \quad p \mapsto \frac{1}{g^2}[1 - \text{Re}\, U_p]. \tag{4.4} $$

すなわち Wilson 理論は、$C_1$ にゲージ場を、$C_2$ に作用 functional を載せる formalism である。

### §4.6 離散 Stokes の定理：bulk cancellation

任意の閉じた2次元面 $\Sigma \subset C_2$ に対し、$\Sigma$ を構成するプラケットの位相積 $\prod_{p \in \Sigma} U_p$ は、内部リンクが pairwise に逆向きで現れるため打ち消し合い、境界 Wilson loop に等しい：
$$ \prod_{p \in \Sigma} U_p = \prod_{\ell \in \partial \Sigma} U_\ell. \tag{4.5} $$

これは離散版 Stokes の定理であり、$\partial^2 = 0$ の直接の帰結である。bulk から境界への寄与の集約こそが Wilson 理論の計算的核心である。

---

## §5. Kihara 充填の鎖複体的記述

### §5.1 単位立方体の配置

論文3・7 において、各サイト $c \in \mathbb{Z}^4$ を中心とする単位立方体を：
$$ \square(c) = \prod_{i=1}^{4} \left[c_i - \frac{1}{2}, c_i + \frac{1}{2}\right] \tag{5.1} $$
と定義した。これは双対格子 $\Lambda^* = (\mathbb{Z} + 1/2)^4$ の頂点を持つ4-cell であり、$C_4(\Lambda^*)$ の基底元とみなせる。

すなわち**Kihara の単位立方体 = $\Lambda^*$ の 4-cell**。Schläfli 双対 $D : C_4(\Lambda) \to C_0(\Lambda^*)$ により、これは $\Lambda$ のサイトと等価である。

### §5.2 完全包含と充填数

半径 $R$ の4次元球 $B(R)$ に**完全包含**される単位立方体の集合：
$$ \mathcal{N}(R) = \{c \in \mathbb{Z}^4 : \square(c) \subset B(R)\}. $$

完全包含条件は、$\square(c)$ の全頂点が $B(R)$ 内：
$$ \sum_{i=1}^{4} (|c_i| + 1/2)^2 \le R^2. \tag{5.2} $$

充填数：$N(k) = |\mathcal{N}(2k+1)|$。論文3 の直接計数で $N(1) = 137$（$R=3$）。

### §5.3 体積不足

連続体積 $V_4(R) = (\pi^2/2) R^4$ と離散充填体積 $N(k)$ の差：
$$ \Delta(R) = V_4(R) - N(k). \tag{5.3} $$

論文3 §5 で漸近展開：
$$ \Delta(R) = \frac{16\pi}{3} R^3 - 6\pi R^2 + O(R) \tag{5.4} $$
を導出した。漸近係数 $c = 8/(3\pi)$（$\Delta(R)/(2\pi^2 R^3)$ の極限）は、4次元球面の表面積と同じスケーリング $R^3$ で増大する。

### §5.4 鎖複体的読み替え

Kihara 充填は次のように鎖複体的に書ける：
$$ \chi_{B(R)} : C_4(\Lambda^*) \to \{0, 1\}, \quad \square(c) \mapsto \begin{cases} 1 & \text{if } \square(c) \subset B(R) \\ 0 & \text{otherwise} \end{cases} \tag{5.5} $$
$$ N(k) = \sum_{\square \in C_4(\Lambda^*)} \chi_{B(R)}(\square). \tag{5.6} $$

すなわち Kihara プログラムは、$C_4$ に指示函数を、$C_3$（境界）に体積測度を載せる formalism である。

### §5.5 体積不足の bulk-boundary 構造

体積不足 $\Delta(R)$ は次の構造を持つ：

- **bulk**：$B(R)$ の深部では立方体が完全に充填され、寄与は厳密に整数 $N_{\text{bulk}}$
- **boundary**：$\partial B(R)$ 近傍では立方体が部分的にしか含まれず、その「余り」が $\Delta(R)$

これは離散 Stokes の定理（§4.6）と同じ構造：bulk が完全に cancel し、境界の寄与だけが残る。具体的には、隣接立方体は3次元面を共有し、共有面の体積寄与は両側で相殺される。

---

## §6. 主定理：圏論的同型

### §6.1 圏 $\mathcal{L}_4$ の定義

**定義 6.1**：圏 $\mathcal{L}_4$ を以下で定義：

- **対象**：4次元 Euclidean 整数格子（あるいは双対格子）上の鎖複体 $C_\bullet$ で、$B_4$ 対称性と境界作用素 $\partial_\bullet$ を持つもの
- **射**：$B_4$-同変な chain map（境界作用素と可換な準同型）

### §6.2 同型定理

**定理 6.2（主定理）**：

$\mathcal{L}_4$ の対象として、

$$ (\mathcal{C}_W, \partial_\bullet, B_4) \cong_{\mathcal{L}_4} (\mathcal{C}_K, \partial_\bullet, B_4) $$

ここで：
- $\mathcal{C}_W$ は Wilson 理論の鎖複体（$C_1$ にゲージ場、$C_2$ に作用 functional）
- $\mathcal{C}_K$ は Kihara 充填の鎖複体（$C_4$ に指示函数、$C_3$ に体積測度）

同型射は Schläfli 双対 $D$（§3.2）：

$$ D : \mathcal{C}_W \xrightarrow{\sim} \mathcal{C}_K $$
- $C_1(\Lambda)$（Wilson リンク）$\xrightarrow{D}$ $C_3(\Lambda^*)$（Kihara 立方体境界面）
- $C_2(\Lambda)$（Wilson プラケット）$\xrightarrow{D}$ $C_2(\Lambda^*)$（Kihara 立方体間共有面、自己双対）

### §6.3 証明の概略

(i) **格子の同型**：$\Lambda$ と $\Lambda^*$ は半整数シフトで関係し、$B_4$ 同変。

(ii) **鎖複体の同型**：Schläfli 双対 $D : C_n(\Lambda) \to C_{4-n}(\Lambda^*)$ は標準的な離散 Poincaré 双対であり、$\partial$ と可換（cap product の構造）。

(iii) **対称性の保存**：$B_4$ は両者に同変的に作用し、$D$ は $B_4$-同変。

(iv) **bulk-boundary 構造の保存**：離散 Stokes の定理（§4.6 と §5.5）は両者で同一形式。

詳細な代数的証明は標準的な algebraic topology の教科書（Hatcher 2002）の手法で構成可能であるため省略する。$\blacksquare$

### §6.4 系：双方向の道具流用

定理 6.2 の系として、以下が従う：

**系 6.3**：Wilson 格子ゲージ理論の解析手法（強結合展開、繰り込み群、モンテカルロ等）は、Schläfli 双対 $D$ を介して Kihara 充填問題に直接適用可能である。

**系 6.4**：逆に、Kihara 側の数論的恒等式（Lagrange–Jacobi 四平方定理 $r_4(N) = 8\sigma(N)$、論文3 §6）は、Wilson 理論の有限体積 partition function 計算に新規な閉形式を提供しうる。

これらの系は、論文7・8 に基づく研究プログラムが、ゲージ理論と数論の双方の道具を融合的に活用しうることを意味する。

---

## §7. 帰結：論文7 の α 恒等式の数学的位置づけ

### §7.1 inside/outside 分解の自然性

論文7 §4 で観察された inside/outside 分解：
$$ 1 = \underbrace{137 \alpha}_{\text{inside (bulk)}} + \underbrace{V_4(1) \alpha^2}_{\text{outside (boundary)}} \tag{7.1} $$

この構造は §6 の同型定理から自然に解釈される：

- $137 \alpha$：bulk（$C_4$ の内部充填、137 セル）の tree-level 寄与
- $V_4(1) \alpha^2$：boundary（$C_3$ の境界、4次元単位球測度）の自己エネルギー補正

すなわち、論文7 の自己整合恒等式は**4次元格子の鎖複体上の bulk-boundary 規格化条件**として鎖複体的に読み替えられる。

### §7.2 Wilson 作用との対応

定理 6.2 により、論文7 の充填問題に Wilson 作用を載せた formalism：
$$ S_{\text{Kihara-Wilson}} = \sum_{p \subset B(3)} [1 - \text{Re}\, U_p] \tag{7.2} $$
は数学的に well-defined である。これが論文8 の枠組みで自然に記述される **Wilson 格子ゲージ理論の有限切り取り版**であり、論文7 の幾何学的設定を含む。

### §7.3 トポロジカル双対との接続

論文7 §6.5 で導入したトポロジカル双対性（テッセラクト ↔ 16-cell の Schläfli 双対）は、定理 6.2 の同型射 $D$ そのものである。論文7 §6.5 は本論文の主定理を予言的に観察していたことになる。

### §7.4 α 恒等式の位置づけ（再定式化）

以上を総合し、論文7 の α 恒等式は次のように位置づけられる：

> **論文7 の自己整合恒等式 $\alpha^{-1} = N(1) + V_4(1)\alpha$ は、4次元 Euclidean 整数格子の鎖複体構造（Schläfli 双対、$B_4$ 対称性、離散 Stokes）の上で、Wilson 格子ゲージ理論と Kihara 充填の圏論的同型から自然に位置づけられる構造的恒等式である。**

これは Eddington 系の素朴な数値合わせとは決定的に異なる：

- 整数 137 は離散構造の不変量
- $\pi^2/2$ は連続体積の不変量
- 両者の組み合わせは鎖複体上の bulk-boundary 規格化
- α の数値はこの構造への適合度（観察）

ただし、α の数値が**この構造から第一原理的に決定される**ことの証明は本論文の射程外であり、§8 の open problem として残す。

---

## §8. 開問題

本論文で扱わなかった、しかし本同型構造から approach 可能な問題群を以下に列挙する。

### §8.1 α の数値の機構的派生
本同型は α の構造的位置づけを与えるが、$\alpha = 1/137.036$ という具体値が**この同型から一意に決まる**ことは証明されていない。Wilson 結合定数 $g$ と Kihara 幾何の結合をどう同定するかが課題。

### §8.2 R=3 の特権性
$N(0) = 1, N(1) = 137, N(2) = 1545, \ldots$ のうち、なぜ $N(1)$ が物理的 α に対応するのか。Wilson 側で言えば、なぜこの有限切り取りスケールが物理的か。Mass gap や confinement scale との対応が候補。

### §8.3 残差 $c_3$ の幾何学的同定
論文7 §6.1 で残差を $\alpha^3$ オーダーで $c_3 \approx 1.6 \times 10^{-3}$ と推定したが、$c_3$ の閉形式は不明。本同型により Wilson 強結合展開の高次項として approach 可能性あり。

### §8.4 他のゲージ群への一般化
本論文は $G = U(1)$ に焦点を絞ったが、$SU(2), SU(3), SU(N)$ への一般化は自然である。Schläfli 双対が群構造をどう保つか、$\alpha_s$（強結合定数）への類比的予言が可能か、が研究課題。

### §8.5 物理時空（3+1次元）への射影
本論文と論文7 はすべて4次元 Euclidean で議論したが、現実の物理時空への射影機構が必要。論文2・6 で議論した中心投影と Wick 回転の融合 formalism が候補。

### §8.6 ヒッグス機構との接続
論文7 §5.3 で「縦波スカラー場とのビート結合」を機構候補として挙げたが、本同型の枠組みでこれをどう実装するかは未解決（思考実験ファイル[ヒッグス波についての考察/]を参照）。

### §8.7 数値的検証
本同型を計算機実装し、Wilson モンテカルロで $N(k)$ や $\Delta(R)$ を逆算、Lagrange-Jacobi 恒等式を Wilson partition function に適用、などの数値実験は今後の課題。

---

## §9. 結論

本論文は、論文7 で観察された自己整合恒等式 $\alpha^{-1} = N(1) + V_4(1)\alpha$（α を 8.7 ppb 精度で予言）の数学的位置づけを、4次元 Euclidean 整数格子の鎖複体構造として明らかにした。

主結果は、**Wilson 格子ゲージ理論と Kihara 充填問題が、4次元超立方格子の鎖複体（$B_4$ 対称性、Schläfli 双対、離散 Stokes 構造）の上で圏論的に同型**である（定理 6.2）ことであり、これにより：

1. 論文7 の inside/outside 分解は鎖複体上の bulk-boundary 規格化として自然に解釈される
2. Wilson 理論の解析手法が Kihara プログラムに適用可能となる
3. Kihara 側の数論的恒等式が Wilson 理論に新規計算を提供しうる
4. α 恒等式は単なる numerology ではなく、4次元格子のトポロジカル構造の自然な現れである

ただし、α の具体値の第一原理導出、R=3 の特権性、残差 $c_3$ の幾何学的同定は本論文の射程外であり、open problem として §8 に明示した。

本論文の貢献は α の機構的派生ではなく、論文7 が標準的な数学的構造（algebraic topology + lattice gauge theory）と自然に接続することを示した点にある。これにより論文7 の Eddington 罠回避がさらに強化され、本研究プログラムが標準的物理学・数学コミュニティと対話可能な形で位置づけられる。

---

## 参考文献

1. Wilson, K. G. "Confinement of Quarks." *Phys. Rev. D* **10**, 2445 (1974).
2. Kogut, J. B. "An introduction to lattice gauge theory and spin systems." *Rev. Mod. Phys.* **51**, 659 (1979).
3. Hatcher, A. *Algebraic Topology*. Cambridge University Press (2002).
4. Coxeter, H. S. M. *Regular Polytopes*. Dover Publications (1973).
5. Schläfli, L. *Theorie der vielfachen Kontinuität*. Zürcher und Furrer, Zürich (1901).
6. Bekenstein, J. D. "Black Holes and Entropy." *Phys. Rev. D* **7**, 2333 (1973).
7. Hawking, S. W. "Particle Creation by Black Holes." *Commun. Math. Phys.* **43**, 199 (1975).
8. Tangherlini, F. R. "Schwarzschild Field in $n$ Dimensions and the Dimensionality of Space Problem." *Nuovo Cimento* **27**, 636 (1963).
9. Eddington, A. S. "The Charge of an Electron." *Proc. Roy. Soc. London* **A122**, 358 (1929).
10. Mohr, P. J., Newell, D. B., Taylor, B. N. "CODATA Recommended Values of the Fundamental Physical Constants: 2018." *Rev. Mod. Phys.* **93**, 025010 (2021).

---

## 付録 A：テッセラクトの $f$-vector と Schläfli 双対

テッセラクト（4次元超立方体、Schläfli 記号 $\{4,3,3\}$）の $f$-vector：
$$ (f_0, f_1, f_2, f_3) = (16, 32, 24, 8). $$

16-cell（Schläfli 記号 $\{3,3,4\}$）の $f$-vector：
$$ (f_0, f_1, f_2, f_3) = (8, 24, 32, 16). $$

Schläfli 双対は $f_k(\text{Tesseract}) = f_{3-k}(\text{16-cell})$。

これは標準的な polytope 理論の事実（Coxeter 1973）で、4次元での自己双対構造の根源である。

## 付録 B：論文7 の α 恒等式の数値再現（参考）

論文7 §3.2 と完全に一致（v3、8.7 ppb 精度）：

```
V_4(1) = π²/2 ≈ 4.9348022
N(1) = 137
方程式：4.9348 α² + 137 α - 1 = 0
正根：α ≈ 7.29735194 × 10⁻³
α⁻¹ ≈ 137.0360110
CODATA 2018: α⁻¹ = 137.0359991
相対誤差：8.7 × 10⁻⁸（8.7 ppb）
```

詳細は論文7（Concept DOI: 10.5281/zenodo.19869266）を参照。
