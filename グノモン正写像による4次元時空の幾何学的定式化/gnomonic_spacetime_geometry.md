# グノモン正写像による4次元時空の幾何学的定式化

**提出者：** 理学研究科 宇宙物理学講座 修士1年  
**提出日：** 2026年4月（最終版・第3回査読対応済）  
**種別：** 中間ゼミ報告（幾何学的考察篇）  
**指導教員：** 宇宙物理学講座 教授、数学講座 准教授

---

## 概要

本報告は，4次元ユークリッド空間をグノモン正写像（中心射影）によって曲率半径 $R$ の5次元超球体表面へ写像したとき，写像後の座標系において生じる誘導計量の変形を純粋に幾何学的観点から考察するものである．特に，(i) 写像後の誘導計量テンソル $g_{\mu\nu}$ の厳密な導出（第5成分の寄与を含む），(ii) アインシュタインテンソルと定数 $\Lambda = 3/R^2$ との幾何学的関係，を段階的に示す．物理的解釈は次年度以降の課題とし，本報告では一切の宇宙論的仮定を排除する．

---

## 1. はじめに

グノモン写像（gnomonic projection）は，球の中心から接平面への中心射影として定義され，その最も重要な性質として「球面上の大円（測地線）が接平面上の直線に写る」ことが知られている [1]．

本報告では，この写像を4次元から5次元へ拡張し，写像後の座標系に生じる計量構造を調べる．この問題設定は純粋に幾何学的であり，物理的解釈（重力，時空，相対性理論等）は本報告の範囲外とする．

---

## 2. 準備：グノモン写像の定義と基本性質

### 2.1 一般次元のグノモン写像

**定義 2.1**（$n$ 次元グノモン写像）[1, 2]  
$n$ 次元ユークリッド空間において，球の中心から距離 $R$（曲率半径 $R$ に対応）の接点を持つ接超平面 $\Pi_R$（接点を原点とする座標 $(x^1, \ldots, x^n) \in \mathbb{R}^n$ を持つ $n$ 次元ユークリッド平面）を考える．$l = \sqrt{R^2 + \sum_{i=1}^n (x^i)^2}$ とおくとき，

$$\Phi: \Pi_R \to S^n(R) \subset \mathbb{R}^{n+1}, \quad (x^1, \ldots, x^n) \mapsto \left(\frac{R x^1}{l}, \ldots, \frac{R x^n}{l},\ \frac{R^2}{l}\right) \tag{2.1}$$

を**グノモン正写像**と呼ぶ．ここで $S^n(R)$ は $\mathbb{R}^{n+1}$ に埋め込まれた半径 $R$ の $n$ 次元球面

$$S^n(R) = \left\{(Y^1, \ldots, Y^{n+1}) \in \mathbb{R}^{n+1} \;\middle|\; \sum_{A=1}^{n+1} (Y^A)^2 = R^2 \right\}$$

である．

**命題 2.1**（大円保存）[1]  
グノモン写像のもとで，接平面上の任意の直線は球面上の大円（測地線）に写る．

---

## 3. 計量テンソルの導出

### 3.1 写像の設定

本節では $n=4$ として，接超平面 $\Pi_R$（座標 $(x^1, x^2, x^3, x^4)$）から半径 $R$ の5次元超球面 $S^4(R) \subset \mathbb{R}^5$ へのグノモン写像を考える．接平面はユークリッド的であるから $x_\mu = x^\mu$（$\mu = 1,2,3,4$）とする．定義 2.1 により写像先の5次元座標は

$$(Y^1, Y^2, Y^3, Y^4, Y^5) = \left(\frac{Rx^1}{l},\ \frac{Rx^2}{l},\ \frac{Rx^3}{l},\ \frac{Rx^4}{l},\ \frac{R^2}{l}\right) \tag{3.1}$$

$$l = \sqrt{R^2 + \sum_{\mu=1}^4 (x^\mu)^2},\quad r^2 = \sum_{\mu=1}^4(x^\mu)^2,\quad l^2 = R^2 + r^2 \tag{3.2}$$

と書く．$(Y^1, \ldots, Y^4)$ を $(X^1, \ldots, X^4)$ と略記する．

### 3.2 引き戻し計量の導出

$S^4(R) \subset \mathbb{R}^5$ の誘導計量は，$\mathbb{R}^5$ の標準内積から

$$dS^2 = \sum_{A=1}^{5} (dY^A)^2 \tag{3.3}$$

によって定まる [6]．**第5成分を含む全5成分の寄与を計算する．**

**第1〜4成分** $(Y^\mu = R x^\mu / l)$：

**補題 3.1**  
$$\frac{\partial Y^\mu}{\partial x^\nu} = \frac{R}{l}\left(\delta^\mu_\nu - \frac{x^\mu x^\nu}{l^2}\right) \tag{3.4}$$

**証明**  
$$\frac{\partial}{\partial x^\nu}\!\left(\frac{Rx^\mu}{l}\right) = \frac{R\delta^\mu_\nu}{l} + Rx^\mu\!\left(-\frac{x^\nu}{l^3}\right) = \frac{R}{l}\!\left(\delta^\mu_\nu - \frac{x^\mu x^\nu}{l^2}\right) \quad \square$$

よって，

$$\sum_{\mu=1}^4(dY^\mu)^2 = \frac{R^2}{l^2}\sum_{\nu,\rho=1}^{4}\left(\delta_{\nu\rho} - \frac{2x_\nu x_\rho}{l^2} + \frac{r^2 x_\nu x_\rho}{l^4}\right)dx^\nu dx^\rho \tag{3.5}$$

**第5成分** $(Y^5 = R^2/l)$：

$$\frac{\partial Y^5}{\partial x^\nu} = -\frac{R^2 x^\nu}{l^3} \tag{3.6}$$

$$\left(dY^5\right)^2 = \frac{R^4}{l^6}\left(\sum_\nu x^\nu dx^\nu\right)^2 = \frac{R^4}{l^6}\,x_\nu x_\rho\,dx^\nu dx^\rho \tag{3.7}$$

**全成分の和：**

$x_\nu x_\rho$ の係数を整理する（(3.5) と (3.7) を足す）：

$$-\frac{2R^2}{l^4} + \frac{R^2 r^2}{l^6} + \frac{R^4}{l^6} = -\frac{2R^2}{l^4} + \frac{R^2(r^2+R^2)}{l^6} = -\frac{2R^2}{l^4} + \frac{R^2 l^2}{l^6} = -\frac{2R^2}{l^4} + \frac{R^2}{l^4} = -\frac{R^2}{l^4}$$

したがって，

$$dS^2 = \frac{R^2}{l^2}\delta_{\nu\rho}\,dx^\nu dx^\rho - \frac{R^2}{l^4}\,x_\nu x_\rho\,dx^\nu dx^\rho$$

$$\boxed{g_{\mu\nu} = \frac{R^2}{l^2}\left(\delta_{\mu\nu} - \frac{x_\mu x_\nu}{l^2}\right)} \tag{3.8}$$

**注意 3.1**（第5成分の不可欠性）  
第5成分 $(dY^5)^2$（式 (3.7)）を省略すると，$x_\nu x_\rho$ の係数は $-2R^2/l^4 + R^2 r^2/l^6$ のまま残り，式 (3.8) は得られない．第5成分の寄与 $+R^4/l^6 = +R^2/l^4 - R^2 r^2/l^6$ がこの差を埋める役割を果たす．

---

## 4. クリストッフェル記号とリーマン曲率

### 4.1 リーマン曲率テンソル

**定理 4.1**（定曲率空間の曲率テンソル）[3, 4, 5]  
$n$ 次元球面 $S^n(R)$ のリーマン曲率テンソルは

$$R_{\mu\nu\rho\sigma} = \frac{1}{R^2}(g_{\mu\rho}g_{\nu\sigma} - g_{\mu\sigma}g_{\nu\rho}) \tag{4.1}$$

で与えられる．

**注意 4.1**  
(4.1) は定曲率空間の標準的な結果である [4, 5]．グノモン座標での直接的な導出（クリストッフェル記号 $\Gamma^\lambda_{\mu\nu}$ からの計算）は付録Aに記す予定であるが，本報告では文献 [3, 4, 5] を引用する．

---

## 5. リッチテンソルとアインシュタインテンソル

### 5.1 リッチテンソル

**命題 5.1**  
$n$ 次元球面 $S^n(R)$ のリッチテンソルは

$$R_{\mu\nu} = \frac{n-1}{R^2}g_{\mu\nu} \tag{5.1}$$

**証明**  
リッチテンソルの定義 $R_{\mu\nu} = g^{\rho\sigma}R_{\mu\rho\nu\sigma}$ に (4.1) を代入する：

$$g^{\rho\sigma}R_{\mu\rho\nu\sigma} = \frac{1}{R^2}g^{\rho\sigma}(g_{\mu\nu}g_{\rho\sigma} - g_{\mu\sigma}g_{\nu\rho})$$

各項を計算する：

- **第1項：** $g^{\rho\sigma}g_{\mu\nu}g_{\rho\sigma} = g_{\mu\nu}\cdot\underbrace{g^{\rho\sigma}g_{\rho\sigma}}_{=\,\delta^\rho{}_\rho\,=\,n} = n\,g_{\mu\nu}$
- **第2項：** $g^{\rho\sigma}g_{\mu\sigma}g_{\nu\rho} = \underbrace{g^{\rho\sigma}g_{\mu\sigma}}_{=\,\delta^\rho_\mu}g_{\nu\rho} = g_{\nu\mu} = g_{\mu\nu}$

したがって，

$$R_{\mu\nu} = \frac{1}{R^2}(n\,g_{\mu\nu} - g_{\mu\nu}) = \frac{n-1}{R^2}g_{\mu\nu} \qquad \square$$

### 5.2 リッチスカラー

$$R_{\text{scalar}} = g^{\mu\nu}R_{\mu\nu} = \frac{n-1}{R^2}\underbrace{g^{\mu\nu}g_{\mu\nu}}_{=\,n} = \frac{n(n-1)}{R^2} \tag{5.2}$$

$n=4$ のとき $R_{\text{scalar}} = 12/R^2$．

### 5.3 アインシュタインテンソル

アインシュタインテンソルの定義 [4] より，

$$G_{\mu\nu} = R_{\mu\nu} - \frac{1}{2}g_{\mu\nu}R_{\text{scalar}} \tag{5.3}$$

$n=4$ では，

$$G_{\mu\nu} = \frac{3}{R^2}g_{\mu\nu} - \frac{1}{2}\cdot\frac{12}{R^2}\,g_{\mu\nu} = \left(\frac{3}{R^2} - \frac{6}{R^2}\right)g_{\mu\nu} = -\frac{3}{R^2}g_{\mu\nu} \tag{5.4}$$

---

## 6. 主定理

**定理 6.1**（主結果）  
4次元ユークリッド空間の接超平面 $\Pi_R$ を，曲率半径 $R$ の超球面 $S^4(R)$ へグノモン正写像したとき，写像後の誘導計量テンソル $g_{\mu\nu}$（式 (3.8)）に対して，

$$G_{\mu\nu} + \Lambda\, g_{\mu\nu} = 0 \tag{6.1}$$

が成立する．ただし，

$$\boxed{\Lambda = \frac{3}{R^2}} \tag{6.2}$$

**証明**  
(5.4) より $G_{\mu\nu} = -\dfrac{3}{R^2}g_{\mu\nu}$ であるから，$\Lambda = \dfrac{3}{R^2}$ とおけば (6.1) が直ちに成立する．$\square$

**注意 6.0**（記法について）  
本定理では $n=4$ の場合を主結果として $\Lambda$（添字なし）で表す．一般次元への拡張は注意 6.2 に $\Lambda_n$ として記した．

**注意 6.1**  
(6.1) は **Euclidean signature の $S^4(R)$ に対する幾何学的恒等式** であり，物理法則として解釈するものではない．Lorentzian signature を持つ de Sitter 時空とは本質的に別物であり，両者は解析接続（Wick 回転 $t \to it$）の関係にある（[5] Appendix D 参照）．

なお，Euclidean 版と Lorentzian 版はともに同じ $\Lambda = 3/R^2$ を持ち，Gibbons–Hawking の Euclidean 量子重力論 [8] における枠組みと一致する．この意味で，Euclidean/Lorentzian の符号の違いはWick 回転という既知の数学的手続きで橋渡しされる**テクニカルな問題**であり，本定式化が一般相対性理論の真空解と矛盾しないことを意味する．物理的含意の詳細は次年度以降に論じる．

**注意 6.2**（一般次元への拡張）  
本定理は $n=4$ に限ったものではない．一般の $n$ 次元球面 $S^n(R)$ に対しては，(5.2)(5.3) より

$$G_{\mu\nu} = \frac{n-1}{R^2}g_{\mu\nu} - \frac{1}{2}\cdot\frac{n(n-1)}{R^2}g_{\mu\nu} = -\frac{(n-1)(n-2)}{2R^2}g_{\mu\nu}$$

となり，$G_{\mu\nu} + \Lambda g_{\mu\nu} = 0$ を満たす $\Lambda$ は

$$\Lambda_n = \frac{(n-1)(n-2)}{2R^2}$$

で与えられる（$n=4$ では $\Lambda_4 = 3/R^2$）．$n=4$ すなわち $3+1$ 次元の場合が de Sitter 真空解と対応する点で物理的に際立つが，その解釈は次年度の課題とする．

---

## 7. 今後の考察予定

本報告の幾何学的結果を出発点として，次年度以降に検討する課題を列挙する（いずれも本報告の範囲外）．

- **Wick 回転の幾何学的解釈：** $t \to it$ という操作をグノモン写像の枠内でどう記述するか．本報告の Euclidean 定式化から Lorentzian 版への橋渡しとなる核心的課題
- **Lorentzian 版グノモン写像：** $l^2 = R^2 + x^2 + y^2 + z^2 - t^2$ とした場合の誘導計量の導出と符号の変化の確認
- **対称条件下の測地線：** $X^2 = X^3 = X^4 = b$，$X^1 = t$ のもとでの測地線方程式
- **対称条件 $b=t$ の幾何学的意味：** $r^2 = 3t^2$ が成立する等距離曲線の解釈
- **物質項の導入：** $T_{\mu\nu} \neq 0$ の場合への拡張
- **Euclidean/Lorentzian の対応：** 定理 6.1 と de Sitter 真空解との関係の精査（[8] 参照）

---

## 8. 結果の整理

| 項目 | 結果 | 式番号 | 出典 |
|------|------|--------|------|
| 写像の定義 | $Y^\mu = Rx^\mu/l$，$Y^5 = R^2/l$ | (3.1) | [1, 2] |
| 誘導計量 | $g_{\mu\nu} = \dfrac{R^2}{l^2}\!\left(\delta_{\mu\nu} - \dfrac{x_\mu x_\nu}{l^2}\right)$ | (3.8) | 本報告 |
| リーマン曲率 | $R_{\mu\nu\rho\sigma} = \dfrac{1}{R^2}(g_{\mu\rho}g_{\nu\sigma}-g_{\mu\sigma}g_{\nu\rho})$ | (4.1) | [3,4,5] |
| リッチテンソル | $R_{\mu\nu} = \dfrac{3}{R^2}g_{\mu\nu}$ | (5.1) | 本報告 |
| リッチスカラー | $R_{\text{scalar}} = \dfrac{12}{R^2}$ | (5.2) | 本報告 |
| アインシュタインテンソル | $G_{\mu\nu} = -\dfrac{3}{R^2}g_{\mu\nu}$ | (5.4) | 本報告 |
| **主定理** | $G_{\mu\nu} + \Lambda g_{\mu\nu} = 0$，$\Lambda = \dfrac{3}{R^2}$ | (6.2) | 本報告 |

---

## 参考文献

[1] Snyder, J. P., *Map Projections: A Working Manual*, U.S. Geological Survey Professional Paper 1395, U.S. Government Printing Office, Washington D.C., 1987, pp. 164–168.

[2] Coxeter, H. S. M., *Introduction to Geometry*, 2nd ed., Wiley, New York, 1969, Chapter 15.

[3] do Carmo, M. P., *Riemannian Geometry*, Birkhäuser, Boston, 1992, Chapter 4.

[4] Misner, C. W., Thorne, K. S., Wheeler, J. A., *Gravitation*, W. H. Freeman, San Francisco, 1973, §13.5.

[5] Wald, R. M., *General Relativity*, University of Chicago Press, Chicago, 1984, Appendix D.

[6] Spivak, M., *A Comprehensive Introduction to Differential Geometry*, Vol. 2, 3rd ed., Publish or Perish, Houston, 1999, Chapter 4.

[7] Lee, J. M., *Introduction to Riemannian Manifolds*, 2nd ed., Springer, New York, 2018, Example 3.14.

[8] Gibbons, G. W. and Hawking, S. W., "Action integrals and partition functions in quantum gravity," *Physical Review D*, **15**, 2752–2756, 1977.  
　　（Euclidean 量子重力論の基礎；Euclidean de Sitter 時空と $\Lambda = 3/R^2$ の関係）

---

## 付録A：クリストッフェル記号の直接計算（次年度予定）

計量 (3.8) からクリストッフェルの公式によって $\Gamma^\lambda_{\mu\nu}$ を陽に計算し，定理 4.1 のリーマン曲率テンソルを直接導出する計算を，次年度の報告において完成させる予定である．現時点では [3, 4] の結果を引用する．

---

## 付録B：$R \to \infty$ の平坦極限

$R \to \infty$ のとき $l \to R$ であるから，

$$g_{\mu\nu} \to \delta_{\mu\nu}, \quad \Lambda = \frac{3}{R^2} \to 0$$

計量は平坦なユークリッド計量に収束し，曲率はゼロに収束する．これはグノモン写像が恒等写像に退化することに対応する（極限操作の詳細については [3] Chapter 4 を参照）．

---

## 付録C：写像の正則性

**命題 C.1**（グノモン正写像の正則性）  
$\Phi: \mathbb{R}^4 \to S^4(R)$ を定義 2.1 のグノモン正写像とする．$r = \sqrt{\sum_\mu (x^\mu)^2}$ とおくとき，

$$Y^5 = \frac{R^2}{l} = \frac{R^2}{\sqrt{R^2 + r^2}}, \quad Y^\mu = \frac{Rx^\mu}{\sqrt{R^2+r^2}}$$

に対して以下が成立する：

1. **$R \to 0$（$r > 0$ 固定）：** $Y^5 = R^2/\sqrt{R^2+r^2} \to 0$，$Y^\mu \to 0$．写像像は原点に縮退するが発散しない．

2. **$R \to \infty$（$r$ 固定）：** $Y^5 = R/\sqrt{1+r^2/R^2} \to \infty$．第5成分が発散する．一方 $Y^\mu \to x^\mu$（有限）．

3. **固定 $R > 0$：** $l = \sqrt{R^2+r^2} \geq R > 0$ であるから分母は常に正であり，$\Phi$ は $\mathbb{R}^4$ 全域で $C^\infty$ 級の正則写像である．

**証明**  
いずれも $l = \sqrt{R^2+r^2}$ の直接計算による．(1) $R\to 0$：$Y^5 = R^2/l \leq R^2/r \to 0$．(2) $R\to\infty$：$Y^5 = R\cdot(1+r^2/R^2)^{-1/2} \to \infty$．(3) $R>0$ のとき $l \geq R > 0$ であるから各成分は $r$ に関して $C^\infty$．$\square$

**注意 C.1**  
命題 C.1-(2) の $R \to \infty$ における $Y^5 \to \infty$ は，第5成分（埋め込み空間の余剰次元方向）が発散することを示す．これは付録Bで示した「$R\to\infty$ で計量 $g_{\mu\nu} \to \delta_{\mu\nu}$（平坦極限）」と整合的である：球面が「無限に大きな球」になる極限で，接平面との差異が消えると同時に，埋め込み次元方向への射影が無限遠に退くことに対応する．
---

## 付録D：Wick 回転不要性——2乗計量の本質的観察

本報告の定式化を精査すると，Lorentzian 計量への拡張において Wick 回転（$t \to it$）が**そもそも不要**であることが見えてくる．以下にその論拠を純粋に数学的事実として記す．

**観察 D.1**（物理計量の2乗性）  
グノモン写像に現れる全ての計量は，$t$ 自身ではなく $t^2$ の形でのみ登場する：

$$l^2 = R^2 + x^2 + y^2 + z^2 \pm t^2, \quad ds^2 = \pm dt^2 + dx^2 + dy^2 + dz^2$$

$t^2,\ x^2,\ y^2,\ z^2$ はいずれも**常に非負の実数**である．符号の違いはこれらの2乗量に付与される**係数 $\epsilon_\mu \in \{+1,-1\}$ の選択**に過ぎない．$t$ 自身は補助的な座標量であり，物理計量として本質的なのは常に $t^2$ である．

**命題 D.1**（Euclidean/Lorentzian の統一的記述）  
符号系 $(\epsilon_1, \epsilon_2, \epsilon_3, \epsilon_4) \in \{+1,-1\}^4$ を選ぶ．

**埋め込み空間の計量**を $\mathbb{R}^5$ の標準計量ではなく，

$$dS^2_{\rm emb} = \sum_{\mu=1}^{4} \epsilon_\mu (dY^\mu)^2 + (dY^5)^2 \tag{D.1}$$

と定義する（第5成分 $Y^5$ の係数は $+1$ に固定）．

接超平面 $\Pi_R$ 上の計量を

$$ds^2_{\rm tan} = \sum_{\mu=1}^{4} \epsilon_\mu (dx^\mu)^2 \tag{D.2}$$

とし，$l^2 = R^2 + \sum_{\mu=1}^4 \epsilon_\mu (x^\mu)^2$ とおく．グノモン写像 $\Phi$ の定義は本文 (2.1) と同じく

$$Y^\mu = \frac{R x^\mu}{l} \quad (\mu=1,2,3,4), \quad Y^5 = \frac{R^2}{l}$$

とする．このとき $\Phi$ によって引き戻された誘導計量は

$$g_{\mu\nu} = \frac{R^2}{l^2}\left(\epsilon_\mu \delta_{\mu\nu} - \frac{\epsilon_\mu \epsilon_\nu x_\mu x_\nu}{l^2}\right) \tag{D.3}$$

で与えられる（添え字の繰り返しに Einstein 規約は用いない）．

**証明**  
(D.1) のもとで $dS^2_{\rm emb}$ を接平面の座標で表す．

$dY^\mu$ は本文補題 3.1 と同様に計算され，

$$\sum_{\mu=1}^4 \epsilon_\mu (dY^\mu)^2 = \frac{R^2}{l^2} \sum_{\nu,\rho=1}^4 \left(\epsilon_\nu \delta_{\nu\rho} - \frac{2\epsilon_\nu \epsilon_\rho x_\nu x_\rho}{l^2} + \frac{\epsilon_\nu \epsilon_\rho r_\epsilon^2 x_\nu x_\rho}{l^4}\right) dx^\nu dx^\rho \tag{D.4}$$

ただし $r_\epsilon^2 = \sum_\mu \epsilon_\mu (x^\mu)^2 = l^2 - R^2$．

第5成分は $Y^5 = R^2/l$ より，

$$(dY^5)^2 = \frac{R^4}{l^6}\left(\sum_\nu \epsilon_\nu x_\nu dx^\nu\right)^2 = \frac{R^4}{l^6}\,\epsilon_\nu \epsilon_\rho\, x_\nu x_\rho\, dx^\nu dx^\rho \tag{D.5}$$

$\epsilon_\nu \epsilon_\rho\, x_\nu x_\rho$ の係数を (D.4) と (D.5) で合算する：

$$-\frac{2R^2}{l^4} + \frac{R^2 r_\epsilon^2}{l^6} + \frac{R^4}{l^6} = -\frac{2R^2}{l^4} + \frac{R^2(r_\epsilon^2 + R^2)}{l^6} = -\frac{2R^2}{l^4} + \frac{R^2 l^2}{l^6} = -\frac{R^2}{l^4}$$

これは $r_\epsilon^2 + R^2 = l^2$ を用いた（本文 3.2 節と同じ相殺機構）．したがって，

$$dS^2_{\rm emb} = \frac{R^2}{l^2}\epsilon_\nu \delta_{\nu\rho}\, dx^\nu dx^\rho - \frac{R^2}{l^4}\epsilon_\nu \epsilon_\rho\, x_\nu x_\rho\, dx^\nu dx^\rho$$

よって (D.3) が成立する．$\square$

**注意 D.2**（相殺機構が機能する条件）  
証明の鍵は $r_\epsilon^2 + R^2 = l^2$，すなわち $l^2$ の定義と $r_\epsilon^2$ の定義が**同じ $\epsilon_\mu$** を用いることである．埋め込み空間の計量 (D.1) を $\epsilon_\mu$ で変えずに接平面の計量のみ変えた場合，この等式が崩れ相殺が機能しない（これが「単純な符号の書き換えでは不十分」な理由であり，埋め込み空間の計量の同時変更が不可欠である）．

**注意 D.3**（Wick 回転との対応）  
Wick 回転 $t \to it$ は，この枠組みでは

$$\epsilon_4 : +1 \to -1$$

という係数符号の変更として記述される（虚数 $i$ は一切不要）：

$$\text{Euclidean}: \quad \epsilon_\mu = (+1,+1,+1,+1) \quad \Rightarrow \quad l^2 = R^2 + x^2 + y^2 + z^2 + t^2$$

$$\text{Lorentzian}: \quad \epsilon_\mu = (+1,+1,+1,-1) \quad \Rightarrow \quad l^2 = R^2 + x^2 + y^2 + z^2 - t^2$$

両者は埋め込み空間の計量 (D.1) と写像の定義を共有する**単一の幾何学的枠組み**の特殊ケースであり，「2つの別理論」ではない．

**注意 D.4**（本付録の位置づけと超曲面の確認）  
命題 D.1 は埋め込み空間の計量を (D.1) と定義した場合の数学的帰結である．

写像の像が乗る超曲面を確認する：

$$\sum_{\mu=1}^4 \epsilon_\mu (Y^\mu)^2 + (Y^5)^2 = \frac{R^2}{l^2}\,r_\epsilon^2 + \frac{R^4}{l^2} = \frac{R^2(r_\epsilon^2 + R^2)}{l^2} = \frac{R^2 \cdot l^2}{l^2} = R^2 \tag{D.6}$$

すなわち写像の像は常に半径 $R$ の超曲面上に乗る．符号選択による超曲面の変化：

$$\epsilon_\mu = (+1,+1,+1,+1): \quad (Y^1)^2+(Y^2)^2+(Y^3)^2+(Y^4)^2+(Y^5)^2 = R^2 \quad \text{（$S^4(R)$ 球面，Euclidean）}$$

$$\epsilon_\mu = (+1,+1,+1,-1): \quad (Y^1)^2+(Y^2)^2+(Y^3)^2-(Y^4)^2+(Y^5)^2 = R^2 \quad \text{（擬球面，Lorentzian）}$$

後者は $\mathbb{R}^{4,1}$ における擬球面であり，de Sitter 超曲面と構造的に対応する．この対応の厳密な検証は次年度以降の課題とする．

---

*本報告は幾何学的事実の整理を目的とし，物理的解釈は行わない．全数式の出典を参考文献に明示した．宇宙物理学講座教授および数学講座准教授による3回の査読指摘を経て改訂した最終版である．*
