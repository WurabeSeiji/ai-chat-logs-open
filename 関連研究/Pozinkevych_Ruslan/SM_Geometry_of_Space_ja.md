> 原文: 英語 | 翻訳: 日本語（機械翻訳）

---

# 素粒子の標準模型と空間の幾何学 ── 補遺および訂正付き

**Standard Model of Elementary Particles and the Geometry of Space with Additions and Corrections**

**著者:** Ruslan Pozinkevych
Faculty of Informations Technologies and Mathematics,
The Eastern European National University, Ukraine,
43021, Lutsk, Potapov str.9
galagut@protonmail.com

**日付:** 2025年6月30日
**HAL Id:** hal-05136544

---

## 概要（Abstract）

**目的:** 本研究の目的は、素粒子物理学の標準模型（SM）の妥当性に対する数学的根拠を提示することである。理論の分析はSMのほぼすべての側面を網羅しているが、ヒッグス粒子については例外となる可能性がある。ヒッグス粒子の追加は、実験データによって後に調整されるか、本稿で示される方程式の因子として組み込まれる可能性がある。SMが長年にわたり証明しようとしてきた空間の対称性は、その不均一性ゆえに数学的に反証される可能性がある。一方で、本稿が提示するアプローチにより、研究者は宇宙規模での質量密度とその分布を計算することが可能となる。著者の見解では、これが宇宙の物質分布が不均一である理由である。四つの既知の基本相互作用の間の単純さと関連性が、我々の知る空間の幾何学に影響を与える要因として考慮される。いわゆる宇宙における空間の「織物（ファブリック）」を可視化することさえ可能であり、すべての既知の素粒子をそのエネルギー準位に基づいて分類することもできる。

**キーワード:** 正方行列、素粒子物理学、標準模型、空間、重力理論

---

## 1. 一般情報

行列の性質に関するよく知られた事実から始める：

$$Q^t Q = QQ^t = I \quad \cdots (1)$$

この性質は正方行列にのみ適用される。すなわち、2×2, 3×3, ...n×n の型の行列である。これにより、二つの行列の外積（クロス積）が単位行列を与えるようにデータをグループ化することが可能となる。以下を確認されたい：

与えられた行列の転置はその逆行列と明らかに等しくない：

$$Q^t \neq Q^{-1} \rightarrow Q^t Q \neq QQ^t \neq I$$

したがって、$QQ^t = I$ を満たすような行列/ベクトルを探す必要がある。

$$\text{diag}(e, g, c, h) \times \text{diag}(1/e, 1/g, 1/c, 1/h) = I \quad \cdots (2)$$

この規則は、定義により $|A|^{n-2}A$ として、すべての正方対角行列に適用される。しかし、本研究の目的は、数学的規則が物理学の世界にどのように適用されるかを示すことである。

---

## 2. 素材と方法

しばらく正方行列の性質に注目しよう。4×4行列の行列式は、余因子とそれぞれの三角行列の積の和として定義される。以下の例を見てみよう：

$$\begin{vmatrix} x-1 & y-2 & z-3 & d-4 \\ 5 & 6 & 7 & 8 \\ 7 & 10 & 11 & 12 \\ 13 & 14 & 15 & 16 \end{vmatrix} \quad \cdots (3)$$

行列(3)の行列式は $16d + 16y - 32z$ に等しい。言い換えれば、方程式(3)の結果は R³空間を定義する3成分ベクトルである。同じ行列が異なる空間を定義することは可能だろうか？ 可能である。以下を見てみよう：

$$\begin{vmatrix} x-1 & y-2 & z-3 & d-4 \\ 5 & 6 & 7 & 8 \\ 7 & 10 & 11 & 2 \\ 13 & 14 & 15 & 16 \end{vmatrix} \quad \cdots (4)$$

同じサイズの行列の行列式は、今度は R⁴空間を定義する4成分ベクトル $16d - 80x + 176y - 112z$ となる。

---

## 3. 方法論

ここで、行列を形成する4つのベクトルの和を相関させ、式(3)と(4)で使用したのと同じ数値を代入して行列式を計算すると：

$$\begin{pmatrix} x-e & y-g & z-h & d-c \\ E_1 = k_e|q|/r^2 & E_1 = mgh & E_1 = h\nu & E_1 = mc^2 \\ E_2 = k_e|q|/r^2 & E_2 = mgh & E_2 = h\nu & E_2 = mc^2 \\ E_3 = k_e|q|/r^2 & E_3 = mgh & E_3 = h\nu & E_3 = mc^2 \end{pmatrix} \quad \cdots (5)$$

最初の場合、結果は $16d + 16y - 32z$ であるが、第二の場合はこの結果が異なり $16d - 80x + 176y - 112z$ となる。エネルギー保存則および一般的な保存則に関して、入力エネルギー量は出力と等しくなければならないが、同じ入力方程式 $(x-e)(y-g)(z-h)(d-c)$ は、使用する次元に応じて数値的に異なる出力を与える。入力として異なるが、エネルギー的に同じ出力値を生み出す一連のエントリーを見つけることは常に可能である。

---

## 4. 証明

様々な次元で測定されるエネルギーの量は異なるが、様々な測定手法を用いれば同じ値となる可能性がある。エネルギー分布は変動しうるため、その密度について議論すべきである。式(3)と(4)を使用しても、余因子を用いて等式化することが可能である：

$$((x-1), -(y-2), +(z-3), -(d-4))$$

より一般的には：

$$((x-e), -(y-g), +(z-h), -(d-c))$$

---

## 5. 結論

標準模型は理論的に自己無撞着であり、いくつかの実験的予測を行うことができると考えられている。ベクトル解析と組み合わせることで、物質/反物質パラドックスのように説明できない物理現象を説明し、一般相対性理論が記述する重力理論を完全に組み込み、暗黒エネルギーによって説明されうる宇宙の加速膨張を説明することが可能となる。本モデルは、観測宇宙論から推論されるすべての必要な特性を持つデータを含む。さらに、波動理論やミクロおよびマクロ物理現象を特徴づけるその他の方法を組み込むことが可能である。一方で、$a_{5,5}$ 要素がヒッグス粒子スカラーであるという追加次元を与えたパターンに対して、余分な次元を加えた式を提示することで、このモデルを拡張することができる（以下の式を参照）：

$$\Delta \begin{vmatrix} h & e & c & g & \cdots \\ & & \cdots & & \\ & \vdots & \ddots & \vdots & \end{vmatrix} = |\vec{E}| = \vec{E} \quad : \quad a_{5,5} = H \text{ (ヒッグス粒子)} \quad \cdots (6)$$

逆に、エネルギー保存則を別の観点から見ることもできる。すなわち、「エネルギー」の和の代わりに、式(2)と同様に、それらのスカラー表現または行列積を使用する：

$$\text{diag}(e, g, c, h) \times \text{diag}(1/e, 1/g, 1/c, 1/h) = I$$

「...ではヒッグス粒子はどこにあるのか？」という疑問が生じるかもしれない。まさにここにある。行列 $A$ の逆行列がヒッグス粒子である。これを踏まえて、「エネルギー保存則」を次のように定式化しよう：「エネルギーの量は一定であり、以下に等しい：

$$A \cdot U = I \quad : \quad U \text{ は行列 } A \text{ の逆行列} \quad \cdots (7)$$

さらに分析すると、$A$ は変数 $a_i$ の行列、$U$ は係数 $b_i$ の行列となりうることがわかる。以下を確認されたい：

$$\sum_{i=0}^{n} a_i b_i - 1 = 0 \quad \cdots (8)$$

---

## 参考文献

[1] Emana, T. K., & Haramaya, E. (2017). FEYNMAN DIAGRAMS OF THE STANDARD MODEL.

[2] Feynman, R. (2018). *Feynman lectures on gravitation*. CRC press.

[3] Aitken, A. C. (2017). *Determinants and matrices*. Read Books Ltd.

[4] Day, J., & Peterson, B. (1988). Growth in Gaussian elimination. *The American mathematical monthly*, 95(6), 489-513.

[5] Huss-Lederman, S., Jacobson, E. M., Tsao, A., Turnbull, T., & Johnson, J. R. (1996). Implementation of Strassen's algorithm for matrix multiplication. In *Proceedings of the 1996 ACM/IEEE Conference on Supercomputing* (pp. 32-es).

[6] Boas, M. L., & Spector, D. (1999). *Mathematical methods in the physical sciences*.

[7] Bapat, R. B. (2010). *Graphs and matrices* (Vol. 27, p. 23). New York: Springer.

[8] Hayter, S. W. (1976). Almost All about Waves by John R. Pierce. *Leonardo*, 9(3), 246-246.

[9] Tobisch, E. (Ed.). (2016). *New Approaches to Nonlinear Waves*. Cham: Springer International Publishing.

[10] Pozinkevych, R. (2020). Ternary Mathematics Principles Truth Tables and Logical Operators 3 D Placement of Logical Elements Extensions of Boolean Algebra. *In other words*, 1(1), 0.

[11] Yuan, T. (2024). Physics Mechanisms Behind Light and Gravity in Universe.

[12] Dolciani, M. P. (1967). Modern introductory analysis. In *Modern introductory analysis* (pp. 660-660).

[13] Hošková-Mayerová, Š., Flaut, E. C., & Maturo, F. (Eds.). (2021). *Algorithms as a Basis of Modern Applied Mathematics*. Springer.

[14] Dowdy, S., Wearden, S., & Chilko, D. (2011). *Statistics for research*. John Wiley & Sons.

[15] Dubovyk, V. P., & Yuryk, I. I. (2011). *Vyshcha matematyka. Zbirnyk zadach: navch. posib.*

[16] Mosquera Rodas, Jhon Jairo. "Systemic Relational Lineality Theory and the Relational Simplification Method for Heterogeneous Equations in Physics and Mathematics." *Available at SSRN 4793404* (2024).
