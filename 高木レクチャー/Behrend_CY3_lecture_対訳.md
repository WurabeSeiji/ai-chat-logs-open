# Kai Behrend「Enumerative Geometry of Calabi–Yau threefolds」対訳

> 第26回高木レクチャー（2026-06-06/07, RIMS）配布講義ノート（Kai Behrend, May 6, 2026）の英日対訳。配布冊子（全7節）からの転記。**§1〜§7 全節 転記完了**。
> 数式は LaTeX で原文を保存。和訳は学習用。出典写真：`~/Downloads/IMG_7067〜7101`。

---

## Contents（目次）

1. Calabi–Yau threefolds … 2
2. Enumerative Geometry … 7
3. Symplectic geometry … 11
4. Critical loci … 13
5. Motivic counts … 18
6. Pardon's classification of curve counting theories … 20
7. A non-commutative Calabi–Yau threefold … 24

---

## Introduction（序）

**EN.** This is a brief tour of some aspects of the enumerative geometry of Calabi–Yau threefolds which I personally find interesting. After briefly recalling what CY3s are, and giving a few examples, we come to the main question, namely why CY3s are uniquely interesting from the point of view of enumerative geometry. We review some of the early highlights of the story, and give a few remarks on the definition of enumerative invariants of various flavours. We briefly get into derived symplectic geometry, as that is the real reason for our interest in CY3s. We explain the fact that moduli spaces associated to CY3s locally look like critical loci, and what consequences this has for the enumerative geometry. We mention the recent result of Pardon, classifying all enumerative theories for CY3s and proving the famous Gromov–Witten Donaldson–Thomas correspondence. We close with a discussion of an example I am particularly fond of, namely a non-commutative incarnation of the mirror quintic.

**JA.** 本稿は、私が個人的に面白いと思うカラビ・ヤウ3次元体（CY3）の数え上げ幾何のいくつかの側面を、手短に巡るものである。まず CY3 とは何かを簡単に復習し、いくつか例を挙げたうえで、主問題——すなわち、数え上げ幾何の観点から見て CY3 がなぜ唯一無二に興味深いのか——に至る。物語の初期のハイライトをいくつか振り返り、種々の流儀の数え上げ不変量の定義についていくつか注意を述べる。続いて**導来シンプレクティック幾何**に簡単に立ち入る——これこそが、我々が CY3 に関心を持つ真の理由だからである。CY3 に付随するモジュライ空間が局所的に**臨界軌跡（critical loci）**のように見えること、そしてそれが数え上げ幾何にどんな帰結をもたらすかを説明する。最後に、すべての CY3 の数え上げ理論を分類し、有名な**グロモフ・ウィッテン＝ドナルドソン・トーマス対応**を証明した Pardon の最近の結果に触れる。締めくくりに、私が特に気に入っている例——**ミラー・クインティックの非可換版**——を論じる。

---

## 1. Calabi–Yau threefolds（カラビ・ヤウ3次元体）

### Definition（定義）

**EN.** A *Calabi–Yau threefold* is a complex algebraic variety $X$, satisfying the following conditions. This means it is the common zero locus of a finite number of homogeneous equations in $n+1$ variables
$$F_1(x_0,\dots,x_n),\ \dots,\ F_k(x_0,\dots,x_n),$$
in $\mathbb{P}^n(\mathbb{C})$.
- (i) $X$ is projective.
- (ii) $X$ is smooth, or non-singular, i.e., a compact complex manifold,
- (iii) $\dim_{\mathbb{C}} X = 3$,
- (iv) there exists a global, nowhere vanishing holomorphic differential 3-form, in other words, a global section of $\omega_X = \Omega^3_X$ (a holomorphic orientation); this is considered part of the structure,
- (v) (strict Calabi–Yau) $H^1(X,\mathcal{O}_X) = H^2(X,\mathcal{O}_X) = 0$, or equivalently, $H^0(\Omega) = H^0(\Omega^2) = 0$. This condition says the Calabi–Yau threefold is not special in the sense that it could decompose, for example, into a product of an elliptic curve (a Calabi–Yau 1-fold) and a K3-surface (a Calabi–Yau 2-fold).

**JA.** **カラビ・ヤウ3次元体**とは、以下の条件を満たす複素代数多様体 $X$ のことである。すなわち $X$ は、$n+1$ 変数の有限個の斉次方程式
$$F_1(x_0,\dots,x_n),\ \dots,\ F_k(x_0,\dots,x_n)$$
の $\mathbb{P}^n(\mathbb{C})$ における共通零点集合である。
- (i) $X$ は射影的である。
- (ii) $X$ は滑らか（非特異）、すなわちコンパクト複素多様体である。
- (iii) $\dim_{\mathbb{C}} X = 3$。
- (iv) 大域的でどこでも消えない正則微分3形式が存在する。言い換えれば $\omega_X = \Omega^3_X$ の大域切断（正則な向き付け）が存在する。これは構造の一部とみなす。
- (v) （狭義カラビ・ヤウ）$H^1(X,\mathcal{O}_X) = H^2(X,\mathcal{O}_X) = 0$、同値に $H^0(\Omega) = H^0(\Omega^2) = 0$。この条件は、CY3 が——例えば楕円曲線（カラビ・ヤウ1次元体）と K3 曲面（カラビ・ヤウ2次元体）の積に分解する、といった意味で——特殊でないことを言う。

**EN.** Because $X$ is projective we can turn $X$ into a Kähler manifold by choosing an appropriate metric. By Calabi's conjecture (Yau's theorem) $X$ admits, in fact, a Ricci-flat Kähler metric. This reduces the holonomy group from $U(3)$ to $SU(3)$.

**JA.** $X$ は射影的なので、適当な計量を選べば $X$ をケーラー多様体にできる。**カラビ予想（ヤウの定理）**により、$X$ は実は**リッチ平坦なケーラー計量**を持つ。これにより、ホロノミー群が $U(3)$ から $SU(3)$ へ縮小する。

### An example: the quintic（例：クインティック）

**EN.** The most famous example of a Calabi–Yau threefold is the quintic. It is defined by a single equation of degree 5 in 5 variables, defining a hypersurface $X$ in $\mathbb{P}^4(\mathbb{C})$. For this hypersurface to be non-singular and strict, the coefficients need to be generic.
$$\sum_{\substack{i_0,i_1,i_2,i_3,i_4\ge 0\\ i_0+i_1+i_2+i_3+i_4=5}} \alpha_{i_0,i_1,i_2,i_3,i_4}\, x_0^{i_0} x_1^{i_1} x_2^{i_2} x_3^{i_3} x_4^{i_4} = 0.$$
The quintic has $\binom{9}{5} = 126$ coefficients, up to scaling only 125. The group $PGL(5)$ of dimension 24 acts on these equations by projective linear coordinate changes, and so up to linear projective coordinate changes there is a 101-dimensional family of quintics. This is reflected in
$$h^{1,2} = \dim H^1(X,T_X) = \dim H^2(X,\Omega_X) = 101,$$
which indicates that all infinitesimal deformations of a quintic as a compact complex manifold are realized as quintics. All non-singular quintics are strict CY3s.

**JA.** CY3 の最も有名な例は**クインティック（5次超曲面）**である。これは5変数・次数5の単一方程式で定義され、$\mathbb{P}^4(\mathbb{C})$ 内の超曲面 $X$ を定める。この超曲面が非特異かつ狭義であるためには、係数が generic（一般的）である必要がある。
$$\sum_{\substack{i_0,\dots,i_4\ge 0\\ i_0+\dots+i_4=5}} \alpha_{i_0,\dots,i_4}\, x_0^{i_0} x_1^{i_1} x_2^{i_2} x_3^{i_3} x_4^{i_4} = 0.$$
クインティックは $\binom{9}{5} = 126$ 個の係数を持ち、スケール変換を除けば 125 個。次元 24 の群 $PGL(5)$ が射影線型座標変換でこれらの方程式に作用するので、線型射影座標変換を除けば、クインティックは **101 次元の族**をなす。これは
$$h^{1,2} = \dim H^1(X,T_X) = \dim H^2(X,\Omega_X) = 101$$
に反映される。これは、コンパクト複素多様体としてのクインティックのすべての無限小変形がクインティックとして実現されることを示す。すべての非特異クインティックは狭義 CY3 である。

> 〔木原メモ〕係数の数え上げ $\binom{9}{5}=126$（＝5変数・次数5の単項式数 $\binom{5+5-1}{5}$）、スケールを割って 125、$PGL(5)$（24次元）を割って 101 次元の族。**「格子点（指数ベクトル $i_0+\dots+i_4=5$）の数え上げ」**が次元を決める＝Ehrhart/格子点計数の典型。α=137 の格子点計数と同じ精神。

### Another example: double covers of $\mathbb{P}^3(\mathbb{C})$（例：$\mathbb{P}^3$ の二重被覆）

**EN.** Another example are double covers of projective space, branched along an octic surface. (These are also hypersurfaces of degree 8 in weighted projective 4-space with weights $(1,1,1,1,4)$.)
$$y^2 = F(x_0, x_1, x_2, x_3)$$
These form a family of dimension 164. The group of coordinate changes in $\mathbb{P}^3(\mathbb{C})$ has dimension 15, giving a dimension of 149 for the family of these double covers. This agrees with
$$h^{1,2} = \dim H^1(X,T_X) = \dim H^2(X,\Omega_X) = 149,$$
so again all complex deformations of such a CY3 are realized as double covers branched along an octic.

**JA.** もう一つの例は、**8次曲面（octic）に沿って分岐する射影空間の二重被覆**である。（これらは重み $(1,1,1,1,4)$ の重み付き射影4次元空間内の次数8の超曲面でもある。）
$$y^2 = F(x_0, x_1, x_2, x_3)$$
これらは次元 164 の族をなす。$\mathbb{P}^3(\mathbb{C})$ の座標変換の群は次元 15 なので、これら二重被覆の族の次元は 149 となる。これは
$$h^{1,2} = \dim H^1(X,T_X) = \dim H^2(X,\Omega_X) = 149$$
と一致し、やはりこの種の CY3 のすべての複素変形が、octic に沿って分岐する二重被覆として実現される。

### The Dwork pencil（ドワーク鉛筆束）

**EN.** This is a 1-dimensional family of quintics with parameter $\psi \in \mathbb{C}$:
$$x_0^5 + x_1^5 + x_2^5 + x_3^5 + x_4^5 - 5\psi\, x_0 x_1 x_2 x_3 x_4 = 0. \tag{1}$$
The quintic $X_\psi$ is non-singular if $\psi \neq \infty$ and $\psi^5 \neq 1$. Consider the group
$$G = \left\{ (\alpha_0,\dots,\alpha_4) \in (\mu_5)^5 \ \middle|\ \prod_{i=0}^{4} \alpha_i = 1 \right\} \Big/ \mu_5,$$
where $\mu_5 \subset \mathbb{C}$ is the group of 5-th roots of unity. The group $G$ has 125 elements, and is isomorphic to $(\mathbb{Z}/5)^3$. It acts on $\mathbb{P}^4(\mathbb{C})$ by
$$(\alpha_0,\dots,\alpha_4) \cdot \langle x_0,\dots,x_4\rangle = \langle \alpha_0 x_0,\dots,\alpha_4 x_4\rangle,$$
leaving Equation (1) unchanged. Hence $G$ acts on the solutions to (1), preserving the holomorphic volume form, and the quotient $X_\psi/G$ is a singular …〔続く〕

**JA.** これはパラメータ $\psi \in \mathbb{C}$ を持つクインティックの **1次元族**である：
$$x_0^5 + x_1^5 + x_2^5 + x_3^5 + x_4^5 - 5\psi\, x_0 x_1 x_2 x_3 x_4 = 0. \tag{1}$$
クインティック $X_\psi$ は、$\psi \neq \infty$ かつ $\psi^5 \neq 1$ のとき非特異である。次の群を考える：
$$G = \left\{ (\alpha_0,\dots,\alpha_4) \in (\mu_5)^5 \ \middle|\ \prod_{i=0}^{4} \alpha_i = 1 \right\} \Big/ \mu_5,$$
ここで $\mu_5 \subset \mathbb{C}$ は1の5乗根のなす群。群 $G$ は **125 個の元**を持ち、$(\mathbb{Z}/5)^3$ に同型である。$G$ は $\mathbb{P}^4(\mathbb{C})$ に
$$(\alpha_0,\dots,\alpha_4) \cdot \langle x_0,\dots,x_4\rangle = \langle \alpha_0 x_0,\dots,\alpha_4 x_4\rangle$$
で作用し、方程式 (1) を不変に保つ。ゆえに $G$ は (1) の解集合に、正則体積形式を保ちながら作用し、商 $X_\psi/G$ は特異な …〔続く〕

> 〔木原メモ〕Dwork pencil ＋ $G\cong(\mathbb{Z}/5)^3$（125元）による商＝**ミラー・クインティックの構成（orbifold）**。前に議論した「Q軸＝有限群（離散ラベル）＝orbifold／discrete torsion」「CY の組合せ・トーリック側」と直結。$(\mathbb{Z}/5)^3$ の有限群作用で CY を作る＝離散内部対称性の幾何化の実例。

---

### The Dwork pencil（続き）— 特異点配置とクレパント解消

**EN.** … CY3. The singularities form along 10 lines, which intersect in 10 points; at every one of these points, 3 lines intersect, each of the lines passes through 3 of the points. 〔図：10直線・10点の配置（各点で3直線が交わり、各直線が3点を通る＝$(10_3)$ 配置）〕

One can resolve the singularities in a way preserving the Calabi–Yau structure (a so-called *crepant* resolution) yielding the 1-parameter family of mirror quintics.

Along the lines the singularity is a transverse $A_4$ surface singularity, with equation $uv = w^5$, which can be resolved by 4 blow-ups in succession, creating for each of the 10 curves 4 ruled surfaces, for a total of 40 ruled surfaces. The singularities at the 10 points are worse; when resolving them, each one creates 6 exceptional divisors. In total we get 100 exceptional divisors, plus the hyperplane this gives a total of 101 independent homology classes in $H_4(X,\mathbb{Q})$, by Poincaré duality 101 classes in $H^2(X,\mathbb{Q})$, and the dimension is $\dim H^2(X,\mathbb{Q}) = h^{1,1}$, yielding $h^{1,1} = 101$.

Together with $h^{1,2} = 1$, this dimension count would justify suspecting the resolved Dwork quotient as the mirror quintic.

**JA.** … は CY3 である。特異点は **10本の直線**に沿って生じ、それらは **10個の点**で交わる。各点では **3直線が交わり**、各直線は**3個の点を通る**（＝$(10_3)$ 配置）。

カラビ・ヤウ構造を保つように特異点を解消できる（いわゆる**クレパント解消（crepant resolution）**）と、**ミラー・クインティックの1パラメータ族**が得られる。

直線に沿った特異点は横断的な **$A_4$ 曲面特異点**で、方程式 $uv = w^5$ を持ち、4回の連続ブローアップで解消でき、10本の各曲線に対し4枚の線織面を作る——計 **40枚の線織面**。10個の点での特異点はより悪く、解消するとそれぞれ **6個の例外因子**を作る。合計で **100個の例外因子**、さらに超平面を加えて $H_4(X,\mathbb{Q})$ の **101個**の独立なホモロジー類、ポアンカレ双対で $H^2(X,\mathbb{Q})$ の101類、その次元 $\dim H^2(X,\mathbb{Q}) = h^{1,1}$、すなわち **$h^{1,1} = 101$**。

$h^{1,2} = 1$ と合わせて、この次元勘定から、解消した Dwork 商を**ミラー・クインティック**と疑う根拠が得られる。

> 〔木原メモ〕**101 = 40（線織面）＋ 60（10点×6）＋ 1（超平面）**＝**例外因子の数え上げ（組合せ）が $h^{1,1}$ を決める**。$(10_3)$ 配置（Desargues 配置）＝純粋な接続組合せ論。これも「数え上げが次元を決める」典型。

### Fermat quintic の不変環 → ミラーの記述

**EN.** Let us consider the quintic mirror for $\psi = 0$. It is not hard to see that the ring of invariants of $G$ acting on the homogeneous coordinate ring of the Fermat quintic
$$\mathbb{C}[x_0,\dots,x_4]/(x_0^5 + \dots + x_4^5)$$
is generated by $x_0^5,\dots,x_4^5$ and $x_0 x_1 x_2 x_3 x_4$. It is therefore isomorphic to
$$\mathbb{C}[y_0,\dots,y_4,z]/(y_0 + \dots + y_4,\ z^5 = y_0 y_1 y_2 y_3 y_4),$$
via the substitutions $y_i = x_i^5$ and $z = x_0 x_1 x_2 x_3 x_4$. Thus the singular scheme is a cyclic five-fold cover of $\mathbb{P}^3$, branched along five planes.

**JA.** $\psi = 0$ のクインティック・ミラーを考える。Fermat クインティックの斉次座標環
$$\mathbb{C}[x_0,\dots,x_4]/(x_0^5 + \dots + x_4^5)$$
に作用する $G$ の**不変環**が、$x_0^5,\dots,x_4^5$ と $x_0 x_1 x_2 x_3 x_4$ で生成されることは容易にわかる。ゆえにそれは、置換 $y_i = x_i^5$, $z = x_0 x_1 x_2 x_3 x_4$ によって
$$\mathbb{C}[y_0,\dots,y_4,z]/(y_0 + \dots + y_4,\ z^5 = y_0 y_1 y_2 y_3 y_4)$$
に同型である。したがって特異スキームは、**5枚の平面に沿って分岐する $\mathbb{P}^3$ の巡回5重被覆**である。

> 〔木原メモ〕ミラー＝$(\mathbb{Z}/5)^3$ 不変環＝**有限群作用（orbifold）で作る CY**。$z^5 = y_0y_1y_2y_3y_4$ の巡回被覆＝離散内部対称性の幾何化。「Q軸＝有限群＝orbifold」と同型の構造。

---

### The Hodge diamond（ホッジ・ダイヤモンド）

**EN.** Keep in mind that for a CY3, with $\Omega^3_X = \mathcal{O}_X$, we have
$$\Omega^2_X = T_X \quad\text{and}\quad \Omega^1_X = \Lambda^2 T_X.$$
By definition, $h^{p,q} = \dim H^q(X, \Omega^p)$, the $q$-th cohomology of the sheaf of holomorphic or algebraic $p$-forms. We arrange these dimensions in a diamond ($h^{3,3}$ at top … $h^{0,0}$ at bottom).

Note that the sum of numbers in row $n$ is the dimension of $H^n(X,\mathbb{Q})$. So the alternating sum of all numbers in the Hodge diamond is the Euler characteristic of $X$:
$$\chi(X) = \sum_{p,q} (-1)^{p+q} h^{p,q}.$$
By Hodge theory, $h^{p,q} = h^{q,p}$, so the Hodge diamond is symmetric with respect to reflection across the central vertical line. By Serre duality,
$$H^q(X, \Omega^p) = H^{3-q}(X, \Lambda^p T_X)^\vee = H^{3-q}(X, \Omega^{3-p})^\vee,$$
so $h^{p,q} = h^{3-p,3-q}$. The Hodge diamond is symmetric with respect to rotation by $180°$. It follows that it is also symmetric with respect to reflection across the central horizontal line. For a strict CY3, with $h^{0,1} = h^{0,2} = 0$, we immediately get the shape of the Hodge diamond.

**JA.** CY3 では $\Omega^3_X = \mathcal{O}_X$ ゆえ、$\Omega^2_X = T_X$ かつ $\Omega^1_X = \Lambda^2 T_X$ であることに注意する。定義により $h^{p,q} = \dim H^q(X, \Omega^p)$（正則／代数的 $p$-形式の層の $q$ 次コホモロジー）。これらの次元を菱形（ダイヤモンド）に並べる（上が $h^{3,3}$、下が $h^{0,0}$）。

第 $n$ 行の数の和が $H^n(X,\mathbb{Q})$ の次元であることに注意。ゆえにダイヤモンド全体の交代和がオイラー標数：
$$\chi(X) = \sum_{p,q} (-1)^{p+q} h^{p,q}.$$
ホッジ理論により $h^{p,q} = h^{q,p}$、ゆえにダイヤモンドは**中央縦線に関して対称**。セール双対により
$$H^q(X, \Omega^p) = H^{3-q}(X, \Lambda^p T_X)^\vee = H^{3-q}(X, \Omega^{3-p})^\vee,$$
ゆえに $h^{p,q} = h^{3-p,3-q}$。ダイヤモンドは **180° 回転に関して対称**で、したがって中央横線に関しても対称。狭義 CY3（$h^{0,1} = h^{0,2} = 0$）では、ダイヤモンドの形がただちに定まる。

**狭義 CY3 のホッジ・ダイヤモンド：**
```
              1
            0   0
          0  h^{1,1}  0
       1 h^{2,1}  h^{2,1} 1
          0  h^{1,1}  0
            0   0
              1
```

**EN.** There are only two significant numbers here. The number
$$h^{2,1} = h^{1,2} = \dim H^1(X, \Omega^2) = \dim H^1(X, T_X)$$
is the dimension of the space of infinitesimal deformations of $X$ as a complex manifold. The number
$$h^{1,1} = \dim H^1(X, \Omega_X) = \dim H^1(X, \Lambda^2 T_X) = H^2(X,\mathbb{Q}) = H_2(X,\mathbb{Q})$$
is essentially the number of independent curve classes in $X$. (Also the Kähler moduli space, but we will not get into that.)

**JA.** ここで意味のある数は2つだけ。
$$h^{2,1} = h^{1,2} = \dim H^1(X, \Omega^2) = \dim H^1(X, T_X)$$
は、複素多様体としての $X$ の**無限小変形の空間の次元**（複素構造モジュライ）。
$$h^{1,1} = \dim H^1(X, \Omega_X) = \dim H^1(X, \Lambda^2 T_X) = H^2(X,\mathbb{Q}) = H_2(X,\mathbb{Q})$$
は、本質的に $X$ の**独立な曲線類の数**（ケーラー・モジュライ空間でもあるが、ここでは立ち入らない）。

**滑らかなクインティックのホッジ・ダイヤモンド**（$h^{1,1}=1$, $h^{2,1}=101$）：
```
          1
        0   0
      0   1   0
   1  101  101  1
      0   1   0
        0   0
          1
```
**解消した Dwork 商（＝ミラー）**（$h^{1,1}=101$, $h^{2,1}=1$）：
```
          1
        0   0
      0  101  0
   1   1   1   1
      0  101  0
        0   0
          1
```

**EN.** The Euler characteristic of the quintic is $6 - 202 = -196$ 〔原文ママ〕, the Euler characteristic of the mirror is $6 + 202 = 208$ 〔原文ママ〕.

**JA.** クインティックのオイラー標数は $6 - 202 = -196$〔配布冊子の記載のまま〕、ミラーのオイラー標数は $6 + 202 = 208$〔同〕。

> 〔木原・要確認メモ〕**標準的にはクインティック $\chi = -200$、ミラー $\chi = +200$**（CY3 は $\chi = 2(h^{1,1}-h^{2,1})$、クインティック $2(1-101)=-200$）。冊子の「6」は、6個の「1」エントリ（4隅 $h^{0,0},h^{3,0},h^{0,3},h^{3,3}$＋$h^{1,1},h^{2,2}$）を $+6$ と数えた値に見えるが、$h^{3,0},h^{0,3}$ は $p+q=3$ で**符号 $-$**。正しくは $+2$ で、$2-202=-200$（ミラーは $2+202$ の符号反転で $+200$）。**冊子の $-196/208$ は符号の取り違えによる誤記と思われる**。当日、Behrend 本人の板書で確認推奨。$202 = 101+101$、ミラー対称性は $h^{1,1}\leftrightarrow h^{2,1}$ を入れ替え $\chi$ の符号を反転。

---

## 2. Enumerative Geometry（数え上げ幾何）

### Why Calabi–Yau threefolds?（なぜ CY3 か）

**EN.** From the enumerative point of view, Calabi–Yau threefolds are particularly interesting because we can expect a finite number of curves of fixed degree and genus.

Consider a curve $C \subset X$, where for now $X$ is an arbitrary complex projective manifold, say of dimension $n$. We assume $C$ is non-singular of genus $g$, and assume that $X \subset \mathbb{P}^n$, so that $C$ has a degree as well, say $d$. The space of infinitesimal deformations of $C$ in $X$ is $H^0(C, N_{C/X})$. This can be approximated, using Riemann–Roch for the curve $C$, as
$$\chi(C, N_{C/X}) = \deg N_{C/X} + (n-1)(1-g),$$
with equality if $H^1(C, N_{C/X}) = 0$. From the short exact sequence
$$0 \to T_C \to T_X|_C \to N_{C/X} \to 0 \tag{2}$$
we see that $\det N_{C/X} \otimes T_C = \det T_X|_C = \omega_X^{-1}|_C$, and hence
$$\deg N_{C/X} = -\deg T_C - \deg(\omega_X|_C) = 2(g-1) - \deg(\omega_X|_C),$$
and so
$$\chi(C, N_{C/X}) = (n-1)(1-g) + 2(g-1) - \deg(\omega_X|_C) = (n-3)(1-g) - \deg(\omega_X|_C).$$
To get this to come out to zero without making special assumptions about the curve or the map, we need
- (i) $n = \dim X = 3$,
- (ii) $\omega_X$ is trivial, the Calabi–Yau condition.

Of course, the 3 is needed to cancel the $3(g-1)$ from Riemann's famous formula for the dimension of the moduli of curves of genus $g$. So in generic situations, one would expect a reasonable answer to the question:
> **How many curves of genus $g$ and degree $d$ exist in the Calabi–Yau threefold $X$?**

This is in contrast, for example, to $X = \mathbb{P}^2$, where questions like how many rational cubics pass through 8 general points have an answer (12 in this case). In the CY3 case, no fixing points or anything else is necessary.

**JA.** 数え上げの観点から、CY3 は特に興味深い。なぜなら、固定した次数と種数の曲線が**有限個**であることが期待できるからである。

曲線 $C \subset X$ を考える（当面 $X$ は任意の複素射影多様体、次元 $n$ とする）。$C$ は種数 $g$ の非特異曲線とし、$X \subset \mathbb{P}^n$ なので $C$ も次数 $d$ を持つ。$X$ 内の $C$ の無限小変形の空間は $H^0(C, N_{C/X})$。これは曲線 $C$ のリーマン・ロッホで近似でき、
$$\chi(C, N_{C/X}) = \deg N_{C/X} + (n-1)(1-g),$$
$H^1(C, N_{C/X}) = 0$ なら等号。短完全列 (2) から $\det N_{C/X} \otimes T_C = \det T_X|_C = \omega_X^{-1}|_C$、ゆえに $\deg N_{C/X} = 2(g-1) - \deg(\omega_X|_C)$、したがって
$$\chi(C, N_{C/X}) = (n-3)(1-g) - \deg(\omega_X|_C).$$
これが、曲線や写像に特別な仮定を置かずにゼロになるには、**(i) $n = \dim X = 3$、(ii) $\omega_X$ が自明（＝カラビ・ヤウ条件）**が必要。**この "3" は、種数 $g$ の曲線のモジュライ次元のリーマンの公式に現れる $3(g-1)$ を打ち消すために必要**。ゆえに generic な状況では、次の問いに合理的な答えが期待できる：
> **CY3 $X$ の中に、種数 $g$・次数 $d$ の曲線は何本あるか？**

これは例えば $X = \mathbb{P}^2$（「8点を通る有理3次曲線は何本か？」に答え 12 がある）と対照的で、CY3 では**点を固定する必要も何もない**。

> 〔木原メモ〕**"なぜ複素3次元（CY3）か"＝曲線モジュライの仮想次元が $(n-3)(1-g)$ で、$n=3$ かつ CY（$\omega_X$ 自明）のときだけゼロ＝有限の数え上げが成立する**。これは「特別な次元で数え上げ／構造が閉じる」典型。あなたの「整数振動数の等方基底が閉じる最小非自明次元＝4」と精神は同型（次元が counting を決める）。ただし CY は複素3次元、あなたは実4次元で別物。

### Clemens conjecture（クレメンス予想）

**EN.** Let $X$ be a sufficiently general/generic quintic threefold. The number of rational curves ($g=0$) of fixed degree is finite. For example the number of lines is
$$n_1 = 2875.$$
The number of conics is
$$n_2 = 609\,250.$$
The conjecture has only been proved for $d \le 11$.

These lines and conics all look like $\mathbb{P}^1 \hookrightarrow \mathcal{O}(-1) \oplus \mathcal{O}(-1)$, if $X$ is generic enough. (There is even an analytic tubular neighbourhood of $\mathbb{P}^1$ in $X$, which is complex analytically isomorphic to an analytic open neighbourhood of the zero section in $\mathcal{O}(-1) \oplus \mathcal{O}(-1)$.) Ideally, all rational curves look like this, but this (the strong version of Clemens conjecture) is false.

In general, the normal bundle of $C \to X$ is a vector bundle of rank 2, $N \to C$, such that $\det N = \omega_C$. For example $N = L_1 \oplus L_2$, such that $L_1 \otimes L_2 \xrightarrow{\sim} \omega_C$. This kind of open CY3 is called a *local curve*, and is very useful for heuristic study of curves in CY3s.

**JA.** $X$ を十分一般的（generic）なクインティック3次元体とする。固定次数の有理曲線（$g=0$）の数は有限。例えば**直線の数は $n_1 = 2875$**、**円錐曲線（2次）の数は $n_2 = 609\,250$**。この予想は $d \le 11$ でのみ証明されている。

これら直線・2次曲線は、$X$ が十分 generic なら、すべて $\mathbb{P}^1 \hookrightarrow \mathcal{O}(-1) \oplus \mathcal{O}(-1)$ のように見える（$X$ 内の $\mathbb{P}^1$ の解析的管状近傍が、$\mathcal{O}(-1) \oplus \mathcal{O}(-1)$ の零切断の近傍と複素解析的に同型でさえある）。理想的には全有理曲線がこう見えてほしいが、これ（クレメンス予想の強い版）は**偽**。

一般には、$C \to X$ の法束は階数2のベクトル束 $N \to C$ で $\det N = \omega_C$。例えば $N = L_1 \oplus L_2$ で $L_1 \otimes L_2 \cong \omega_C$。この種の開 CY3 を **局所曲線（local curve）** と呼び、CY3 内の曲線の発見的研究に非常に有用。

> 〔木原メモ〕**2875（直線）・609250（2次）** ＝有名な数え上げ数。法束が $\mathcal{O}(-1)\oplus\mathcal{O}(-1)$（"$(-1,-1)$ 曲線"）＝局所モデル。後の DT/位相的頂点は、この局所構造の数え上げに帰着。

### The story of degree $d$ rational curves on a quintic（Candelas ら 1991）

**EN.** The following was discovered by Candelas, de la Ossa, Green, and Parkes in 1991. Let $n_d$, $d = 1, 2, \dots$ be the number of rational curves of degree $d$ on the generic quintic 3-fold $X$. Form the generating function
$$\partial_t F(t) = \frac{5}{2}t^2 + \frac{25}{6}t + \frac{25}{6} + \frac{1}{(2\pi i)^2} \sum_{d=1}^{\infty} d\, n_d\, \mathrm{Li}_2\!\left(e^{2\pi i d t}\right),$$
and consider the two periods:
$$y_0(z) = \sum_{n=0}^{\infty} \frac{(5n)!}{(n!)^5} z^n,$$
$$y_1(z) = y_0(z)\log(z) + 5 \sum_{n=1}^{\infty} \frac{(5n)!}{(n!)^5}\left(\sum_{j=n+1}^{5n}\frac{1}{j}\right) z^n.$$
Form the mirror map $t(z) = y_1(z)/y_0(z)$; then the above generating function after changing coordinates from $t$ to $z$ using this quotient of periods, is a solution to the Picard–Fuchs equation
$$\left[\left(z\frac{d}{dz}\right)^4 - 5^5 z\left(z\frac{d}{dz}+\tfrac15\right)\left(z\frac{d}{dz}+\tfrac25\right)\left(z\frac{d}{dz}+\tfrac35\right)\left(z\frac{d}{dz}+\tfrac45\right)\right]\left(y_0(z)\,\partial_t F(t(z))\right) = 0.$$
This determines the $n_d$ recursively.

What goes on here? There is the mirror quintic $\hat{X}$, which is the resolved Dwork pencil quotient, and the enumerative geometry of $X$ is related to the Hodge theory (periods) of $\hat{X}$. Explaining this would take us too far afield; we just wanted to give a hint at how complicated the enumerative geometry of CY3s turns out to be.

**JA.** 以下は **Candelas, de la Ossa, Green, Parkes が 1991 年に発見**した。$n_d$（$d=1,2,\dots$）を generic クインティック3次元体 $X$ 上の次数 $d$ の有理曲線の数とする。母関数
$$\partial_t F(t) = \frac{5}{2}t^2 + \frac{25}{6}t + \frac{25}{6} + \frac{1}{(2\pi i)^2} \sum_{d=1}^{\infty} d\, n_d\, \mathrm{Li}_2\!\left(e^{2\pi i d t}\right)$$
を作り、2つの周期
$$y_0(z) = \sum_{n=0}^{\infty} \frac{(5n)!}{(n!)^5} z^n,\qquad y_1(z) = y_0(z)\log(z) + 5 \sum_{n=1}^{\infty} \frac{(5n)!}{(n!)^5}\left(\sum_{j=n+1}^{5n}\frac{1}{j}\right) z^n$$
を考える。**ミラー写像** $t(z) = y_1(z)/y_0(z)$ を作ると、上の母関数は、この周期の商で $t$ から $z$ へ座標変換した後、**ピカール・フックス方程式**
$$\left[\left(z\tfrac{d}{dz}\right)^4 - 5^5 z\textstyle\prod_{k=1}^{4}\left(z\tfrac{d}{dz}+\tfrac{k}{5}\right)\right]\left(y_0(z)\,\partial_t F(t(z))\right) = 0$$
の解になる。これが $n_d$ を**再帰的に**決定する。

何が起きているか。**ミラー・クインティック $\hat{X}$**（＝解消した Dwork 鉛筆束の商）があり、$X$ の数え上げ幾何が $\hat{X}$ のホッジ理論（周期）と関係する。これを説明すると話が逸れすぎるので、CY3 の数え上げ幾何がいかに複雑かのヒントを与えるに留める。

> 〔木原メモ〕**ミラー対称性の核心**：A 模型（$X$ の曲線数え上げ $n_d$）＝ B 模型（$\hat{X}$ の周期＝ピカール・フックス）。周期 $y_0=\sum\frac{(5n)!}{(n!)^5}z^n$（超幾何級数）。**$\hat{X}$＝$(\mathbb{Z}/5)^3$ orbifold（前出）**。「数え上げ ↔ 周期（解析）」の双対は、あなたの「外部観測＝内部和」「離散 ↔ 連続」の双対と精神的に共鳴（ただし別物・厳密対応ではない）。

### How to make sense of the numbers（数の意味づけ — Fermat の例）

**EN.** As we had seen, even for the single most prominent example of a CY3 (the quintic), the finiteness of $n_{g,d}$ is only known for $g=0$ when $d \le 11$. So it is a success to even prove that there is a well-defined number $n_{g,d}$ that one might have a chance to compute. (In other words, what did Candelas et al. actually compute?)

Let us explain this kind of phenomenon for the case of lines in a non-generic quintic. Probably, for higher genus or degree all quintics are non-generic, so this is most likely a typical situation. Consider the Fermat quintic
$$x_0^5 + x_1^5 + x_2^5 + x_3^5 + x_4^5 = 0.$$
It is not generic with respect to lines. In fact, it has **50 pencils of lines**: Given $\langle a,b,c\rangle \subset \mathbb{P}^2$ lying on the Fermat curve $a^5 + b^5 + c^5 = 0$, we have
$$\mathbb{P}^1 \to \mathbb{P}^4,\qquad \langle u, v\rangle \mapsto \langle u, \zeta u, a v, b v, c v\rangle.$$
We have 10 choices of where to put the $u$ vs. the $v$ part, and 5 choices for the root of $-1$, $\zeta$. There are **375 lines** which are on two of these components; they are given by setting one of $a, b, c$ equal to zero. The graph looks like this: 〔図〕

**JA.** 既に見たように、CY3 の最も顕著な唯一の例（クインティック）でさえ、$n_{g,d}$ の有限性は $g=0$ かつ $d \le 11$ でしか知られていない。だから、計算できる見込みのある well-defined な数 $n_{g,d}$ が存在することを証明するだけでも成功なのである。（言い換えれば、Candelas らは実際に何を計算したのか？）

この種の現象を、**non-generic なクインティック上の直線**の場合で説明しよう。おそらく高種数・高次数では全クインティックが non-generic なので、これが典型的状況だろう。Fermat クインティック
$$x_0^5 + x_1^5 + x_2^5 + x_3^5 + x_4^5 = 0$$
を考える。これは直線に関して generic でない。実際、**50 個の直線の鉛筆束（pencils）**を持つ：Fermat 曲線 $a^5 + b^5 + c^5 = 0$ 上の $\langle a,b,c\rangle \subset \mathbb{P}^2$ に対し、
$$\mathbb{P}^1 \to \mathbb{P}^4,\qquad \langle u, v\rangle \mapsto \langle u, \zeta u, a v, b v, c v\rangle.$$
$u$ と $v$ の置き場所の選び方が **10 通り**、$-1$ の根 $\zeta$ の選び方が **5 通り**（計 50）。これらの成分のうち2つに乗る**直線は 375 本**で、$a,b,c$ のいずれかを 0 とおくことで与えられる。グラフは次のようになる：〔図〕

> 〔木原メモ〕**Fermat（最も対称な特異クインティック）では、有限個（2875）であるべき直線が "50 鉛筆束（連続族）" に退化**＝generic で有限な数え上げが、対称点では連続族に化ける（後の §3-6 で「仮想的に数える＝Behrend 函数・DT」で救済）。**$50 = 10\times5$（配置10×5乗根）、$375$（2成分に乗る線）** ＝対称配置上の純組合せ計数。$\zeta$＝5乗根＝$(\mathbb{Z}/5)$ 構造。あなたの「対称（等方）基底では退化／許容状態の離散選択」と同じ匂い（対称性が数え上げを退化させる）。

---

### モジュライ空間の二つの記述（Grassmann 束と仮想計数 2875）

**EN.** The moduli space of lines on the Fermat quintic consists of 50 genus 6 smooth curves, which intersect in 375 ordinary double points, each component intersects 15 others. Away from the intersections this moduli space is not smooth, but carries a scheme structure transversally given by $\mathbb{C}[t]/t^2$; in other words the scheme structure is the critical locus of $t^3$. At the nodes, the scheme structure is the critical locus of $t_1^2 t_2^2$.

This is quite complicated, but we can describe it a different way, namely as the space of lines in $\mathbb{P}^4(\mathbb{C})$, which happen to lie on $X \subset \mathbb{P}^4(\mathbb{C})$. Let $Y$ be the space of lines in $\mathbb{P}^4(\mathbb{C})$. (This is a Grassmannian variety, it has dimension 6.) Over $Y$ we have a vector bundle whose fibre over a line $L \subset \mathbb{P}^4(\mathbb{C})$, given by a 2-dimensional subspace $W \subset \mathbb{C}^5$, is $\mathrm{Sym}^5 W^*$. Call this vector bundle $\mathrm{Sym}^5 W^*$ over $Y$ (it has rank 6). The Fermat equation to $W$ gives us a global section $s \in \Gamma(Y, \mathrm{Sym}^5 W^*)$. The zero locus of this section is the space of lines on the Fermat quintic
$$M_{0,1}(X) = Z(s) \subset Y.$$

For a generic quintic this section will have a discrete set of zeros, for the Fermat quintic, it has the complicated structure described above. The *virtual fundamental class* of $M_{0,1}(X)$ is the localized top Chern class of the bundle. The point is that the rank of the bundle is equal to the dimension of $Y$. It is a homology class of degree 0, which pushes forward into $Y$ to give the virtual Poincaré dual of the top Chern class of the bundle $\mathrm{Sym}^5 W^*$. The virtual number of points on the Fermat quintic is
$$\int_{[M_{0,1}(X)]^{\mathrm{virt}}} 1 = \deg e(\mathrm{Sym}^5 W^*) \cap [Y] = 2875.$$

Main advances: Fulton-MacPherson style intersection theory, upgraded in the 1990s by Li-Tian and Behrend-Fantechi to be intrinsic to the moduli space $M_{0,1}(X)$ without embedding into an ambient $Y$, which may not exist, or maybe useless for the definition of a virtual class.

**JA.** Fermat クインティック上の直線のモジュライ空間は、**種数 6 の滑らかな曲線 50 本**からなり、それらは **375 個の通常二重点（ODP）**で交わる。各成分は他の 15 本と交わる。交点を離れたところではこのモジュライ空間は滑らかでなく、横断的に $\mathbb{C}[t]/t^2$ で与えられるスキーム構造を持つ。言い換えればスキーム構造は $t^3$ の臨界軌跡である。節点では、スキーム構造は $t_1^2 t_2^2$ の臨界軌跡である。

これは非常に複雑だが、別の方法で記述できる。すなわち $\mathbb{P}^4(\mathbb{C})$ 内の直線の空間のうち、たまたま $X \subset \mathbb{P}^4(\mathbb{C})$ 上に乗るもの、として。$Y$ を $\mathbb{P}^4(\mathbb{C})$ 内の直線の空間とする。（これは**グラスマン多様体**で、次元 6 を持つ。）$Y$ 上にベクトル束があり、その直線 $L \subset \mathbb{P}^4(\mathbb{C})$（2 次元部分空間 $W \subset \mathbb{C}^5$ で与えられる）上のファイバーは $\mathrm{Sym}^5 W^*$ である。このベクトル束を $Y$ 上の $\mathrm{Sym}^5 W^*$ と呼ぶ（階数 6）。Fermat 方程式は $W$ に対し大域切断 $s \in \Gamma(Y, \mathrm{Sym}^5 W^*)$ を与える。この切断の零点軌跡が Fermat クインティック上の直線の空間
$$M_{0,1}(X) = Z(s) \subset Y$$
である。

generic なクインティックではこの切断は離散的な零点集合を持つが、Fermat クインティックでは上述の複雑な構造を持つ。$M_{0,1}(X)$ の**仮想基本類（virtual fundamental class）**は束の局所化された最高次チャーン類である。要点は、**束の階数が $Y$ の次元に等しい**ことである。それは次数 0 のホモロジー類で、$Y$ へ押し出すと束 $\mathrm{Sym}^5 W^*$ の最高次チャーン類の仮想ポアンカレ双対を与える。Fermat クインティック上の点の仮想個数は
$$\int_{[M_{0,1}(X)]^{\mathrm{virt}}} 1 = \deg e(\mathrm{Sym}^5 W^*) \cap [Y] = 2875.$$

主要な進展：Fulton-MacPherson 流の交叉理論が、1990 年代に Li-Tian と Behrend-Fantechi によって、（存在しないか、仮想類の定義に役立たないかもしれない）周囲空間 $Y$ への埋め込みを経ずに、モジュライ空間 $M_{0,1}(X)$ に内在的（intrinsic）なものへと格上げされた。

> 〔木原メモ〕**同じ数 2875 を二通りに数える**：(a) 各直線を点として数える（generic）／(b) Grassmann 束 $Y$（dim 6）上の階数 6 の束 $\mathrm{Sym}^5 W^*$ の最高次チャーン類 $\int_Y e(\mathrm{Sym}^5 W^*)$。**「rank = dim」だから 0 次元の仮想点が出る**＝過剰交叉を Euler 類で救う構図。退化して連続族（50 鉛筆束）になっても、仮想基本類は依然 2875 を返す。**束の階数とベース次元の一致＝あなたの「自由度の数合わせ（許容状態の離散選択）」と同型の発想**（個数は配置の組合せでなく特性類で決まる）。

---

## 3　シンプレクティック幾何（Symplectic geometry）

**EN.** Consider again the tangent space to the space of curves in $X$ at a given curve $C \to X$. This is $H^0(C, N_{C/X})$, normal vector fields to $C$ in $X$. Now $N_{C/X}$ is of rank two. The existence of $H^1(C, N_{C/X})$ is the main obstruction to the moduli space being of the expected dimension, i.e., just a finite set of points.

Consider again the short exact sequence (2). It implies
$$\det N_{C/X} = \omega_C,$$
if $X$ is Calabi-Yau. Since $N_{C/X}$ is of rank two, we have
$$N^\vee_{C/X} \otimes \det N_{C/X} = N_{C/X}.$$
and hence
$$N^\vee_{C/X} \otimes \omega_C = N_{C/X}.$$
Consider Serre-duality for the curve $C$:
$$H^1(C, N_{C/X})^\vee = H^0(C, N^\vee_{C/X} \otimes \omega_C) = H^0(C, N_{C/X}).$$

**The tangent space to the space of curves in $X$ is the dual space of the obstruction space.**

**JA.** 与えられた曲線 $C \to X$ における、$X$ 内の曲線の空間への接空間を再び考える。これは $H^0(C, N_{C/X})$、すなわち $X$ 内の $C$ への法線ベクトル場である。さて $N_{C/X}$ は階数 2 である。$H^1(C, N_{C/X})$ の存在こそが、モジュライ空間が期待次元（＝単なる有限個の点）になることへの主要な障害（obstruction）である。

短完全列 (2) を再び考えると、$X$ がカラビ・ヤウならば
$$\det N_{C/X} = \omega_C$$
が従う。$N_{C/X}$ は階数 2 なので、
$$N^\vee_{C/X} \otimes \det N_{C/X} = N_{C/X},$$
したがって
$$N^\vee_{C/X} \otimes \omega_C = N_{C/X}.$$
曲線 $C$ のセール双対性を考えると：
$$H^1(C, N_{C/X})^\vee = H^0(C, N^\vee_{C/X} \otimes \omega_C) = H^0(C, N_{C/X}).$$

**$X$ 内の曲線の空間への接空間は、障害空間の双対空間である。**

> 〔木原メモ〕**接空間 $T^0=H^0$ と障害空間 $T^1=H^1$ が CY 条件で互いに双対**になる（$\det N=\omega_C$ → Serre 双対）。これが「$-1$-シフト・シンプレクティック構造」の芽。**自由度（接）と拘束（障害）が対をなして釣り合う**＝あなたの「内部自由度と拘束の対称な対合」と同じ骨格（CY 条件が双対性を強制する）。

**EN.** This exhibits a much deeper symmetry of the space of curves of geometry of curves being zero. In fact, the expected dimension of the space of curves in a CY3 is zero. The obstructions should be considered, with a symplectic form of degree $-1$, and so we get an isomorphism of the tangent vectors in degree 1, and the (shifted) dual. This turns out to be a symplectic complex of $M$ with the (shifted) dual. The complex is (skew-symmetric and closed).

$$
\begin{array}{ccc}
T^0_M & \longrightarrow & T^1_M \\
\downarrow & & \downarrow \\
(T^1_M)^\vee & \longrightarrow & (T^0_M)^\vee
\end{array}
\qquad
\begin{array}{c}
T^\bullet_M \\
\downarrow \\
(T^\bullet_M)^\vee \\
\\
(T^\bullet_M)^\vee[-1]
\end{array}
$$

**Theorem 3.1 (shifted Darboux theorem)** *The local geometry of a $-1$-shifted symplectic space with $-1$-shifted symplectic structure is given by the (derived) scheme theoretic critical locus of a holomorphic or polynomial function.*

Given a function $f : \mathbb{C}^n \to \mathbb{C}$, the critical locus $\mathrm{Crit}(f)$ is defined as the common zero locus of the $n$ functions
$$\frac{\partial f}{\partial x_1}, \dots, \frac{\partial f}{\partial x_n} : \mathbb{C}^n \to \mathbb{C}.$$

The number of equations being equal to the number of variables already indicates that the expected dimension is zero. More precisely, the tangent space and the obstruction space at the origin are, respectively, the kernel and the cokernel of the matrix of partial derivatives of these functions namely;
$$0 \to T^0_{\mathrm{Crit}\,f}|_0 \to \mathbb{C}^n \xrightarrow{\frac{\partial^2 f}{\partial x_i \partial x_j}|_0} \mathbb{C}^n \to T^1_{\mathrm{Crit}\,f}|_0 \to 0$$
The key is that $T^0$ (tangents) and $T^1$ (obstructions) are duals of each other. This proves that the Hessian matrix of second partials is symmetric. This is the main observation in proving that a critical locus is a $-1$-shifted symplectic (derived) scheme.

A large part of research in this area in recent years has focused on finding functions whose critical loci describe moduli spaces of curves in CY3s, at least locally.

**JA.** これは曲線の空間の、はるかに深い対称性を示す——CY3 内の曲線の空間の期待次元はゼロである。障害は次数 $-1$ のシンプレクティック形式とともに考えるべきで、その結果、次数 1 の接ベクトルと（シフトした）双対との同型が得られる。これは $M$ の複体が（シフトした）双対と一致する**シンプレクティック複体**であることを意味する（歪対称かつ閉）。

（上図：$T^0_M \to T^1_M$ と双対 $(T^1_M)^\vee \to (T^0_M)^\vee$ が対をなし、複体 $T^\bullet_M$ がその双対 $(T^\bullet_M)^\vee[-1]$ と同型。）

**定理 3.1（シフト・ダルブー定理）** *$-1$-シフト・シンプレクティック構造をもつ $-1$-シフト・シンプレクティック空間の局所幾何は、正則関数または多項式関数の（導来）スキーム論的臨界軌跡で与えられる。*

関数 $f : \mathbb{C}^n \to \mathbb{C}$ に対し、臨界軌跡 $\mathrm{Crit}(f)$ は $n$ 個の関数
$$\frac{\partial f}{\partial x_1}, \dots, \frac{\partial f}{\partial x_n} : \mathbb{C}^n \to \mathbb{C}$$
の共通零点軌跡として定義される。

**方程式の本数が変数の本数に等しい**ことが、既に期待次元ゼロを示している。より正確には、原点における接空間と障害空間は、それぞれこれらの関数の偏導関数行列の核と余核であり、
$$0 \to T^0_{\mathrm{Crit}\,f}|_0 \to \mathbb{C}^n \xrightarrow{\frac{\partial^2 f}{\partial x_i \partial x_j}|_0} \mathbb{C}^n \to T^1_{\mathrm{Crit}\,f}|_0 \to 0.$$
要点は、$T^0$（接）と $T^1$（障害）が互いに双対であることである。これは**2 階偏微分のヘッセ行列が対称である**ことを証明する。これが、臨界軌跡が $-1$-シフト・シンプレクティック（導来）スキームであることを証明する際の主要な観察である。

近年この分野の研究の大部分は、CY3 内の曲線のモジュライ空間を（少なくとも局所的に）記述するような、臨界軌跡をもつ関数を見つけることに注がれてきた。

> 〔木原メモ〕**全部が一つの「ポテンシャル関数 $f$ の臨界点 $\mathrm{Crit}(f)=\{\nabla f=0\}$」に帰着する**（シフト・ダルブー）。**方程式数＝変数数**で期待次元 0、**ヘッセ行列の対称性＝接と障害の双対性**。これはあなたの離散基礎の「許容状態＝あるポテンシャル／作用の停留点」「対称な 2 階構造（ヘッセ）」と直接同型の発想。**数え上げ＝停留点の重み付き計数**という構図が以降ずっと効く。

---

## 4　臨界軌跡（Critical loci）

**EN.** Let $M$ be a smooth, non-compact complex manifold, and let $f : M \to \mathbb{C}$ be a holomorphic function. We consider the critical locus $X = \mathrm{Crit}(f) \subset M$. For the purposes of this discussion, let us assume that $X$ is compact.

The space $X$ can be realized as the intersection of two Lagrangian submanifolds (of complementary dimensions) within the cotangent bundle $\Omega_M$ (which is a complex symplectic manifold). Specifically, it is the intersection of the zero section and the graph of the differential $df$:
$$
\begin{array}{ccc}
X & \longrightarrow & M \\
\downarrow & & \downarrow \Gamma_{df} \\
M & \xrightarrow{\ 0\ } & \Omega_M
\end{array}
$$

Because $X$ is compact, the virtual intersection number is well-defined:
$$\#^{\mathrm{virt}}(X) = \mathcal{I}_{\Omega_M}(M, \Gamma_{df}) = \int_{[X]^{\mathrm{virt}}} 1$$
where $[X]^{\mathrm{virt}} \in H_0(X, \mathbb{Q})$ is the virtual fundamental class of the intersection scheme $X$.

**Theorem 4.1 (Singular Gauß-Bonnet)** *The intersection number is given by the Euler characteristic of $X$ weighted by a constructible function $\mu : X \to \mathbb{Z}$, namely:*
$$\mathcal{I}_{\Omega_M}(M, \Gamma_{df}) = \chi(X, \mu).$$
*For a point $P \in X$, the value $\mu(P)$ is the Milnor number of $f$ at $P$:*
$$\mu(P) = (-1)^{\dim M}(1 - \chi(F_P)).$$
*Here $F_P$ is the Milnor fibre of $f$ at the point $P$, which is the intersection of a nearby fibre of $f$ with a small ball around $P$.*

**JA.** $M$ を滑らかな非コンパクト複素多様体とし、$f : M \to \mathbb{C}$ を正則関数とする。臨界軌跡 $X = \mathrm{Crit}(f) \subset M$ を考える。この議論のため、$X$ はコンパクトと仮定する。

空間 $X$ は、複素シンプレクティック多様体である余接束 $\Omega_M$ の内部で、（相補次元の）二つのラグランジュ部分多様体の交わりとして実現できる。具体的には、零切断と微分 $df$ のグラフの交わりである：
$$
\begin{array}{ccc}
X & \longrightarrow & M \\
\downarrow & & \downarrow \Gamma_{df} \\
M & \xrightarrow{\ 0\ } & \Omega_M
\end{array}
$$

$X$ がコンパクトなので、**仮想交叉数**が well-defined である：
$$\#^{\mathrm{virt}}(X) = \mathcal{I}_{\Omega_M}(M, \Gamma_{df}) = \int_{[X]^{\mathrm{virt}}} 1,$$
ここで $[X]^{\mathrm{virt}} \in H_0(X, \mathbb{Q})$ は交叉スキーム $X$ の仮想基本類である。

**定理 4.1（特異ガウス・ボンネ）** *交叉数は、構成可能関数 $\mu : X \to \mathbb{Z}$ で重み付けされた $X$ のオイラー標数で与えられる：*
$$\mathcal{I}_{\Omega_M}(M, \Gamma_{df}) = \chi(X, \mu).$$
*点 $P \in X$ に対し、値 $\mu(P)$ は $P$ における $f$ のミルナー数である：*
$$\mu(P) = (-1)^{\dim M}(1 - \chi(F_P)).$$
*ここで $F_P$ は点 $P$ における $f$ のミルナーファイバーで、$f$ の近傍ファイバーと $P$ のまわりの小球との交わりである。*

> 〔木原メモ〕**$X=\mathrm{Crit}(f)$ ＝ 余接束 $\Omega_M$ 内の「零切断 ∩ $df$ のグラフ」というラグランジアン交叉**。数え上げ $\#^{\mathrm{virt}}(X)$ が、各点の局所データ（ミルナー数 $\mu(P)$）を**重みにした Euler 標数 $\chi(X,\mu)$** に等しい（特異ガウス・ボンネ）。$\mu(P)=(-1)^{\dim M}(1-\chi(F_P))$ は**純粋に局所・組合せ的な整数**。あなたの「大域的な数え上げ＝局所重みの離散和」の理念と同じ。$\mathbf{Behrend\ 函数\ \mu}$ がこの講演の主役。

### 孤立特異点と例 $f=x^2+y^2$

**EN.** This singular Gauß-Bonnet theorem is a result of microlocal geometry, from the 1970s. Main ingredient in the proof is the microlocal index theorem of Kashiwara and MacPherson. Also involved is the determination of the characteristic variety of the perverse sheaf of vanishing cycles.

**Isolated Singularities (dim $X = 0$).** In the case where the critical locus consists of isolated singularities, the surface $f^{-1}(0)$ near $P$ forms a cone over the link of the singularity. The link bounds the Milnor fibre, which supports the vanishing cycles.

**Example:** $f(x,y) = x^2 + y^2$.

For the function $f(x,y) = x^2 + y^2$, the Milnor number is $\mu(P) = 1$. The critical locus is $X = \mathrm{Crit}(f) = \{P\}$, which is an isolated singularity. It is defined by setting the two partials $2x = 0$ and $2y = 0$, so it is the origin with multiplicity 1. The intersection number is 1.

Near the point $P$, the surface $f^{-1}(0)$ forms a cone over the link of the singularity. (The link is the intersection of $f^{-1}(0)$ with a small ball around $P$.) The link consists of two circles, the cone over the link looks like a classical cone. The cone is contractible. The Milnor fibre $F_P$ is a manifold with boundary, and that boundary is the link itself. It is a hyperboloid. It has Euler characteristic 0. The Milnor number is therefore
$$\mu(P) = (-1)^{\dim M}(1 - \chi(F_P)) = 1,$$
which agrees with the intersection number.

**JA.** この特異ガウス・ボンネ定理は、1970 年代のマイクロローカル幾何の結果である。証明の主要な材料は柏原・MacPherson のマイクロローカル指数定理である。また、消滅サイクルの偏屈層（perverse sheaf）の特性多様体の決定も関わる。

**孤立特異点（$\dim X = 0$）。** 臨界軌跡が孤立特異点からなる場合、$P$ の近くの曲面 $f^{-1}(0)$ は特異点のリンク上の錐をなす。リンクはミルナーファイバーを縁取り、ミルナーファイバーが消滅サイクルを担う。

**例：** $f(x,y) = x^2 + y^2$。

関数 $f(x,y) = x^2 + y^2$ に対し、ミルナー数は $\mu(P) = 1$。臨界軌跡は $X = \mathrm{Crit}(f) = \{P\}$ で、孤立特異点である。これは二つの偏微分 $2x = 0$, $2y = 0$ を満たす点として定義され、重複度 1 の原点である。交叉数は 1。

点 $P$ の近くで、曲面 $f^{-1}(0)$ は特異点のリンク上の錐をなす。（リンクは $f^{-1}(0)$ と $P$ まわりの小球の交わり。）リンクは**二つの円**からなり、その上の錐は古典的な錐に見える。錐は可縮である。ミルナーファイバー $F_P$ は境界をもつ多様体で、その境界はリンク自身である。それは**双曲面（hyperboloid）**で、オイラー標数 0 を持つ。したがってミルナー数は
$$\mu(P) = (-1)^{\dim M}(1 - \chi(F_P)) = 1$$
となり、交叉数と一致する。

> 〔木原メモ〕最も単純な孤立特異点 $x^2+y^2$：リンク＝2 円、ミルナーファイバー＝双曲面（$\chi=0$）、よって $\mu=(-1)^2(1-0)=1$。**$\mu$ が「消滅サイクルの数」を数える整数**である最初の具体例。次の $x^2+y^3$ で $\mu=2$ になる（非自明）。

---

### 例 $f=x^2+y^3$（三葉結び目・$\mu=2$）

**EN.** **Example:** $f(x,y) = x^2 + y^3$.

This time, the partial derivatives are $2x = 0$, and $3y^2 = 0$, so this is an isolated point of multiplicity 2. The intersection number is 2.

The link is a trefoil knot. The singularity of $f^{-1}(0)$ at the origin $P$ is the cone over the trefoil knot. The Milnor fibre is again a surface with boundary, but this time the boundary is 1 circle, in the shape of a trefoil knot. (In the sketch, the interior of the sphere is visible in the central triangle of the trefoil.) The homotopy type of the Milnor fibre is that of a sphere with three points removed, which is the same as a plane with two points removed, or a bouquet of two circles or a figure eight. The Euler characteristic of the Milnor fibre is $-1$, and the Milnor number at $P$ is $\mu(P) = (-1)^2(1 - (-1)) = 2$, indicating that there are 2 vanishing cycles.

Again the Milnor number reproduces the intersection number.

In fact, for isolated singularities, $\dim X = 0$, it is a classical theorem of Milnor that the Milnor fibre has the homotopy type of a bouquet of spheres. The multiplicity of the singularity equals the number of these spheres. In the case where the actual dimension of the critical locus equals the expected dimension, the singular Gauß-Bonnet theorem is Milnor's theorem.

**JA.** **例：** $f(x,y) = x^2 + y^3$。

今回は偏微分が $2x = 0$, $3y^2 = 0$ なので、これは**重複度 2** の孤立点である。交叉数は 2。

リンクは**三葉結び目（trefoil knot）**である。原点 $P$ における $f^{-1}(0)$ の特異点は三葉結び目上の錐である。ミルナーファイバーは再び境界をもつ曲面だが、今回は境界が三葉結び目の形をした **1 個の円**である。（図では、球の内部が三葉の中央三角形に見えている。）ミルナーファイバーのホモトピー型は**3 点を除いた球面**で、これは 2 点を除いた平面、あるいは**二つの円のブーケ（八の字）**と同じである。ミルナーファイバーのオイラー標数は $-1$ で、$P$ におけるミルナー数は $\mu(P) = (-1)^2(1 - (-1)) = 2$、すなわち**消滅サイクルが 2 個**あることを示す。

再びミルナー数が交叉数を再現する。

実際、孤立特異点（$\dim X = 0$）に対しては、ミルナーファイバーが**球面のブーケ**のホモトピー型を持つことがミルナーの古典的定理である。特異点の重複度がこれらの球面の個数に等しい。臨界軌跡の実次元が期待次元に等しい場合、特異ガウス・ボンネ定理はミルナーの定理に帰着する。

> 〔木原メモ〕$x^2+y^3$：リンク＝三葉結び目、ミルナーファイバー＝八の字（$\chi=-1$）、$\mu=(-1)^2(1-(-1))=2$。**$\mu$＝消滅サイクルの個数＝特異点の重複度**。Milnor の定理「孤立特異点 → 球面のブーケ、重複度＝球面の数」。あなたの「離散的な巻き数・サイクルの数で不変量が決まる」と同じ語法（連続的な曲面が、結び目・球面の個数という整数に還元される）。

### 滑らかな場合（Smooth Case）と定理 4.2・4.3

**EN.** **The Smooth Case.** Suppose instead that $X$ is smooth. The various obstruction spaces form a vector bundle over $X$. Since at each point of $X$ the obstruction space is the dual space to the tangent space, the obstruction bundle is also the dual of the tangent bundle. It follows that
$$\mathcal{I}_{\Omega_M}(M, \Gamma_{df}) = \int_{[X]^{\mathrm{virt}}} 1 = \int_{[X]} c_{\mathrm{top}}(\Omega_X) = (-1)^{\dim X} \int_{[X]} c_{\mathrm{top}}(T_X) = (-1)^{\dim X} \chi(X),$$
by the classical Gauß-Bonnet theorem.

On the other hand, for the critical locus of a function to be non-singular, the function has to essentially look like:
$$f(x_1, \dots, x_n) = x_1^2 + \dots + x_k^2.$$
Let's change coordinates, assuming $k$ is even
$$f(x_1, \dots, x_n) = x_1 x_2 + \dots + x_{k-1} x_k.$$
The Milnor fibre is
$$F_0 = \{(x; y) \in \mathbb{C}^k \times \mathbb{C}^{n-k} \mid x_1 x_2 + \dots + x_{k-1} x_k = \delta,\ \|x\|^2 + \|y\|^2 \le \epsilon^2\}.$$
It is acted upon by $S^1$:
$$\theta \cdot (x; y) = (\theta x_1, \theta^{-1} x_2, \dots, \theta x_{k-1}, \theta^{-1} x_k, \theta y_1, \dots, \theta y_n).$$
This action of $S^1$ on $F_0$ is free, forcing the Euler characteristic of $F_0$ to be 0. We see that
$$\mu(P) = (-1)^{\dim M} = (-1)^{\dim X}.$$
The same result holds when $k$ is odd.

**Theorem 4.2** *The constructible function $\mu : X \to \mathbb{Z}$ is intrinsic to the scheme $X$, if $X$ is smooth it takes the value $\mu(P) = (-1)^{\dim X}$, for all $P \in X$.*

We see that
$$\mathcal{I}_{\Omega_M}(M, \Gamma_{df}) = \int_{[X]^{\mathrm{vir}}} 1 = (-1)^{\dim X} \chi(X) = \chi(X, \mu).$$

**JA.** **滑らかな場合。** 代わりに $X$ が滑らかと仮定する。様々な障害空間は $X$ 上のベクトル束をなす。$X$ の各点で障害空間は接空間の双対空間なので、障害束も接束の双対である。したがって、古典的ガウス・ボンネ定理により
$$\mathcal{I}_{\Omega_M}(M, \Gamma_{df}) = \int_{[X]^{\mathrm{virt}}} 1 = \int_{[X]} c_{\mathrm{top}}(\Omega_X) = (-1)^{\dim X} \int_{[X]} c_{\mathrm{top}}(T_X) = (-1)^{\dim X} \chi(X).$$

他方、関数の臨界軌跡が非特異であるためには、関数は本質的に
$$f(x_1, \dots, x_n) = x_1^2 + \dots + x_k^2$$
の形でなければならない。$k$ を偶数と仮定して座標変換すると
$$f(x_1, \dots, x_n) = x_1 x_2 + \dots + x_{k-1} x_k.$$
ミルナーファイバーは
$$F_0 = \{(x; y) \in \mathbb{C}^k \times \mathbb{C}^{n-k} \mid x_1 x_2 + \dots + x_{k-1} x_k = \delta,\ \|x\|^2 + \|y\|^2 \le \epsilon^2\}.$$
これには $S^1$ が作用する：
$$\theta \cdot (x; y) = (\theta x_1, \theta^{-1} x_2, \dots, \theta x_{k-1}, \theta^{-1} x_k, \theta y_1, \dots, \theta y_n).$$
この $F_0$ への $S^1$ 作用は**自由**で、$F_0$ のオイラー標数を 0 に強制する。したがって
$$\mu(P) = (-1)^{\dim M} = (-1)^{\dim X}.$$
$k$ が奇数のときも同じ結果が成り立つ。

**定理 4.2** *構成可能関数 $\mu : X \to \mathbb{Z}$ はスキーム $X$ に内在的である。$X$ が滑らかなら、すべての $P \in X$ に対し値 $\mu(P) = (-1)^{\dim X}$ をとる。*

したがって
$$\mathcal{I}_{\Omega_M}(M, \Gamma_{df}) = \int_{[X]^{\mathrm{vir}}} 1 = (-1)^{\dim X} \chi(X) = \chi(X, \mu).$$

> 〔木原メモ〕**滑らかな場合**：障害束＝接束の双対 → 仮想計数＝$(-1)^{\dim X}\chi(X)$。標準形 $x_1 x_2+\cdots$ のミルナーファイバーに**自由 $S^1$ 作用**があるので $\chi(F_0)=0$、ゆえ $\mu=(-1)^{\dim X}$（符号だけ）。**Behrend 函数 $\mu$ はスキーム $X$ に内在的（埋め込み $M,f$ に依らない）**＝定理 4.2 の核心。あなたの「不変量は内在的で、外部の埋め込みに依存しない」と完全に同じ主張。

**EN.** **Additive nature of $\#^{\mathrm{virt}}(X)$.** The singular Gauß-Bonnet theorem can be stitched together to give rise to an obstruction theory, dual to the deformation theory. For any scheme $Y$ with an obstruction theory, we have
$$\#^{\mathrm{virt}}(Y) = \int_{[Y]^{\mathrm{vir}}} 1 = \chi(Y, \mu_Y).$$

**Theorem 4.3** One of the main implications of this theorem is that the intersection number is *motivic*, i.e.,

(i) the intersection number makes sense for non-compact schemes:
$$\#^{\mathrm{virt}}(Y) = \chi(Y, \mu),$$

(ii) the intersection number is additive over stratifications:
$$\chi(Y, \mu_Y) = \chi(X \setminus Z, \mu_Y) + \chi(Z, \mu_Y),$$
if $Z \subset Y$ is closed.

This is unusual for intersection numbers, it is only true for Lagrangian intersections.

**JA.** **$\#^{\mathrm{virt}}(X)$ の加法性。** 特異ガウス・ボンネ定理は貼り合わせて、変形理論に双対な障害理論を生み出せる。障害理論をもつ任意のスキーム $Y$ に対し、
$$\#^{\mathrm{virt}}(Y) = \int_{[Y]^{\mathrm{vir}}} 1 = \chi(Y, \mu_Y).$$

**定理 4.3** この定理の主要な含意の一つは、交叉数が**モチーフ的（motivic）**であること、すなわち：

(i) 交叉数が**非コンパクト**スキームに対しても意味をもつ：
$$\#^{\mathrm{virt}}(Y) = \chi(Y, \mu),$$

(ii) 交叉数が**層別化（stratification）に関して加法的**：
$$\chi(Y, \mu_Y) = \chi(X \setminus Z, \mu_Y) + \chi(Z, \mu_Y),$$
ただし $Z \subset Y$ は閉。

これは交叉数としては異例であり、**ラグランジアン交叉に対してのみ成り立つ**。

> 〔木原メモ〕**仮想計数のモチーフ性**：(i) 非コンパクトでも定義可、(ii) 層別化に対し**加法的** $\chi(Y)=\chi(Y\setminus Z)+\chi(Z)$。普通の交叉数では成り立たず、**ラグランジアン交叉だけの特権**。あなたの離散基礎の「全体＝部分の素朴な足し算（切断＝はさみ関係）で組める」という加法的構造そのもの。次の §5 がこれを母関数に載せる。

---

## 5　モチーフ的計数（Motivic counts）

**EN.** **Motivic critical loci.** Let $K(\mathrm{Var})$ be the Grothendieck group of $\mathbb{C}$-varieties modulo scissor relations: $[Y] = [Y \setminus Z] + [Z]$, whenever $Z \to Y$ is a closed immersion. There exists a lift, the motivic virtual count of critical loci:
$$
\begin{array}{ccc}
& & K(\mathrm{Var}) \\
& \Phi \nearrow & \downarrow \chi \\
\text{critical loci} & \xrightarrow{\ \#^{\mathrm{virt}}\ } & \mathbb{Z}
\end{array}
$$
Here, $\Phi$ is defined by
$$\Phi(M, f) = -q^{-\dim M / 2} [\phi_f],$$
where $q = [\mathbb{C}]$ is the motivic weight of the affine line, and $[\phi_f]$ are the motivic vanishing cycles of Denef-Loeser (2000), from their work on motivic integration, a motivic version of Milnor fibres.

**JA.** **モチーフ的臨界軌跡。** $K(\mathrm{Var})$ を、はさみ関係 $[Y] = [Y \setminus Z] + [Z]$（$Z \to Y$ が閉埋め込みのとき）で割った $\mathbb{C}$-多様体のグロタンディーク群とする。臨界軌跡のモチーフ的仮想計数というリフトが存在する：
$$
\begin{array}{ccc}
& & K(\mathrm{Var}) \\
& \Phi \nearrow & \downarrow \chi \\
\text{臨界軌跡} & \xrightarrow{\ \#^{\mathrm{virt}}\ } & \mathbb{Z}
\end{array}
$$
ここで $\Phi$ は
$$\Phi(M, f) = -q^{-\dim M / 2} [\phi_f]$$
で定義される。$q = [\mathbb{C}]$ はアフィン直線のモチーフ的重み、$[\phi_f]$ は Denef-Loeser（2000）の**モチーフ的消滅サイクル**（モチーフ的積分の研究より、ミルナーファイバーのモチーフ版）である。

> 〔木原メモ〕整数値の Euler 標数 $\chi$ を、**多様体のグロタンディーク群 $K(\mathrm{Var})$ への持ち上げ $\Phi$**（モチーフ的計数）に格上げ。$q=[\mathbb{C}]$ を変数とし、$q^{1/2}\to-1$ で整数に落ちる。**はさみ関係（加法性）を母関数係数のまま保つ**のが鍵。あなたの「離散的な重み付き分配を、変数 $q$ で母関数化する」発想と同型。

### 定理 5.1（$\mathrm{Hilb}^n(\mathbb{C}^3)$ と 3 次元分割）

**EN.** **Example:** $\mathrm{Hilb}^n(\mathbb{C}^3)$. The Hilbert scheme $\mathrm{Hilb}^n(\mathbb{C}^3)$ is the scheme parametrizing subschemes of length $n$ of affine 3-space. The affine 3-space $\mathbb{C}^3$ is the Calabi-Yau threefold, so that fits with our theme, although $\mathrm{Hilb}^n(\mathbb{C}^3)$ is a critical scheme. Nevertheless, $\mathrm{Hilb}^n(\mathbb{C}^3)$ is an example of a $(-1)$-shifted symplectic scheme. In fact, we can view $\mathrm{Hilb}^n(\mathbb{C}^3)$ as the scheme of three commuting $n \times n$-matrices, namely, respectively, the $x$, $y$ and $z$ coordinates of the subscheme, together with a generating vector. Then $\mathrm{Hilb}^n(\mathbb{C}^3)$ is the critical locus of
$$f : (M_{n \times n}(\mathbb{C})^3 \times \mathbb{C}^n)^{\mathrm{stab}} / GL_n \longrightarrow \mathbb{C}$$
$$(A, B, C, v) \longmapsto \mathrm{tr}([A, B]C).$$

**Theorem 5.1** *For the motivic weight of $\mathrm{Hilb}^n(\mathbb{C}^3)$ we have*
$$\sum_{n=0}^\infty \Phi(\mathrm{Hilb}^n(\mathbb{C}^3)) t^n = \prod_{m=1}^\infty \prod_{k=1}^m \frac{1}{1 - q^{k+1-\frac{m}{2}} t^m}.$$

If we specialize to $q^{\frac{1}{2}} \to -1$, we get
$$\sum_{n=0}^\infty \#^{\mathrm{virt}}(\mathrm{Hilb}^n(\mathbb{C}^3)) t^n = \prod_{m=1}^\infty \left(\frac{1}{1 - (-t)^m}\right)^m.$$
This is (up to signs) the generating function for 3-dimensional partitions.
$$\sum_{n=0}^\infty \#\{\text{3D partitions of } n\} t^n = \prod_{m=1}^\infty \left(\frac{1}{1 - t^m}\right)^m.$$

**JA.** **例：** $\mathrm{Hilb}^n(\mathbb{C}^3)$。ヒルベルトスキーム $\mathrm{Hilb}^n(\mathbb{C}^3)$ は、アフィン 3 空間の長さ $n$ の部分スキームをパラメトライズするスキームである。アフィン 3 空間 $\mathbb{C}^3$ はカラビ・ヤウ 3 次元体なので我々のテーマに合う。ただし $\mathrm{Hilb}^n(\mathbb{C}^3)$ は臨界スキームである。にもかかわらず、$\mathrm{Hilb}^n(\mathbb{C}^3)$ は $(-1)$-シフト・シンプレクティックスキームの例である。実際、$\mathrm{Hilb}^n(\mathbb{C}^3)$ を**3 つの可換な $n \times n$ 行列**（部分スキームの $x, y, z$ 座標）と生成ベクトルの組のスキームとみなせる。すると $\mathrm{Hilb}^n(\mathbb{C}^3)$ は
$$f : (M_{n \times n}(\mathbb{C})^3 \times \mathbb{C}^n)^{\mathrm{stab}} / GL_n \longrightarrow \mathbb{C}$$
$$(A, B, C, v) \longmapsto \mathrm{tr}([A, B]C)$$
の臨界軌跡である。

**定理 5.1** *$\mathrm{Hilb}^n(\mathbb{C}^3)$ のモチーフ的重みについて、*
$$\sum_{n=0}^\infty \Phi(\mathrm{Hilb}^n(\mathbb{C}^3)) t^n = \prod_{m=1}^\infty \prod_{k=1}^m \frac{1}{1 - q^{k+1-\frac{m}{2}} t^m}.$$

$q^{\frac{1}{2}} \to -1$ と特殊化すると、
$$\sum_{n=0}^\infty \#^{\mathrm{virt}}(\mathrm{Hilb}^n(\mathbb{C}^3)) t^n = \prod_{m=1}^\infty \left(\frac{1}{1 - (-t)^m}\right)^m.$$
これは（符号を除いて）**3 次元分割（3D partitions）の母関数**である。
$$\sum_{n=0}^\infty \#\{n \text{ の 3 次元分割}\} t^n = \prod_{m=1}^\infty \left(\frac{1}{1 - t^m}\right)^m.$$

> 〔木原メモ〕**ここが講演とあなたの研究の最接近点**。$\mathrm{Hilb}^n(\mathbb{C}^3)$＝**3 つの可換行列 $[A,B]=0$ ＋生成ベクトル**＝ポテンシャル $\mathrm{tr}([A,B]C)$ の臨界軌跡。その仮想計数の母関数が $\prod_m (1-(-t)^m)^{-m}$＝**3 次元ヤング図形（plane partition）の母関数**（MacMahon）。**$\mathbb{C}^3$ 上の点の数え上げ＝立方体の箱積み＝格子点計数**。あなたの「離散格子・整数配列・箱の積み上げ」と**文字通り同じ組合せ対象**。可換条件 $[A,B]C$ のトレース＝あなたの「位相配列の整合条件」と構造的に響く。

---

### 定理 5.2（一般の CY3・$\mathrm{Hilb}^n Y$）と quintic の指数 $-200$

**EN.** **Motivic curve counts.** Using Donaldson-Thomas theory, we can make sense of motivic curve counts in Calabi-Yau threefolds. Our moduli space will be the Hilbert scheme of 1-dimensional proper subschemes $Z$, with fixed Euler characteristic and holomorphic Euler characteristic. Call this the DT invariant of $Z$, $\#(X)$ is constructed as follows.

First, the moduli space is locally the critical locus of a function, so it will have motivic vanishing cycles from these critical loci, they get combined into global invariants by an orientation of the moduli space. In the extreme case where the curve part of the subscheme is empty, we are back in the case of $\mathrm{Hilb}^n(X)$.

**Theorem 5.2** *The motivic weight of the Hilbert scheme of the Calabi-Yau threefold $X$, compact or not, is given by*
$$\sum_{n=0}^\infty \Phi(\mathrm{Hilb}^n Y) t^n = \left(\prod_{m=1}^\infty \prod_{k=1}^m \frac{1}{1 - q^{k-2-\frac{m}{2}} t^m}\right)^{[Y]}.$$
This theorem uses the power structure on $K(\mathrm{Var})$. If we specialize $q^{\frac{1}{2}} = -1$, we get
$$\sum_{n=0}^\infty \#^{\mathrm{virt}}(\mathrm{Hilb}^n Y) t^n = \left(\prod_{m=1}^\infty \frac{1}{1 - (-t)^m}\right)^{-200}.$$
For example, the quintic gives:

**JA.** **モチーフ的曲線計数。** ドナルドソン・トーマス理論を用いて、カラビ・ヤウ 3 次元体内のモチーフ的曲線計数に意味を与えられる。モジュライ空間は、固定したオイラー標数と正則オイラー標数をもつ 1 次元固有部分スキーム $Z$ のヒルベルトスキームである。これを $Z$ の DT 不変量と呼ぶ。

まず、モジュライ空間は局所的に関数の臨界軌跡なので、これらの臨界軌跡からモチーフ的消滅サイクルを持ち、それらはモジュライ空間の向き付けによって大域的不変量に統合される。部分スキームの曲線部分が空という極端な場合には、$\mathrm{Hilb}^n(X)$ の場合に戻る。

**定理 5.2** *コンパクトとは限らないカラビ・ヤウ 3 次元体 $X$ のヒルベルトスキームのモチーフ的重みは、*
$$\sum_{n=0}^\infty \Phi(\mathrm{Hilb}^n Y) t^n = \left(\prod_{m=1}^\infty \prod_{k=1}^m \frac{1}{1 - q^{k-2-\frac{m}{2}} t^m}\right)^{[Y]}.$$
この定理は $K(\mathrm{Var})$ 上の**べき構造（power structure）**を用いる。$q^{\frac{1}{2}} = -1$ と特殊化すると、
$$\sum_{n=0}^\infty \#^{\mathrm{virt}}(\mathrm{Hilb}^n Y) t^n = \left(\prod_{m=1}^\infty \frac{1}{1 - (-t)^m}\right)^{-200}.$$
例えばクインティックは指数 $-200$ を与える。

> 〔木原メモ〕一般 CY3 では母関数が $\big(\prod(1-(-t)^m)^{-1}\big)^{\chi(Y)}$＝**指数が CY3 のオイラー標数**。**クインティックの $\chi=-200$ がそのまま指数に乗る**（§1 の Hodge ダイヤモンドの $\pm 200$ が再登場）。3D 分割の母関数の「$\chi$ 乗」＝局所（$\mathbb{C}^3$）の計数を大域に貼り合わせる。$[Y]$ 乗のべき構造＝あなたの「局所重みのテンソル積／指数化で大域を組む」発想と同型。

---

## 6　Pardon による曲線計数理論の分類

**EN.** **6 Pardon's classification of curve counting theories.** We have been vague about the actual moduli spaces used to define curve counts in CY3s. Briefly, some of the more important ones are the following.

**Gromov-Witten invariants.** Here we model curves in $X$ by maps $f : C \to X$ from an abstract curve $C$ to $X$. To compactify the moduli space we allow $C$ to acquire nodal singularities. This gives rise to $\overline{M}_{g,n}(X, d)$, the space of stable maps of genus $g$ and degree $d$. Because of the symmetries of the moduli stack the invariants
$$GW_{g,d}(X, d) \in \mathbb{Q}.$$
These are rational numbers, because of the $(-1)$-shifted symplectic structure. The moduli stack may exhibit, $\overline{M}_g(X, d)$ is actually an algebraic stack. The integral or degree of the virtual fundamental class makes sense and we get the Gromov-Witten invariants of the virtual fundamental class. If $X$ is compact, the integral makes sense, but they are not proper, so wrong at the boundary. If $X$ is not compact, we can still make sense of GW invariants if there is a $\mathbb{C}^*$-action on $X$, which induces a $\mathbb{C}^*$-action on $\overline{M}_g(X, d)$ with compact fixed locus. One uses the Atiyah-Bott localization formula to define invariants in this case.

**Donaldson-Thomas invariants.** Here we model curves as subschemes $Z \subset X$, and we consider the moduli space of ideal sheaves of 1-dimensional subschemes $Z \subset X$. This moduli space is $(-1)$-shifted symplectic, and therefore the subscheme $Z_1(X)$. This moduli space make sense whether or not $X$ is compact, because in the non-compact case we can use the weighted Euler characteristic as definition. They are also defined by construction.

**Pandharipande-Thomas PT-invariants.** These are a variation of DT-invariants, which remove the zero-dimensional fuzz which 1-dimensional subschemes tend to acquire under deformations. Technically, these are defined via moduli spaces of certain derived category objects on $X$.

**Gopakumar-Vafa invariants.** These are certain delicate recombinations of Gromov-Witten invariants, which turn out to be integers, although this was not known for a long time. They can also be defined as certain sheaf counting invariants, at least in the genus zero case.

**Pardon's classification.** This is a precise theorem, once it has been defined what a curve counting theory is, and what the equivariant local curves are. Every curve counting theory on Calabi-Yau threefolds is entirely determined by its values on so-called *equivariant local curves* $x_{g,m}$, where $g, m \in \mathbb{Z}_{\ge 0}$. Assigning rational numbers arbitrarily to these $x_{g,m}$ defines a curve counting theory for Calabi-Yau threefolds.

**JA.** **6　Pardon による曲線計数理論の分類。** CY3 における曲線計数を定義するのに使われる実際のモジュライ空間について、これまで曖昧にしてきた。簡潔に、より重要なものをいくつか挙げる。

**グロモフ・ウィッテン（GW）不変量。** ここでは $X$ 内の曲線を、抽象曲線 $C$ から $X$ への写像 $f : C \to X$ でモデル化する。モジュライ空間をコンパクト化するため、$C$ が節点特異点を獲得することを許す。これにより種数 $g$・次数 $d$ の**安定写像**の空間 $\overline{M}_{g,n}(X, d)$ が生じる。モジュライスタックの対称性のため、不変量は
$$GW_{g,d}(X, d) \in \mathbb{Q}$$
である。$(-1)$-シフト・シンプレクティック構造のため、これらは**有理数**である。$\overline{M}_g(X, d)$ は代数的スタックである。仮想基本類の積分（次数）が意味をもち、GW 不変量を得る。$X$ がコンパクトなら積分は意味をもつ。$X$ が非コンパクトでも、$X$ 上に $\mathbb{C}^*$ 作用があって $\overline{M}_g(X, d)$ 上にコンパクトな固定点軌跡をもつ $\mathbb{C}^*$ 作用を誘導するなら、GW 不変量に意味を与えられる。この場合 Atiyah-Bott 局所化公式で不変量を定義する。

**ドナルドソン・トーマス（DT）不変量。** ここでは曲線を部分スキーム $Z \subset X$ としてモデル化し、1 次元部分スキーム $Z \subset X$ のイデアル層のモジュライ空間を考える。このモジュライ空間は $(-1)$-シフト・シンプレクティックである。$X$ がコンパクトかどうかによらず意味をもつ。非コンパクトの場合は重み付きオイラー標数を定義として使えるからである。これらは構成上、整数である。

**パンダリパンデ・トーマス（PT）不変量。** これは DT 不変量の変種で、1 次元部分スキームが変形のもとで獲得しがちな 0 次元の「もや（fuzz）」を除去する。技術的には、$X$ 上のある導来圏の対象のモジュライ空間を介して定義される。

**ゴパクマール・ヴァファ（GV）不変量。** これは GW 不変量のある精妙な組み替えで、長い間知られていなかったが**整数**になることが判明している。少なくとも種数 0 の場合には、ある層計数不変量としても定義できる。

**Pardon の分類。** これは、曲線計数理論とは何か、同変局所曲線とは何かが定義されれば、精密な定理である。**カラビ・ヤウ 3 次元体上のすべての曲線計数理論は、いわゆる同変局所曲線 $x_{g,m}$（$g, m \in \mathbb{Z}_{\ge 0}$）上の値で完全に決定される。** これらの $x_{g,m}$ に有理数を任意に割り当てることが、カラビ・ヤウ 3 次元体に対する曲線計数理論を定義する。

> 〔木原メモ〕**4 つの計数理論**：GW（写像・有理数）／DT（イデアル層・整数）／PT（DT から 0 次元のもやを除く）／GV（GW の整数組み替え）。**DT/GV が整数、GW が有理数**——あなたの「真の基礎は整数、連続値は近似・派生」という離散基礎の規律と同じ温度感。**Pardon の分類＝全理論が「同変局所曲線 $x_{g,m}$」という離散的な基底ベクトル（$g,m\in\mathbb{Z}_{\ge0}$）上の値で決まる**＝**無限自由度の理論が、可算離散の生成元で張られる**。あなたの「許容状態＝離散基底で全体を張る」と構造的に同型。

### 同変局所曲線 $x_{g,m}$ と e-不変量

**EN.** Then we define an *enumerative problem* to be given by a single cohomology class
$$e(X/S) \in H^{\dim S}_*(\mathfrak{Z}(X/S), \mathbb{Q}).$$
Then $e(X/S) \in H^{\dim S}_*(\mathfrak{Z}(X/S), \mathbb{Q})$ in cohomology with compact supports of the relative cycle space of a single CY3 family $X/S$.

Given an enumerative theory $e$ and an enumerative problem $\alpha(X/S)$ we get a number
$$(e(X/S), \alpha(X/S)) \in \mathbb{Q},$$
the e-invariant of $\alpha(X/S)$. For example, $e$ could be GW-theory, and $\alpha(Q/s) = 1_{\mathfrak{Z}_{0,d}(Q)}$, where $Q$ is a quintic.

The number
$$(GW, 1_{\mathfrak{Z}_{0,d}(Q)})$$
is the Gromov-Witten invariant of curves of genus $g$ and degree $d$ on the quintic $Q$.

The equivariant local curve elements $x_{g,m}$ are certain enumerative problems. Start with a smooth and proper curve of genus $g$, with a vector bundle $E$ of rank 2 over it, such that $\det E = \omega_C$. This is a non-compact Calabi-Yau 3-fold over $C$, with $E$ giving the local Calabi-Yau structure. By judicious choice of $E$ we can assume that the cycle space fibres of $E$ over $C$, which are of degree $m$ over $C$ is smooth. The threefold $\mathrm{tot}\,E$ has an action by $\mathbb{C}^*$ on it, simply by rescaling the fibres of $E$ over $C$. We can use this action to make a twisted family parametrized by $\mathbb{P}^N$, with fibre $\mathrm{tot}\,E$. The relative cycle space of this family is parametrized by $\mathbb{P}^N$, associated to the induced $\mathbb{C}^*$-action on $\mathbb{P}^N$ is the twist of $\mathfrak{Z}_m(\mathrm{tot}\,E)$ over $\mathbb{P}^N$ defines the local curve element $\mathfrak{Z}_m(\mathrm{tot}\,E)/\mathbb{P}^N$. Gysin pushforward of $1_{\mathbb{P}^N}$ defines the local curve element
$$x_{g,m} \in H^{2N}_c(\mathfrak{Z}_m(\mathrm{tot}\,E)/\mathbb{P}^N).$$

One can prove that the invariant $(e, x_{g,m})$ is an equivariant version of the enumerative theory $e$ associated to the local curve over $E$. For example,
$$\langle GW(u), x_{g,m}\rangle = \frac{1}{m}\left(2\sin\left(\frac{mu}{2}\right)\right)^{2g-2} \in \mathbb{Q}((u)).$$
But better is the disconnected theory, which is
$$\langle GW(u), x_{g,m}\rangle = \sum_{\sum i \cdot k_i = m} \prod_{i=1}^m \frac{1}{k_i!}\left(\frac{1}{i}\left(2\sin\left(\frac{iu}{2}\right)\right)^{2g-2}\right)^{k_i}.$$
The Gromov-Witten invariant is, in fact, a formal function of a parameter $u$:
$$GW(u)(X/S) = \sum_{h \ge 0} [\overline{M}_h(X/S)]^{\mathrm{vir}} u^{2h-2}.$$
This formal variable is necessary because for a fixed target curve $C$ and covering degree $m$, there are infinitely many possible genera of source curves. Gromov-Witten theory is inherently a sum over genus of the source curve. By contrast,
$$\langle GV, x_{g,m}\rangle = \begin{cases} 1 & \text{if } m = 1 \\ 0 & \text{if } m > 1 \end{cases}.$$
Or, for the disconnected GV invariants:
$$\langle GV, x_{g,m}\rangle = 1.$$
The PT-invariants are again a power series, this time
$$PT(q)(X/S) = \sum_n [\mathcal{I}_n(X/S)]^{\mathrm{vir}} q^n.$$
And
$$\langle PT(q), x_{g,m}\rangle = \sum_{\sum i \cdot k_i = m} \prod_{i=1}^m \frac{1}{k_i! \cdot q^{i k_i}}\left((-1)^{g-1}(-q)^{i(1-g)}(1 - (-q)^i)^{2g-2}\right)^{k_i}.$$
The PT and GW formulas become equal under the substitution $-q = e^{iu}$.

**JA.** そして**数え上げ問題（enumerative problem）**を、単一のコホモロジー類
$$e(X/S) \in H^{\dim S}_*(\mathfrak{Z}(X/S), \mathbb{Q})$$
——単一の CY3 族 $X/S$ の相対サイクル空間のコンパクト台コホモロジーにおける類——で与えられるものとして定義する。

数え上げ理論 $e$ と数え上げ問題 $\alpha(X/S)$ が与えられると、数
$$(e(X/S), \alpha(X/S)) \in \mathbb{Q}$$
——$\alpha(X/S)$ の **e-不変量**——を得る。例えば $e$ を GW 理論、$\alpha(Q/s) = 1_{\mathfrak{Z}_{0,d}(Q)}$（$Q$ はクインティック）とすると、数
$$(GW, 1_{\mathfrak{Z}_{0,d}(Q)})$$
はクインティック $Q$ 上の種数 $g$・次数 $d$ の曲線の GW 不変量である。

**同変局所曲線元 $x_{g,m}$** はある数え上げ問題である。種数 $g$ の滑らかで固有な曲線から始め、その上に $\det E = \omega_C$ を満たす階数 2 のベクトル束 $E$ をとる。これは $C$ 上の非コンパクトなカラビ・ヤウ 3 次元体で、$E$ が局所 CY 構造を与える。$E$ を巧みに選べば、$C$ 上で次数 $m$ の $E$ のサイクル空間ファイバーが滑らかと仮定できる。3 次元体 $\mathrm{tot}\,E$ は、$C$ 上の $E$ のファイバーをスケールするだけで $\mathbb{C}^*$ 作用をもつ。この作用を使って $\mathbb{P}^N$ でパラメトライズされた、ファイバー $\mathrm{tot}\,E$ の捻れ族を作れる。$1_{\mathbb{P}^N}$ の Gysin 押し出しが局所曲線元
$$x_{g,m} \in H^{2N}_c(\mathfrak{Z}_m(\mathrm{tot}\,E)/\mathbb{P}^N)$$
を定義する。

不変量 $(e, x_{g,m})$ が、$E$ 上の局所曲線に付随する数え上げ理論 $e$ の同変版であることを証明できる。例えば、
$$\langle GW(u), x_{g,m}\rangle = \frac{1}{m}\left(2\sin\left(\frac{mu}{2}\right)\right)^{2g-2} \in \mathbb{Q}((u)).$$
だがより良いのは非連結（disconnected）理論で、
$$\langle GW(u), x_{g,m}\rangle = \sum_{\sum i \cdot k_i = m} \prod_{i=1}^m \frac{1}{k_i!}\left(\frac{1}{i}\left(2\sin\left(\frac{iu}{2}\right)\right)^{2g-2}\right)^{k_i}.$$
GW 不変量は実際、パラメータ $u$ の形式的関数である：
$$GW(u)(X/S) = \sum_{h \ge 0} [\overline{M}_h(X/S)]^{\mathrm{vir}} u^{2h-2}.$$
この形式変数が必要なのは、固定した標的曲線 $C$ と被覆次数 $m$ に対し、源曲線の種数が無限通りあり得るからである。GW 理論は本質的に源曲線の種数についての和である。対照的に、
$$\langle GV, x_{g,m}\rangle = \begin{cases} 1 & (m = 1) \\ 0 & (m > 1) \end{cases}.$$
あるいは非連結 GV 不変量では：
$$\langle GV, x_{g,m}\rangle = 1.$$
PT 不変量も再びべき級数で、今度は
$$PT(q)(X/S) = \sum_n [\mathcal{I}_n(X/S)]^{\mathrm{vir}} q^n.$$
そして
$$\langle PT(q), x_{g,m}\rangle = \sum_{\sum i \cdot k_i = m} \prod_{i=1}^m \frac{1}{k_i! \cdot q^{i k_i}}\left((-1)^{g-1}(-q)^{i(1-g)}(1 - (-q)^i)^{2g-2}\right)^{k_i}.$$
**PT と GW の公式は、置換 $-q = e^{iu}$ のもとで一致する。**

> 〔木原メモ〕**MNOP（GW=DT/PT 同値）の核心がこの 1 行**：$-q=e^{iu}$ で PT と GW の公式が一致。**GV が最も単純（$m=1$ で 1、それ以外 0／非連結なら常に 1）**＝整数理論が基底 $x_{g,m}$ 上で最も「素」な値をとる。$2\sin(mu/2)$ の偶数べき＝あなたの好きな三角・周期構造。**$\sum i\,k_i=m$ の分割和**（被覆の分岐の組合せ）＝ここでも plane partition と同じ「分割を数える」骨格。同じ幾何を $u$（連続・GW）と $q$（離散・PT）で見ているのが MNOP。

---

## 7　非可換カラビ・ヤウ 3 次元体（A non-commutative Calabi-Yau threefold）

### 量子 Fermat クインティック $Q$

**EN.** **7 A non-commutative Calabi-Yau threefold.** We start with quantum projective 4-space: the non-commutative graded algebra
$$\mathbb{P}^4_q = \mathbb{C}\langle t_0, \dots, t_4\rangle / (t_i t_j = q^{n_{ij}} t_j t_i),$$
where $q \in \mathbb{C}$, $q = \sqrt[5]{1}$ is a fixed fifth root of unity. Here $N = (n_{ij}) \in M_{5 \times 5}(\mathbb{F}_5)$ is a skew-symmetric matrix. To fix the formulas take
$$N = \begin{pmatrix} 0 & 1 & 1 & 1 & -1 \\ -1 & 0 & 1 & -1 & 1 \\ 1 & -1 & 0 & 1 & -1 \\ -1 & 1 & -1 & 0 & 1 \\ 1 & -1 & 1 & -1 & 0 \end{pmatrix},$$
this is, in fact, generic.

Note that the $t_i^5$ are central elements. Therefore, it makes sense to pass to the quotient to obtain the **Quantum Fermat Quintic**
$$Q = \mathbb{C}\langle t_0, \dots, t_4\rangle_q / (t_0^5 + \dots + t_4^5).$$
This is a graded algebra. We think of
$$Q \rightsquigarrow \mathbb{P}^4_q.$$

**Non-commutative projective schemes.** $Q$ is a non-commutative projective scheme (in the sense of Artin-Zhang). According to Artin-Zhang there is a correspondence
$$(\text{graded } \mathbb{C}\text{-algebras } S) \longleftrightarrow (\text{triples } (\mathscr{C}, \mathscr{O}, (1))).$$
Here, $\mathscr{C}$ is an abelian category, $\mathscr{O} \in \mathscr{C}$ is an object, and $(1)$ is an auto-equivalence. The correspondence sends
$$S \mapsto \mathrm{Proj}\,S = (\mathrm{qgr}(S), S, \mathrm{shift}),$$
where $\mathrm{qgr}(S)$ is the category of tails of finitely generated graded $S$-modules. In the reverse direction,
$$\mathrm{Hom}_\mathscr{C}(\mathscr{O}, \mathscr{O}(n)) \leftarrow (\mathscr{C}, \mathscr{O}, (1)),$$
where multiplication is defined by $a \cdot b = a(1) \circ b$. With enough conditions on triples this gives equivalence of categories (On algebra side up to finite modules).

**JA.** **7　非可換カラビ・ヤウ 3 次元体。** **量子射影 4 空間**から始める：非可換次数付き代数
$$\mathbb{P}^4_q = \mathbb{C}\langle t_0, \dots, t_4\rangle / (t_i t_j = q^{n_{ij}} t_j t_i),$$
ここで $q \in \mathbb{C}$、$q = \sqrt[5]{1}$ は固定した 1 の 5 乗根。$N = (n_{ij}) \in M_{5 \times 5}(\mathbb{F}_5)$ は**歪対称行列**。公式を固定するため上の $N$ をとる（これは実際 generic）。

$t_i^5$ が**中心元**であることに注意する。したがって商に移って**量子 Fermat クインティック**
$$Q = \mathbb{C}\langle t_0, \dots, t_4\rangle_q / (t_0^5 + \dots + t_4^5)$$
を得るのが意味をもつ。これは次数付き代数である。$Q \rightsquigarrow \mathbb{P}^4_q$ と考える。

**非可換射影スキーム。** $Q$ は（Artin-Zhang の意味で）非可換射影スキームである。Artin-Zhang によれば対応
$$(\text{次数付き } \mathbb{C}\text{-代数 } S) \longleftrightarrow (\text{三つ組 } (\mathscr{C}, \mathscr{O}, (1)))$$
がある。ここで $\mathscr{C}$ はアーベル圏、$\mathscr{O} \in \mathscr{C}$ は対象、$(1)$ は自己同値。対応は
$$S \mapsto \mathrm{Proj}\,S = (\mathrm{qgr}(S), S, \mathrm{shift})$$
を送る。$\mathrm{qgr}(S)$ は有限生成次数付き $S$-加群の tails の圏。逆方向は
$$\mathrm{Hom}_\mathscr{C}(\mathscr{O}, \mathscr{O}(n)) \leftarrow (\mathscr{C}, \mathscr{O}, (1)),$$
乗法は $a \cdot b = a(1) \circ b$ で定義される。三つ組に十分な条件を課せばこれは圏同値を与える（代数側は有限加群を除いて）。

> 〔木原メモ〕**幾何を捨てて代数だけで CY3 を作る**：可換多項式環の代わりに**非可換歪可換環**（$t_i t_j = q^{n_{ij}} t_j t_i$、$q=$ 1 の 5 乗根、$N$ は $\mathbb{F}_5$ 上の歪対称行列）。**Proj＝点集合ではなく圏 $(\mathscr{C},\mathscr{O},(1))$**（Artin-Zhang）。あなたの「空間＝点でなく関係（位相配列・整合条件）で定義する」「$\mathbb{F}_5$／5 乗根の離散構造」と直接呼応。歪対称行列 $N$＝シンプレクティック形式の離散版。

### 定理 7.1（Kanazawa）と Frobenius 代数の層

**EN.** **Theorem 7.1 (Kanazawa)** *For the quantum Fermat quintic $Q$ (any $N = (n_{ij})$), the category $\mathrm{qgr}(Q)$ satisfies*
*(i) it has global dimension 3,*
*(ii) it is a Calabi-Yau 3 category if and only if $\binom{1}{1}$ is an eigenvector of $N$.*

Here the two conditions essentially mean that
(i) $\mathrm{Ext}^i(E, F) = 0$ $\forall i > 3$,
(ii) $\mathrm{Ext}^i(E, F)^\vee = \mathrm{Ext}^{3-i}(F, E)$.

Or, in words,
(i) $Q$ is smooth of dimension 3,
(ii) $Q$ is a Calabi-Yau threefold.

Because of this, moduli spaces of objects in $\mathrm{qgr}(Q)$ should admit Donaldson-Thomas schemes with $(-1)$-shifted symplectic structures. We were not able to construct it using techniques from non-commutative projective geometry.

**Sheaves of Frobenius algebras.** An important observation is that $Q$ has a central (in particular commutative) subalgebra over which it is finite:
$$\mathbb{C}[t_0^5, \dots, t_4^5]/(t_0^5 + \dots + t_4^5) \hookrightarrow \mathbb{C}\langle t_0, \dots, t_4\rangle_q/(t_0^5 + \dots + t_4^5).$$
We have
$$\mathbb{C}[t_0^5, \dots, t_4^5]/(t_0^5 + \dots + t_4^5) = \mathbb{C}[x_0, \dots, x_4]/(x_0 + \dots + x_4) = X \to \mathbb{P}^4.$$
This is a hyperplane $\mathbb{P}^3 \cong X \to \mathbb{P}^4$. Via (3), the quantum Fermat quintic $Q$, locally free turns into a non-commutative sheaf $\mathscr{A}$ of $\mathscr{O}_X$-algebras over $X$, a graded free module of rank 625. In fact, the 5-Veronese subalgebra of $\mathbb{C}\langle t_0, \dots, t_4\rangle_q$ is a graded free module over $\mathbb{C}[t_0^5, \dots, t_4^5]$ on the basis $t^{\vec{k}}$, where $\sum k_i = 5$, $0 \le k_i \le 4$. Therefore,
$$\mathscr{A} \cong \mathscr{O}_X \oplus \mathscr{O}_X(-1)^{\oplus 121} \oplus \mathscr{O}_X(-2)^{\oplus 381} \oplus \mathscr{O}_X(-3)^{\oplus 121} \oplus \mathscr{O}_X(-4)$$
as $\mathscr{O}_X$-module (not as algebra). Multiplication in $\mathscr{A}$ composed with projection $\mathrm{tr} : \mathscr{A} \to \mathscr{O}_X(-4)$ defines a perfect pairing
$$\mathscr{A} \otimes_{\mathscr{O}_X} \mathscr{A} \to \mathscr{O}_X(-4) = \omega_X,$$
given by $a \otimes b \mapsto \mathrm{tr}(ab)$. This pairing is symmetric if and only if $t^{\vec{k}} t^{4-\vec{k}} = t^{4-\vec{k}} t^{\vec{k}}$, which happens if and only if $\binom{1}{1}$ is an eigenvector of $N$.

**JA.** **定理 7.1（Kanazawa）** *量子 Fermat クインティック $Q$（任意の $N=(n_{ij})$）に対し、圏 $\mathrm{qgr}(Q)$ は次を満たす：*
*(i) 大域次元 3 をもつ、*
*(ii) $\binom{1}{1}$ が $N$ の固有ベクトルであるとき、かつそのときに限りカラビ・ヤウ 3 圏である。*

ここで二つの条件は本質的に、
(i) $\mathrm{Ext}^i(E, F) = 0$ ($\forall i > 3$)、
(ii) $\mathrm{Ext}^i(E, F)^\vee = \mathrm{Ext}^{3-i}(F, E)$
を意味する。言葉で言えば、
(i) $Q$ は次元 3 で滑らか、
(ii) $Q$ はカラビ・ヤウ 3 次元体である。

このため、$\mathrm{qgr}(Q)$ の対象のモジュライ空間は $(-1)$-シフト・シンプレクティック構造をもつ DT スキームを許すはずである。我々は非可換射影幾何の技法ではそれを構成できなかった。

**Frobenius 代数の層。** 重要な観察は、$Q$ がそれ上で有限であるような中心（特に可換）部分代数をもつことである：
$$\mathbb{C}[t_0^5, \dots, t_4^5]/(t_0^5 + \dots + t_4^5) \hookrightarrow \mathbb{C}\langle t_0, \dots, t_4\rangle_q/(t_0^5 + \dots + t_4^5).$$
すると
$$\mathbb{C}[t_0^5, \dots, t_4^5]/(t_0^5 + \dots + t_4^5) = \mathbb{C}[x_0, \dots, x_4]/(x_0 + \dots + x_4) = X \to \mathbb{P}^4.$$
これは超平面 $\mathbb{P}^3 \cong X \to \mathbb{P}^4$ である。(3) を介して、量子 Fermat クインティック $Q$ は $X$ 上の $\mathscr{O}_X$-代数の非可換層 $\mathscr{A}$、すなわち階数 625 の次数付き自由加群になる。実際、$\mathbb{C}\langle t_0, \dots, t_4\rangle_q$ の 5-Veronese 部分代数は、$\sum k_i = 5$, $0 \le k_i \le 4$ なる基底 $t^{\vec{k}}$ の上で $\mathbb{C}[t_0^5, \dots, t_4^5]$ 上の次数付き自由加群である。したがって $\mathscr{O}_X$-加群として（代数としてではなく）
$$\mathscr{A} \cong \mathscr{O}_X \oplus \mathscr{O}_X(-1)^{\oplus 121} \oplus \mathscr{O}_X(-2)^{\oplus 381} \oplus \mathscr{O}_X(-3)^{\oplus 121} \oplus \mathscr{O}_X(-4).$$
$\mathscr{A}$ における乗法を射影 $\mathrm{tr} : \mathscr{A} \to \mathscr{O}_X(-4)$ と合成すると、完全対 pairing
$$\mathscr{A} \otimes_{\mathscr{O}_X} \mathscr{A} \to \mathscr{O}_X(-4) = \omega_X,\qquad a \otimes b \mapsto \mathrm{tr}(ab)$$
を定義する。この対は、$t^{\vec{k}} t^{4-\vec{k}} = t^{4-\vec{k}} t^{\vec{k}}$ のとき、かつそのときに限り対称で、それは $\binom{1}{1}$ が $N$ の固有ベクトルのとき、かつそのときに限り起こる。

> 〔木原メモ〕**CY 条件が「$\binom{1}{1}$ が歪対称行列 $N$ の固有ベクトル」という離散的・代数的条件に翻訳される**（Kanazawa）。$\mathscr{A}$ の階数分解 $1+121+381+121+1=625=5^4$＝**Veronese 基底の組合せ計数**（$\sum k_i=5$, $0\le k_i\le4$）。トレース対 $\mathrm{tr}(ab):\mathscr{A}\otimes\mathscr{A}\to\omega_X$ の対称性＝Frobenius／CY の本質。あなたの「対称双線形形式・トレース・$\mathbb{F}_5$ 上の組合せ」と全面的に同型。

### 定義 7.2（Frobenius 代数の層）・定理 7.3（Liu）・量子クインティックのミラー

**EN.** **Definition 7.2** Let $X$ be a smooth scheme. A locally free sheaf of $\mathscr{O}_X$-algebras $\mathscr{A}$, with a symmetric perfect pairing $\mathscr{A} \otimes \mathscr{A} \to \omega_X$ is a sheaf of *Frobenius algebras* over $X$.

If a sheaf of algebras $\mathscr{A}$ over $\mathscr{O}_X$ has finite global dimension $n = \dim X$, it has a dualizing bimodule $\omega_\mathscr{A} = \mathscr{H}om_{\mathscr{O}_X}(\mathscr{A}, \omega_X)$ such that
$$\mathrm{Ext}^i_\mathscr{A}(\mathscr{F}, \mathscr{G}) = \mathrm{Ext}^{n-i}_\mathscr{A}(\mathscr{G}, \omega_\mathscr{A} \otimes_\mathscr{A} \mathscr{F})^\vee,$$
for all $\mathscr{F}, \mathscr{G} \in \mathrm{Coh}(\mathscr{A})$, where $\mathrm{Coh}(\mathscr{A})$ is the category of left $\mathscr{A}$-modules which are coherent $\mathscr{O}_X$-modules. A symmetric pairing identifies $\omega_\mathscr{A} = \mathscr{H}om_{\mathscr{O}_X}(\mathscr{A}, \omega_X) = \mathscr{A}$ as $\mathscr{A}$-bimodule. So $\mathrm{Coh}(\mathscr{A})$ becomes a Calabi-Yau $n$-category. We remark that in our situation $\mathrm{qgr}(Q)$ and $\mathrm{Coh}(\mathscr{A})$ are equivalent. So we may study $\mathrm{Coh}(\mathscr{A})$ instead.

**The quintic mirror.** This quantum Fermat quintic is, in fact, a non-commutative resolution of the singular scheme giving rise to the quintic mirror (for which, we recall, the projective scheme associated to the graded ring
$$W = \mathbb{C}[x_0, \dots, x_4, y]/(x_0^5 + \dots + x_4^5 - x_0 \dots x_4)).$$

The quantum Fermat quintic may thus be considered as a non-commutative version of the quintic mirror.

**Theorem 7.3 (Liu)** *$M^{s,h}(X, \mathscr{A})$ carries a virtual fundamental class $[M^{s,h}(X, \mathscr{A})]^{\mathrm{vir}} \in A_0(M^{s,h}(X, \mathscr{A}))$.*

**Definition 7.4** Suppose $h$ is chosen such that $M^{s,h}(X, \mathscr{A})$ is proper, i.e., assume that semi-stable implies stable. For example, if $\mathscr{A} \otimes \mathbb{C}(X)$ is a division algebra and we consider sheaves of dimension $\dim X$ and rank $\mathrm{rk}\,\mathscr{A}$. Then we define
$$DT(M^{s,h}(X, \mathscr{A})) := \int_{[M^{s,h}(X, \mathscr{A})]^{\mathrm{vir}}} 1 \in \mathbb{Z}.$$
If (i) and (ii) are satisfied, also
$$DT(\mathrm{Hilb}^h(X, \mathscr{A})) = \int_{[\mathrm{Hilb}^h(X, \mathscr{A})]^{\mathrm{vir}}} 1.$$

**JA.** **定義 7.2** $X$ を滑らかなスキームとする。$\mathscr{O}_X$-代数の局所自由層 $\mathscr{A}$ で、対称完全対 $\mathscr{A} \otimes \mathscr{A} \to \omega_X$ をもつものを、$X$ 上の **Frobenius 代数の層**という。

$\mathscr{O}_X$ 上の代数の層 $\mathscr{A}$ が有限大域次元 $n = \dim X$ をもつなら、双対化双加群 $\omega_\mathscr{A} = \mathscr{H}om_{\mathscr{O}_X}(\mathscr{A}, \omega_X)$ をもち、すべての $\mathscr{F}, \mathscr{G} \in \mathrm{Coh}(\mathscr{A})$ に対し
$$\mathrm{Ext}^i_\mathscr{A}(\mathscr{F}, \mathscr{G}) = \mathrm{Ext}^{n-i}_\mathscr{A}(\mathscr{G}, \omega_\mathscr{A} \otimes_\mathscr{A} \mathscr{F})^\vee$$
を満たす。ここで $\mathrm{Coh}(\mathscr{A})$ は、連接 $\mathscr{O}_X$-加群である左 $\mathscr{A}$-加群の圏。対称対は $\omega_\mathscr{A} = \mathscr{H}om_{\mathscr{O}_X}(\mathscr{A}, \omega_X) = \mathscr{A}$ を $\mathscr{A}$-双加群として同一視する。したがって $\mathrm{Coh}(\mathscr{A})$ はカラビ・ヤウ $n$-圏になる。我々の状況では $\mathrm{qgr}(Q)$ と $\mathrm{Coh}(\mathscr{A})$ が同値なので、代わりに $\mathrm{Coh}(\mathscr{A})$ を調べてよい。

**クインティックのミラー。** この量子 Fermat クインティックは、実際、クインティックのミラーを生じる特異スキームの**非可換解消（non-commutative resolution）**である（その特異スキームとは、次数付き環
$$W = \mathbb{C}[x_0, \dots, x_4, y]/(x_0^5 + \dots + x_4^5 - x_0 \dots x_4)$$
に付随する射影スキームを思い出されたい）。量子 Fermat クインティックはこうしてクインティックのミラーの非可換版とみなせる。

**定理 7.3（Liu）** *$M^{s,h}(X, \mathscr{A})$ は仮想基本類 $[M^{s,h}(X, \mathscr{A})]^{\mathrm{vir}} \in A_0(M^{s,h}(X, \mathscr{A}))$ をもつ。*

**定義 7.4** $M^{s,h}(X, \mathscr{A})$ が固有（proper）になるように $h$ を選ぶ、すなわち半安定が安定を含意すると仮定する。例えば $\mathscr{A} \otimes \mathbb{C}(X)$ が斜体（division algebra）で、次元 $\dim X$・階数 $\mathrm{rk}\,\mathscr{A}$ の層を考える場合。このとき
$$DT(M^{s,h}(X, \mathscr{A})) := \int_{[M^{s,h}(X, \mathscr{A})]^{\mathrm{vir}}} 1 \in \mathbb{Z}.$$
(i), (ii) が満たされれば、
$$DT(\mathrm{Hilb}^h(X, \mathscr{A})) = \int_{[\mathrm{Hilb}^h(X, \mathscr{A})]^{\mathrm{vir}}} 1.$$

> 〔木原メモ〕**非可換 $Q$ ＝ クインティック・ミラー（$\sum x_i^5 - x_0\cdots x_4$）の特異点の非可換解消**。$\mathrm{qgr}(Q)\simeq\mathrm{Coh}(\mathscr{A})$ で**幾何的な層の圏に翻訳**して DT 不変量（整数）を定義（Liu）。あなたの「特異・退化した対象を、離散代数データで解消して整数不変量を取り出す」と同じ運び。Frobenius／対称対＝あなたの内積構造そのもの。

### $Z_Y(t)$ と $Z_Q(t)$ の計算・系 7.6（最終答）

**EN.** **Computation of $Z_Y(t)$ for the commutative quintic.** Let us review how to compute the weighted Euler characteristic of $\mathrm{Hilb}^n(Y)$, where $Y$ is a commutative quintic. Define
$$\mathrm{Hilb}^n(Y|P) \subset \mathrm{Hilb}^n Y,$$
the punctual Hilbert Scheme, to consist of all subschemes of $Y$ of length $n$, supported at $P \in Y$. Then
$$Z_Y(t) = \sum_n \chi(\mathrm{Hilb}^n Y, \nu_{\mathrm{Hilb}^n Y}) t^n$$
and
$$Z_{Y|P}(t) = \sum_n \chi(\mathrm{Hilb}^n(Y|P), \nu_{\mathrm{Hilb}^n Y}) t^n.$$
By a stratification argument, one proves that
$$Z_Y(t) = Z_{Y|P}(t)^{\chi(Y)}.$$
This reduces to the study of the punctual Hilbert scheme. Now we use that
$$\mathrm{Germ}(\mathrm{Hilb}^n(Y|P), \mathrm{Hilb}^n Y) = \mathrm{Germ}(\mathrm{Hilb}^n(\mathbb{C}^3|0), \mathrm{Hilb}^n \mathbb{C}^3).$$
Hence also,
$$\chi(\mathrm{Hilb}^n(Y|P), \nu_{\mathrm{Hilb}^n Y}) = \chi(\mathrm{Hilb}^n(\mathbb{C}^3|0), \nu_{\mathrm{Hilb}^n \mathbb{C}^3}).$$
This reduces to the punctual Hilbert scheme of $\mathbb{C}^3$. We can use the $\mathbb{C}^*$-action, and Property (iii) of $\nu$. We get
$$\chi(\mathrm{Hilb}^n(\mathbb{C}^3|0), \nu_{\mathrm{Hilb}^n \mathbb{C}^3}) = (-1)^n \cdot \#\{\text{3D partitions of } n\}.$$
Thus,
$$Z_{Y|P}(t) = Z_{\mathbb{C}^3|0}(t) = M(-t),$$
where $M(t) = \prod_{m=1}^\infty \frac{1}{(1 - t^m)^m}$ is the MacMahon function. Finally, we deduce
$$Z_Y(t) = M(-t)^{-200}.$$

**JA.** **可換クインティックに対する $Z_Y(t)$ の計算。** $Y$ を可換クインティックとして、$\mathrm{Hilb}^n(Y)$ の重み付きオイラー標数の計算法を復習する。**点状（punctual）ヒルベルトスキーム**
$$\mathrm{Hilb}^n(Y|P) \subset \mathrm{Hilb}^n Y$$
を、$P \in Y$ に台をもつ長さ $n$ の $Y$ の部分スキーム全体として定義する。すると
$$Z_Y(t) = \sum_n \chi(\mathrm{Hilb}^n Y, \nu_{\mathrm{Hilb}^n Y}) t^n,\qquad Z_{Y|P}(t) = \sum_n \chi(\mathrm{Hilb}^n(Y|P), \nu_{\mathrm{Hilb}^n Y}) t^n.$$
層別化の議論により、
$$Z_Y(t) = Z_{Y|P}(t)^{\chi(Y)}$$
が証明される。これは点状ヒルベルトスキームの研究に帰着する。さらに
$$\mathrm{Germ}(\mathrm{Hilb}^n(Y|P), \mathrm{Hilb}^n Y) = \mathrm{Germ}(\mathrm{Hilb}^n(\mathbb{C}^3|0), \mathrm{Hilb}^n \mathbb{C}^3)$$
を使うと、
$$\chi(\mathrm{Hilb}^n(Y|P), \nu_{\mathrm{Hilb}^n Y}) = \chi(\mathrm{Hilb}^n(\mathbb{C}^3|0), \nu_{\mathrm{Hilb}^n \mathbb{C}^3}).$$
これは $\mathbb{C}^3$ の点状ヒルベルトスキームに帰着する。$\mathbb{C}^*$ 作用と $\nu$ の性質 (iii) を使えば、
$$\chi(\mathrm{Hilb}^n(\mathbb{C}^3|0), \nu_{\mathrm{Hilb}^n \mathbb{C}^3}) = (-1)^n \cdot \#\{n \text{ の 3 次元分割}\}.$$
したがって
$$Z_{Y|P}(t) = Z_{\mathbb{C}^3|0}(t) = M(-t),$$
ここで $M(t) = \prod_{m=1}^\infty \frac{1}{(1 - t^m)^m}$ は**MacMahon 関数**。最終的に
$$Z_Y(t) = M(-t)^{-200}.$$

> 〔木原メモ〕**可換クインティックの DT 母関数 $Z_Y(t)=M(-t)^{-200}$**：$M(t)=\prod(1-t^m)^{-m}$＝MacMahon 関数＝3D 分割の母関数、指数 $-200=\chi(Y)$。**大域は局所（点状＝$\mathbb{C}^3$ の原点での 3D 分割）の $\chi(Y)$ 乗**——§5 の構図の完成形。点状ヒルベルトスキームの芽が $\mathbb{C}^3$ のそれと一致＝**局所普遍性**。あなたの「局所の離散計数（箱積み）を $\chi$ 乗で大域に持ち上げる」と完全一致。

**EN.** **Computation of $Z_Q(t)$.** We will now indicate how to compute
$$Z_Q(t) = \sum_n \chi(\mathrm{Hilb}^n(X, \mathscr{A}), \nu) t^n.$$
We start by noting that finite length $\mathscr{A}$-modules have 0-dimensional support in $X$. Therefore we can study the situation locally over $X = \{x_0 + \dots + x_4 = 0\} \subset \mathbb{P}^4$. Let us localize setting $x_0 = 1$. Then $x_1 + \dots + x_4 = -1$. Let us introduce new inhomogeneous variables
$$u_i = \frac{t_0^4 t_i}{t_0^5} = \frac{t_0^4 t_i}{x_0}.$$
We are now in the affine situation with
$$X_0 = \mathbb{C}[x_1, \dots, x_4]/(x_1 + \dots + x_4 = -1),$$
and
$$A_0 = \mathbb{C}[u_1, \dots, u_4]/(u_1^5 + \dots + u_4^5 = -1),\qquad u_i u_j = q^{\bar{n}_{ij}} u_j u_i,$$
where $x_i = u_i^5$. The new commuting relations are given by $\bar{n}_{ij} = n_{ij} - n_{i0} - n_{0j}$. They form a matrix $\bar{N} \in M_{4 \times 4}(\mathbb{F}_5)$, skew-symmetric, $\bar{N}\mathbf{1} = 0$.
$$\bar{N} = \begin{pmatrix} 0 & -2 & -1 & -2 \\ 2 & 0 & -1 & -1 \\ 1 & 1 & 0 & -2 \\ 2 & 1 & 2 & 0 \end{pmatrix}.$$

Point modules are representations of $A$ on $\mathbb{C}$. Under such a representation, $u_1, \dots, u_4$ turn into numbers (which commute). So the non-trivial commutation relations among the $u_i$ force that at most one $u_i$ is non-zero. Say $u_2 = u_3 = u_4 = 0$ and $u_1^5 = -1$, so $u_1 = -q^i$, where $i \in \mathbb{F}_5$. This discussion shows that there are 5 point modules $S_0, \dots, S_4$ supported at $(1 : -1 : 0 : 0 : 0) \in X$. There are $\binom{5}{2} = 10$ such points in $X \subset \mathbb{P}^4$. This proves that for $Q = (X, \mathscr{A})$, we get a total of 50 point modules for $Q = (X, \mathscr{A})$. This proves that
$$DT(\mathrm{Hilb}^1(X, \mathscr{A})) = 50.$$

**JA.** **$Z_Q(t)$ の計算。** $Z_Q(t) = \sum_n \chi(\mathrm{Hilb}^n(X, \mathscr{A}), \nu) t^n$ の計算法を示そう。まず有限長 $\mathscr{A}$-加群は $X$ 内で 0 次元台をもつことに注意する。よって $X = \{x_0 + \dots + x_4 = 0\} \subset \mathbb{P}^4$ 上で局所的に調べられる。$x_0 = 1$ と局所化すると $x_1 + \dots + x_4 = -1$。新しい非斉次変数 $u_i = t_0^4 t_i / t_0^5 = t_0^4 t_i / x_0$ を導入する。いまアフィンな状況にあり、
$$X_0 = \mathbb{C}[x_1, \dots, x_4]/(x_1 + \dots + x_4 = -1),$$
$$A_0 = \mathbb{C}[u_1, \dots, u_4]/(u_1^5 + \dots + u_4^5 = -1),\qquad u_i u_j = q^{\bar{n}_{ij}} u_j u_i,$$
ここで $x_i = u_i^5$。新しい交換関係は $\bar{n}_{ij} = n_{ij} - n_{i0} - n_{0j}$ で与えられ、上の歪対称行列 $\bar{N} \in M_{4 \times 4}(\mathbb{F}_5)$（$\bar{N}\mathbf{1} = 0$）をなす。

**点加群（point modules）**は $A$ の $\mathbb{C}$ 上の表現である。そのような表現のもとで $u_1, \dots, u_4$ は数になる（可換）。よって $u_i$ 間の非自明な交換関係は、**高々一つの $u_i$ だけが非零**であることを強制する。$u_2 = u_3 = u_4 = 0$、$u_1^5 = -1$ とすると $u_1 = -q^i$（$i \in \mathbb{F}_5$）。この議論から、$(1 : -1 : 0 : 0 : 0) \in X$ に台をもつ **5 個の点加群 $S_0, \dots, S_4$** があることがわかる。$X \subset \mathbb{P}^4$ 内にはそのような点が $\binom{5}{2} = 10$ 個ある。よって $Q = (X, \mathscr{A})$ に対し合計 **50 個の点加群**を得る。これは
$$DT(\mathrm{Hilb}^1(X, \mathscr{A})) = 50$$
を証明する。

> 〔木原メモ〕**$50$ が再登場**——§2 の Fermat の「50 鉛筆束」と同じ数が、非可換側では「**50 個の点加群**」として出る（$\binom{5}{2}=10$ 点 × 5 個 $S_0,\dots,S_4$）。非可換性（$u_iu_j=q^{\bar n_{ij}}u_ju_i$）が「**高々一つの $u_i$ だけ非零**」を強制＝離散的選択則。$\mathbb{F}_5$ 上の歪対称行列 $\bar N$、$\bar N\mathbf 1=0$。あなたの「交換関係が許容状態を離散選択に絞る」と完全同型。

**EN.** Consider $\mathscr{A}$ near $P = (1 : -1 : 0 : 0 : 0)$. From general principles we would expect that (assuming all simple $\mathscr{A}$-modules at $P$ are point modules), that
$$\mathrm{Germ}(\mathrm{Hilb}^n(\mathscr{A}|P), \mathrm{Hilb}^n(\mathscr{A})) \cong \mathrm{Germ}\left(\prod_{|\vec{d}|=n} M^s(Q, \vec{d}, v)|0, \prod_{|\vec{d}|=n} M^s(Q, \vec{d}, v)\right),$$
where $(Q, f)$ is the Ext-quiver of $S = S_0 \oplus \dots \oplus S_4$, with potential $f$, $\vec{d}$ is a dimension vector, and $v$ a framing. This expectation relies on a theorem of Toda, which says that on a commutative Calabi-Yau 3-fold $Y$,
$$\mathrm{Germ}(\mathcal{M}^{ss}_\omega|P, \mathcal{M}^{ss}_\omega) \cong \mathrm{Germ}(\mathcal{M}_Q|0, \mathcal{M}_Q).$$
Here, $\mathcal{M}_\omega$ is a stack of Gieseker semi-stable sheaves on $Y$, $\mathcal{M}_\omega|P$ fixes the associated polystable sheaf $\bigoplus_i \mathscr{F}_i^{\oplus k_i}$, $\mathcal{M}_Q$ is the stack of representations of the Ext-quiver of $\bigoplus \mathscr{F}_i$, with potential, and dimension vector $\vec{k}$; and $\mathcal{M}_Q|0$ are the nilpotent representations.

**Theorem 7.5 (Liu)** *The expectation holds. The quiver the following.* The vertices correspond to the point modules $S_0, \dots, S_4$, the arrows to the basic extensions between the $S_i$. For example, $a_i \in \mathrm{Ext}^1(S_i, S_{i-2})$.

So $a_i$ is a 2-dimensional representation of $A_0$,
$$S_{i-2} \xrightarrow{\ a_i\ } S_i,$$
given by
$$u_1 = \begin{pmatrix} -q^{i-2} & 0 \\ 0 & -q^i \end{pmatrix},\quad u_2 = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix},\quad u_2 u_1 = \begin{pmatrix} 0 & -q^i \\ 0 & 0 \end{pmatrix},\quad u_3 = u_4 = 0.$$
We calculate
$$u_1 u_2 = q^{n_{12}} u_2 u_1 = q^{-2} u_2 u_1$$
is satisfied. Similarly, we have $b_i \in \mathrm{Ext}^1(S_i, S_{i-1})$ and $c_i \in \mathrm{Ext}^1(S_i, S_{i-2})$.
The potential is
$$f = \left(\sum q^{i-1} b_i\right)\left(\sum a_i\right)\left(\sum c_i\right) - q^{-1}\left(\sum q^{i-1} b_i\right)\left(\sum c_i\right)\left(\sum a_i\right).$$
Analytically locally near $P = (1, -1, 0, 0, 0)$ we have that $\mathscr{A} \cong J(Q, f)$, the Jacobi algebra of the quiver with potential. Under this isomorphism, commutation relations among $u_2, u_3, u_4$ correspond to $\sum a_i, \sum b_i, \sum c_i$. Commutation relations among $u_2, u_3, u_4$ give relations among $\sum a_i, \sum b_i, \sum c_i$. In fact, there are 15 relations, for example $\partial_{a_i}$ gives
$$q^{i+2} c_{i+2} b_{i+3} = q^{-1} q^i b_{i+1} c_{i+3},$$
so for $i = 0$ this gives $q^2 c_2 b_3 = b_1 c_3$, which we visualize as 〔図：5 頂点 $S_0, \dots, S_4$ のクイバー〕. The framing vector is $\vec{1} = (1, 1, 1, 1, 1)$.

**JA.** $P = (1 : -1 : 0 : 0 : 0)$ の近くの $\mathscr{A}$ を考える。一般原理から（$P$ におけるすべての単純 $\mathscr{A}$-加群が点加群と仮定して）、
$$\mathrm{Germ}(\mathrm{Hilb}^n(\mathscr{A}|P), \mathrm{Hilb}^n(\mathscr{A})) \cong \mathrm{Germ}\left(\prod_{|\vec{d}|=n} M^s(Q, \vec{d}, v)|0, \prod_{|\vec{d}|=n} M^s(Q, \vec{d}, v)\right)$$
が期待される。ここで $(Q, f)$ は $S = S_0 \oplus \dots \oplus S_4$ の **Ext-クイバー**（ポテンシャル $f$ 付き）、$\vec{d}$ は次元ベクトル、$v$ は枠付け（framing）。この期待は **Toda の定理**に依拠する：可換カラビ・ヤウ 3 次元体 $Y$ 上で
$$\mathrm{Germ}(\mathcal{M}^{ss}_\omega|P, \mathcal{M}^{ss}_\omega) \cong \mathrm{Germ}(\mathcal{M}_Q|0, \mathcal{M}_Q).$$
ここで $\mathcal{M}_\omega$ は $Y$ 上の Gieseker 半安定層のスタック、$\mathcal{M}_\omega|P$ は付随する polystable 層 $\bigoplus_i \mathscr{F}_i^{\oplus k_i}$ を固定し、$\mathcal{M}_Q$ は $\bigoplus \mathscr{F}_i$ の Ext-クイバー（ポテンシャル・次元ベクトル $\vec{k}$ 付き）の表現のスタック、$\mathcal{M}_Q|0$ は冪零表現。

**定理 7.5（Liu）** *期待は成り立つ。クイバーは次のとおり。* 頂点は点加群 $S_0, \dots, S_4$ に、矢印は $S_i$ 間の基本拡大に対応する。例えば $a_i \in \mathrm{Ext}^1(S_i, S_{i-2})$。

すなわち $a_i$ は $A_0$ の 2 次元表現 $S_{i-2} \xrightarrow{a_i} S_i$ で、上の $u_1, u_2$（と $u_2 u_1$、$u_3 = u_4 = 0$）で与えられる。計算すると $u_1 u_2 = q^{n_{12}} u_2 u_1 = q^{-2} u_2 u_1$ が満たされる。同様に $b_i \in \mathrm{Ext}^1(S_i, S_{i-1})$ と $c_i \in \mathrm{Ext}^1(S_i, S_{i-2})$ がある。ポテンシャルは
$$f = \left(\sum q^{i-1} b_i\right)\left(\sum a_i\right)\left(\sum c_i\right) - q^{-1}\left(\sum q^{i-1} b_i\right)\left(\sum c_i\right)\left(\sum a_i\right).$$
$P = (1, -1, 0, 0, 0)$ の近くで解析的局所的に $\mathscr{A} \cong J(Q, f)$（ポテンシャル付きクイバーのヤコビ代数）。この同型のもとで $u_2, u_3, u_4$ 間の交換関係が $\sum a_i, \sum b_i, \sum c_i$ 間の関係に対応する。実際 15 個の関係があり、例えば $\partial_{a_i}$ は $q^{i+2} c_{i+2} b_{i+3} = q^{-1} q^i b_{i+1} c_{i+3}$ を与え、$i = 0$ では $q^2 c_2 b_3 = b_1 c_3$ となる（5 頂点のクイバーとして図示される）。枠付けベクトルは $\vec{1} = (1, 1, 1, 1, 1)$。

> 〔木原メモ〕**非可換 CY3 の局所構造＝5 頂点のクイバー（$S_0,\dots,S_4$）＋ポテンシャル $f$ のヤコビ代数 $J(Q,f)$**（Toda／Liu）。頂点＝点加群、矢印 $a_i,b_i,c_i$＝$\mathrm{Ext}^1$、関係 15 個＝$\partial f=0$。**$\mathbb{Z}/5$ 巡回対称の有限クイバー**で全幾何が符号化される＝あなたの「離散ネットワーク（位相配列）＋整合条件で空間を定義」の最も明示的な実例。ポテンシャルの臨界点（$\partial_{a_i}f=0$）＝§3 の $\mathrm{Crit}(f)$ 構図の局所版。

### 系 7.6（最終答）

**EN.** **Corollary 7.6** *The punctual partition function is given by*
$$Z(\mathscr{A}|P)(t) = Z(Q, f, \vec{1})(t, \dots, t).$$

Let us call this $Z(Q, f)(t)$.

So the 10 special points of the form $(1 : -1 : 0 : 0 : 0)$ contribute $Z(Q, f)(t)^{10}$ to the total partition function. There is a (complicated) box counting problem giving $Z(Q, f)(t)$ but we were not able to get a closed formula.

Generically, away from the 10 special points, we have
$$\mathscr{A} \cong M_{5 \times 5}(\mathscr{O}_X(\sqrt[5]{x_3}, \sqrt[5]{x_1})),\qquad \text{if } x_1 \ne 0, x_2 \ne 0.$$
So $\mathscr{A}$ is Morita equivalent to a commutative algebra. To study modules, ignore $M_{5 \times 5}$ up to rescaling the length by 5.

We get the final answer
$$\boxed{Z(X, \mathscr{A})(t) = Z(Q, f)(t)^{10} \cdot M(-t^5)^{-50}.}$$

**JA.** **系 7.6** *点状分配関数は*
$$Z(\mathscr{A}|P)(t) = Z(Q, f, \vec{1})(t, \dots, t)$$
*で与えられる。*

これを $Z(Q, f)(t)$ と呼ぶ。

すると $(1 : -1 : 0 : 0 : 0)$ の形の **10 個の特殊点**が、全分配関数に $Z(Q, f)(t)^{10}$ を寄与する。$Z(Q, f)(t)$ を与える（複雑な）箱数え上げ問題があるが、閉じた公式は得られなかった。

generic には、10 個の特殊点を離れたところで、
$$\mathscr{A} \cong M_{5 \times 5}(\mathscr{O}_X(\sqrt[5]{x_3}, \sqrt[5]{x_1})),\qquad (x_1 \ne 0,\ x_2 \ne 0).$$
よって $\mathscr{A}$ は可換代数と**森田同値（Morita equivalent）**である。加群を調べるには、長さを 5 でスケールし直して $M_{5 \times 5}$ を無視する。

最終的な答えを得る：
$$\boxed{Z(X, \mathscr{A})(t) = Z(Q, f)(t)^{10} \cdot M(-t^5)^{-50}.}$$

> 〔木原メモ〕**講演の到達点**：非可換クインティックの DT 分配関数 $Z(X,\mathscr{A})(t)=Z(Q,f)(t)^{10}\cdot M(-t^5)^{-50}$。構造が**「10 個の特殊点（クイバー局所・箱数え上げ）」×「generic 部分（森田同値で可換に帰着、$M(-t^5)^{-50}$）」**にきれいに分離。$M(-t^5)$＝MacMahon の $t\to t^5$、指数 $-50$＝50 点加群（$\binom{5}{2}\times5$）。**$10$ と $50$ と $5^4=625$ と $\mathbb{F}_5$ が全部効いた純組合せ的最終式**。「特殊点では離散クイバー、generic では連続的可換に崩れる」＝あなたの「対称点での退化と一般点での平坦化」と同じ分離構造。$Z(Q,f)$ が閉形式未到達＝この理論の最前線の未解決部分。

---

## 対訳まとめ — Behrend 講演の全体像

§1〜§7 を一望すると、この講演は次の一本の筋である：

- **§1–§2**：CY3 上の曲線の数え上げは、素朴な次元勘定では有限（2875 直線）に見えるが、Fermat のような対称点では連続族（50 鉛筆束）に**退化**する。これを救うのが**仮想基本類**。
- **§3–§4**：救済の鍵は、モジュライ空間が局所的に**ポテンシャル $f$ の臨界軌跡 $\mathrm{Crit}(f)$**（$-1$-シフト・シンプレクティック）であること。数え上げ＝**Behrend 函数 $\mu$ で重み付けた Euler 標数 $\chi(X,\mu)$**（特異ガウス・ボンネ）。$\mu$ は局所的・整数的・内在的。
- **§5**：これを**モチーフ化**して母関数に載せると、$\mathrm{Hilb}^n(\mathbb{C}^3)$＝**3 次元分割（plane partition）の母関数**、一般 CY3 では指数が $\chi$（クインティック $-200$）。
- **§6**：実際の理論は GW（有理数）／DT・GV（整数）／PT。**MNOP（Pardon）**＝全理論が離散基底 $x_{g,m}$ 上の値で決まり、$-q=e^{iu}$ で GW=PT。
- **§7**：幾何を全部捨てた**非可換クインティック $Q$**（$\mathbb{F}_5$・歪対称行列・5 頂点クイバー）でも、同じ DT 機械が走り、最終式 $Z(Q,f)^{10}\cdot M(-t^5)^{-50}$ に到達。

> 〔木原メモ・総括〕**この講演の背骨＝「連続的な数え上げ問題を、離散的・組合せ的・整数的な対象（3D 分割・クイバー・$\mathbb{F}_5$・MacMahon 関数）に還元して厳密に数える」**。あなたの離散基礎（整数・サイクル・格子・$\pi$-free、連続幾何は創発）と**方法論として強く共鳴**する。ただし規律として：これは確立済みの代数幾何（DT 理論）であり、あなたの中心投影／球面投影シリーズとは**別系統**。共鳴は「発想の隣接」として記録するに留め、既存論文へ概念注入はしない（[[feedback_foundational_buildup_discipline]]）。**6/6–6/7 の Behrend I・II の予習資料として完成。**
