# 関係位相差からの離散動力学と Noether 型局所保存則の導出
## ――零閉包・等振幅・有限位数回帰から連続保存則が現れるまで

**木原範昭**

## 要旨

前論文「物理学の対称性は、本当に最初から与える必要があるのか」では、背景時空、既成のゲージ群、連続力学を出発点に置かず、複素波の零閉包

$$
\sum_{a=1}^{M}X_a^2=0,
$$

有限位数回帰

$$
U^N=I,
$$

simplex 閉包、および自己無撞着条件という少数の公理から、相対論的構造および標準模型へ接続する豊かな対称構造を導出・整理した。一方、動力学への接続は主要な未解決課題として残った。

本研究は、その未解決問題を検討する過程で独立に得られた。出発点は既存の Noether 理論ではない。零閉包において時間軸を特権化せず、既導出の等振幅条件

$$
X_i=Ae^{i\phi_i},
\qquad
A=\mathrm{const.},
$$

の下で、隣接する二状態間の相対位相差が次の位相を決める最小自己写像を構成した。

等振幅複素波の二点間 current を

$$
\boxed{
J_{ij}
=
\operatorname{Im}(X_i^*X_j)
=
A^2\sin(\phi_j-\phi_i)
}
$$

と定義する。これは相対位相のみに依存し、辺の反転で符号が反転し、$2\pi$ 周期を持つ最小の非自明な反対称関係量である。

この current を用いて局所位相更新を

$$
\boxed{
\phi_i^{(n+1)}
=
\phi_i^{(n)}
+
\kappa\sum_{j\sim i}J_{ij}^{(n)}
}
$$

と置く。自己無撞着固定点

$$
X^{(n+1)}=X^{(n)}
$$

を要求すると、

$$
\kappa\sum_{j\sim i}J_{ij}
=
2\pi m_i,
\qquad
m_i\in\mathbb Z
$$

を得る。通常 sector $m_i=0$ では、

$$
\boxed{
\sum_{j\sim i}J_{ij}=0
}
$$

となり、これは頂点ごとの厳密な離散 continuity equation である。

さらに読出し分解能 $N$ を高くし、辺の位相差が小さい連続読出し極限では、

$$
J_{ij}
=
A^2\sin(\phi_j-\phi_i)
\longrightarrow
A^2\partial_\mu\phi\,\Delta q^\mu,
$$

したがって

$$
\boxed{
J^\mu=A^2\partial^\mu\phi
}
$$

を得る。離散 divergence の連続極限は、

$$
\boxed{
\partial_\mu J^\mu=0
}
$$

である。

本研究の中心結果は、Noether 型局所保存則を連続作用原理から出発せず、零閉包、等振幅、相対位相、simplex 閉包、自己無撞着、および有限回帰を持つ離散自己写像の連続読出し極限として導出したことである。さらに、$m_i\neq0$ sector は保存則の破れではなく、量子化された winding source / defect として現れる。

---

## 1. 研究の動機

前論文 [1] では、

$$
\sum_nX_n^2=0,
\qquad
U^N=I
$$

を中心とする零閉包、有限回帰、simplex 閉包、自己無撞着という少数の公理から、背景時空や既成の対称群を最初から与えずに、理論物理学が要求する多くの対称構造を導出・整理した。

そこに残った主要課題が動力学への接続である。

本研究の発見順序は明確である。最初に既存の離散 Noether 理論を調査して本模型へ移植したのではない。零閉包において時間軸を特権化する必要がないこと、また既導出の等振幅条件の下では状態変化の本体が振幅ではなく位相差にあることを再検討した。

その結果、

$$
\boxed{
\text{現在の関係位相差が次の位相を決める}
}
$$

という最小自己写像が、動力学の候補として自然に現れた。

この自己写像から局所離散 current を構成し、自己無撞着固定点を要求すると、局所 current の流入と流出が一致することが分かった。さらにその高分解能極限を取ると、偏微分形式

$$
\partial_\mu J^\mu=0
$$

が得られる。

この導出を得た後に、Noether の原論文、離散変分法、difference equations、離散場理論における Noether 型定理との関係を調査した。

したがって本研究の導出順序は、

$$
\boxed{
\text{零閉包}
\to
\text{等振幅}
\to
\text{相対位相 current}
\to
\text{自己無撞着離散写像}
\to
\text{局所離散保存則}
\to
N\to\infty
\to
\partial_\mu J^\mu=0
}
$$

である。

---

## 2. 基礎構造

### 2.1 零閉包

基礎式を

$$
\boxed{
C(X):=\sum_{a=1}^{M}X_a^2=0
}
\tag{1}
$$

とする。

この式では、特定の成分に時間・空間・内部量という名称を最初から与えない。

任意の一成分 $X_k$ を中心投影軸として右辺へ移せば、

$$
\boxed{
\sum_{a\ne k}X_a^2=-X_k^2
}
\tag{2}
$$

である。

$t$ を選べば $-t^2$ が閉包量として見え、$x$ を選べば $-x^2$ が同じ役割に見える。

従って、零閉包そのものには時間軸の特権性はない。

### 2.2 有限位数回帰

位相写像は、

$$
\boxed{
U^N=I
}
\tag{3}
$$

を満たす有限回帰として扱う。

ここで $N$ は基礎存在の絶対的な粒子数・波数ではない。頂点数としても読めるが、読出し側からは分解能として機能する。

$N$ を大きくすることは、新しい物理自由度を追加することではなく、同じ閉包構造をより細かく読むことに対応する。

### 2.3 等振幅

既導出の等振幅構造 [1] を用い、

$$
\boxed{
X_i=Ae^{i\phi_i},
\qquad
A=\mathrm{const.}
}
\tag{4}
$$

とする。

従って、状態変化の自由度は振幅ではなく位相に集約される。

---

## 3. 相対位相から current を作る

隣接する二状態 $i,j$ を考える。

$$
X_i=Ae^{i\phi_i},
\qquad
X_j=Ae^{i\phi_j}.
$$

二状態の積は、

$$
X_i^*X_j
=
A^2e^{i(\phi_j-\phi_i)}.
\tag{5}
$$

この虚部は、

$$
\operatorname{Im}(X_i^*X_j)
=
A^2\sin(\phi_j-\phi_i).
\tag{6}
$$

そこで辺 $i\to j$ に沿う離散 current を、

$$
\boxed{
J_{ij}
:=
\operatorname{Im}(X_i^*X_j)
=
A^2\sin(\phi_j-\phi_i)
}
\tag{7}
$$

と定義する。

この current は、

$$
J_{ji}
=
A^2\sin(\phi_i-\phi_j)
=
-J_{ij}
\tag{8}
$$

なので反対称である。

また、

$$
J_{ij}
$$

は絶対位相には依存せず、相対位相差

$$
\phi_j-\phi_i
$$

だけに依存する。

全位相を

$$
\phi_i\to\phi_i+\alpha
$$

と平行移動しても、

$$
J_{ij}
$$

は不変である。

従ってこの current は、無名な等振幅複素関係系に対する最小の局所関係量である。

---

## 4. なぜ $\sin(\Delta\phi)$ なのか

一般に、辺 current が相対位相差だけに依存するとして、

$$
J_{ij}=F(\phi_j-\phi_i)
\tag{9}
$$

と置く。

必要条件は、

1. 絶対位相に依存しない
2. $2\pi$ 周期を持つ
3. 辺反転で符号反転する
4. $\Delta\phi=0$ で current は 0
5. 最小の非自明 Fourier mode を採る

である。

条件 2,3 より $F$ は奇関数の周期関数であり、

$$
F(\Delta\phi)
=
\sum_{n=1}^{\infty}c_n\sin(n\Delta\phi)
\tag{10}
$$

と書ける。

このうち最小の非自明 mode は $n=1$ である。

従って、

$$
F(\Delta\phi)
\propto
\sin(\Delta\phi).
$$

比例係数は等振幅二乗 $A^2$ に吸収され、

$$
\boxed{
F(\Delta\phi)
=
A^2\sin(\Delta\phi)
}
\tag{11}
$$

を得る。

従って式 (7) は、相対位相のみ、反対称性、周期性、最小性を同時に満たす最小 current である。

---

## 5. 位相差が次の位相を決める最小自己写像

頂点 $i$ の次の位相を、その頂点へ接続された辺 current の総和で更新する。

$$
\boxed{
\phi_i^{(n+1)}
=
\phi_i^{(n)}
+
\kappa
\sum_{j\sim i}
J_{ij}^{(n)}
}
\tag{12}
$$

ここで $\kappa$ は位相更新の単位を合わせる共通係数であり、独立な物理フィッティングパラメータではない。読出し単位の選択に吸収できる。

式 (12) は、

$$
\boxed{
\text{現在の関係位相差}
\longrightarrow
\text{次の位相}
}
$$

という最小自己写像である。

重要なのは、絶対位相ではなく、関係位相のみが次状態を決めることである。

---

## 6. 自己無撞着固定点から離散保存則が出る

自己無撞着条件を、

$$
\boxed{
X_i^{(n+1)}=X_i^{(n)}
}
\tag{13}
$$

とする。

複素位相では、

$$
e^{i\phi_i^{(n+1)}}
=
e^{i\phi_i^{(n)}}
$$

なので、

$$
\phi_i^{(n+1)}-\phi_i^{(n)}
=
2\pi m_i,
\qquad
m_i\in\mathbb Z.
\tag{14}
$$

式 (12) を代入すると、

$$
\boxed{
\kappa
\sum_{j\sim i}J_{ij}
=
2\pi m_i
}
\tag{15}
$$

を得る。

### 6.1 通常 sector

局所的に winding を持たない sector、

$$
m_i=0
$$

では、

$$
\boxed{
\sum_{j\sim i}J_{ij}=0
}
\tag{16}
$$

である。

頂点 $i$ における離散 divergence を、

$$
\boxed{
(\operatorname{div}_dJ)_i
:=
\sum_{j\sim i}J_{ij}
}
\tag{17}
$$

と定義すれば、

$$
\boxed{
\operatorname{div}_dJ=0
}
\tag{18}
$$

である。

これは局所離散 continuity equation である。

### 6.2 winding sector

$m_i\ne0$ なら、

$$
\boxed{
(\operatorname{div}_dJ)_i
=
\frac{2\pi m_i}{\kappa}
}
\tag{19}
$$

となる。

従って $m_i\neq0$ は保存則の破れではなく、量子化された局所 source / defect として現れる。

---

## 7. simplex 上の離散 current

simplex 構造では、current $J_{ij}$ は向き付き辺上の 1-cochain として解釈できる。

辺反転で、

$$
J_{ji}=-J_{ij}
$$

なので、向き付き 1-cochain の条件を満たす。

頂点 $i$ での離散余微分は、

$$
(\delta J)_i
=
\sum_{j\sim i}J_{ij}.
\tag{20}
$$

従って通常 sector では、

$$
\boxed{
\delta J=0
}
\tag{21}
$$

である。

これは discrete co-closed current である。

---

## 8. 連続極限

読出し分解能を高くし、

$$
N\to\infty
$$

とする。

隣接頂点間の位相差が小さいとき、

$$
\Delta_\mu\phi
=
\phi(q+\Delta q^\mu)-\phi(q)
$$

に対して、

$$
\sin(\Delta_\mu\phi)
=
\Delta_\mu\phi
+
O((\Delta_\mu\phi)^3).
\tag{22}
$$

従って、

$$
J^\mu_N
=
\frac{A^2}{\Delta q^\mu}
\sin(\Delta_\mu\phi)
\tag{23}
$$

と定義すると、

$$
J^\mu_N
=
A^2
\frac{\Delta_\mu\phi}{\Delta q^\mu}
+
O((\Delta q^\mu)^2).
\tag{24}
$$

連続極限で、

$$
\boxed{
J^\mu
=
A^2\partial^\mu\phi
}
\tag{25}
$$

を得る。

---

## 9. 離散 divergence から偏微分保存則へ

離散 divergence を、

$$
(\operatorname{div}_dJ)(q)
=
\sum_\mu
\frac{
J^\mu(q+\Delta q^\mu)-J^\mu(q)
}{
\Delta q^\mu
}
\tag{26}
$$

とする。

通常 sector では、

$$
\boxed{
\operatorname{div}_dJ=0.
}
\tag{27}
$$

$N\to\infty$、$\Delta q^\mu\to0$ で、

$$
\frac{
J^\mu(q+\Delta q^\mu)-J^\mu(q)
}{
\Delta q^\mu
}
\longrightarrow
\partial_\mu J^\mu.
$$

従って、

$$
\boxed{
\partial_\mu J^\mu=0
}
\tag{28}
$$

を得る。

これが本論文の中心定理である。

---

## 10. Noether 型保存則の導出鎖

導出全体を一行で書けば、

$$
\boxed{
X_i=Ae^{i\phi_i}
\Longrightarrow
J_{ij}=A^2\sin(\phi_j-\phi_i)
\Longrightarrow
X' = X
\Longrightarrow
\operatorname{div}_dJ=0
\overset{N\to\infty}{\Longrightarrow}
\partial_\mu J^\mu=0
}
\tag{29}
$$

である。

重要なのは、この式に特権的時間軸が存在しないことである。

$\mu$ は任意の読出し方向であり、

$$
t
$$

はその一方向を更新軸として読むときに付けられる名称である。

---

## 11. Noether の定理との関係

Noether の原定理は、連続対称性を持つ変分問題と保存則を結び付ける [2]。

離散系についても、離散変分原理に基づく discrete Noether theorem [3]、difference equations に対する Noether 型定理 [4]、simplicial/cochain 離散場理論における保存則 [5] などが存在する。さらに、有限差分を局所連続性の基礎として離散時空上の保存則と連続極限を論じる研究もある [6]。

本研究は、これらの先行研究を出発点として構成したものではない。

式 (7)、式 (12)、式 (16)、式 (28) の導出を独立に得た後、その位置づけを確認するために先行研究を調査した。

既存の主要な離散 Noether 理論は、作用または離散 Lagrangian とその対称性から保存則を構成する。

これに対し本研究は、

$$
\boxed{
\text{零閉包}
+
\text{等振幅}
+
\text{相対位相}
+
\text{自己無撞着}
}
$$

から離散 current を構成し、その固定点条件から局所保存則を得る。

従って導出方向は異なる。

---

## 12. 時間とは何か

本模型で基礎的なのは、

$$
\boxed{
\text{任意軸間の位相差}
}
$$

である。

連続極限では、

$$
d\phi
=
\sum_\mu
\partial_\mu\phi\,dq^\mu.
\tag{30}
$$

そのうち一方向を観測上の更新方向として $t$ と読めば、

$$
\partial_t\phi
$$

が時間変化に見える。

しかし、基礎式では $t$ に特権性はない。

従って、

$$
\boxed{
\text{時間発展が位相を変えるのではなく、
位相写像の一方向を時間として読んでいる}
}
$$

と解釈できる。

---

## 13. 離散が基礎、連続は高分解能極限

本模型の順序は、

$$
\boxed{
\text{離散関係位相}
\to
\text{離散 current}
\to
\text{自己無撞着固定点}
\to
\text{局所離散保存則}
\to
N\to\infty
\to
\text{偏微分保存則}
}
\tag{31}
$$

である。

従って連続微分方程式は基礎ではない。

有限回帰可能な離散構造を高分解能で読んだ極限として現れる。

### 表1　離散構造と連続読出し

| 離散側 | 連続極限 | 物理的読出し |
|---|---|---|
| $U^N=I$ | 連続位相写像 | 有限回帰の高分解能表示 |
| $\phi_j-\phi_i$ | $\partial_\mu\phi\,dq^\mu$ | 位相差から位相勾配 |
| $J_{ij}=A^2\sin(\phi_j-\phi_i)$ | $J^\mu=A^2\partial^\mu\phi$ | 離散 current から局所 current |
| $\sum_jJ_{ij}=0$ | $\partial_\mu J^\mu=0$ | 離散 continuity から連続 continuity |
| $m_i\neq0$ | source / defect | winding による量子化 source |
| 任意読出し軸 | $t$ を一方向として選択 | 時間に基礎的特権性なし |

---

## 14. 動力学への接続

本研究で最も重要な帰結は、動力学の最小形が得られたことである。

$$
\boxed{
\phi_i^{(n+1)}
=
\phi_i^{(n)}
+
\kappa
\sum_{j\sim i}
A^2\sin(\phi_j^{(n)}-\phi_i^{(n)})
}
\tag{32}
$$

である。

この式は、

$$
\boxed{
\text{位相差が次の位相を決める}
}
$$

という関係動力学である。

従って動力学の問題は、外部から力を加える問題ではなく、

$$
\boxed{
\text{零閉包と有限回帰を保つ自己無撞着位相写像を分類すること}
}
$$

に変わる。

この分類から、連続極限における有効方程式を導くことが次段階の課題である。

---

## 15. winding source の帰結

式 (15) より、

$$
\kappa\sum_jJ_{ij}
=
2\pi m_i.
$$

従って、

$$
m_i=0
$$

なら source-free sector、

$$
m_i\neq0
$$

なら quantized source sector である。

連続極限では、

$$
\boxed{
\partial_\mu J^\mu
=
\rho_{\mathrm{defect}}
}
\tag{33}
$$

という形が現れ、

$$
\rho_{\mathrm{defect}}
$$

は離散 winding に由来する。

従って局所 source は保存則の破れではなく、非自明な位相回帰の読出しとして現れる。

---

## 16. 主結果

本稿の主結果は次の三点である。

第一に、等振幅複素関係系の相対位相から、

$$
\boxed{
J_{ij}=A^2\sin(\phi_j-\phi_i)
}
$$

という最小離散 current を構成した。

第二に、その current が次の位相を決める最小自己写像に自己無撞着固定点を要求すると、

$$
\boxed{
\operatorname{div}_dJ=0
}
$$

が得られることを示した。

第三に、その高分解能極限として、

$$
\boxed{
\partial_\mu J^\mu=0
}
$$

を導出した。

従って、

$$
\boxed{
\text{Noether 型局所保存則は、
離散自己無撞着位相動力学の連続極限として現れる}
}
$$

が本稿の中心結論である。

---

## 17. 結論

本研究では、前論文で残された動力学問題に対し、相対位相差が次の位相を決める最小自己写像を導入した。

これは外部から力学を追加する構成ではない。

零閉包、有限回帰、等振幅、simplex 閉包、自己無撞着という既存の構造から、

$$
J_{ij}
=
A^2\sin(\phi_j-\phi_i)
$$

という離散 current が自然に現れ、その current による局所位相更新の固定点条件から、

$$
\sum_jJ_{ij}=0
$$

という離散局所保存則が得られる。

さらに、

$$
N\to\infty
$$

の高分解能極限において、

$$
J^\mu=A^2\partial^\mu\phi
$$

および

$$
\partial_\mu J^\mu=0
$$

を得る。

従って、本模型では連続 Noether 型保存則は出発公理ではない。

$$
\boxed{
\text{離散関係位相}
\to
\text{離散 current}
\to
\text{自己無撞着}
\to
\text{局所保存}
\to
\text{連続極限}
}
$$

として現れる。

また、

$$
m_i\neq0
$$

sector は量子化された source / defect を与える。

本研究により、前論文で未解決だった「動力学への接続」は、単なる将来課題ではなく、

$$
\boxed{
\text{位相差が次の位相を決める自己無撞着離散写像}
}
$$

として具体的な数学形式を得た。

---

## 参考文献

1. 木原範昭，「物理学の対称性は、本当に最初から与える必要があるのか」，前論文.
2. E. Noether, “Invariante Variationsprobleme,” *Nachrichten von der Gesellschaft der Wissenschaften zu Göttingen, Mathematisch-Physikalische Klasse*, 1918, pp. 235–257.
3. J. E. Marsden and M. West, “Discrete Mechanics and Variational Integrators,” *Acta Numerica*, 10 (2001), 357–514. DOI: 10.1017/S096249290100006X.
4. V. A. Dorodnitsyn, “Noether-type theorems for difference equations,” *Applied Numerical Mathematics*, 39 (2001), 307–321. DOI: 10.1016/S0168-9274(00)00041-6.
5. M. Skopenkov, “Discrete Field Theory: Symmetries and Conservation Laws,” *Mathematical Physics, Analysis and Geometry*, 26 (2023), Article 19. DOI: 10.1007/s11040-023-09459-4.
6. S. R. Totorica, “Symmetries and conservation laws in discrete spacetime,” arXiv:2506.02119 (2025), preprint.

---

## 著者注

本研究の中心導出は、既存の離散 Noether 文献を参照して構成したものではない。前論文で残された動力学問題を、零閉包、等振幅、相対位相、任意軸への写像、自己無撞着から検討する過程で独立に得た。

その後、既存研究との重複と差異を確認するために、Noether の原論文、離散変分理論、difference equations、離散場理論に関する先行研究を調査した。
