# 離散零閉包からの Noether 保存則と関係位相動力学
## ――零閉包を厳密保存する離散自己写像と $N\to\infty$ 連続場方程式

**木原範昭**

## 要旨

前論文「物理学の対称性は、本当に最初から与える必要があるのか」では、背景時空や既成の対称群を最初から置かず、

$$
\sum_{a=1}^{M}X_a^2=0,
\qquad
U^N=I,
$$

simplex 閉包および自己無撞着という少数の条件から、多数の対称構造を導出・整理した。しかし、Noether の定理に対応する局所保存則と、次状態を決める動力学への接続は残った。

本論文では、等振幅

$$
X_i=Ae^{i\phi_i},
\qquad
A=\mathrm{const.}
$$

の下で、関係位相差から

$$
J_{ij}
=
A^2\sin(\phi_j-\phi_i)
$$

という向き付き離散 current を得る。さらに、

$$
S_N[\phi]
=
-A^2\sum_{\langle ij\rangle}
\cos(\phi_j-\phi_i)
$$

という離散作用を構成すると、

$$
\frac{\partial S_N}{\partial\phi_i}
=
-\sum_{j\sim i}J_{ij}
$$

であり、停留条件は

$$
\sum_{j\sim i}J_{ij}=0
$$

となる。

動力学については、単純な位相更新をそのまま用いるのではなく、零閉包

$$
C(\phi)
=
\sum_i e^{2i\phi_i}
=
0
$$

を厳密に保存する制約付き自己写像を構成する。零閉包は実二条件

$$
C_R(\phi)=\sum_i\cos2\phi_i=0,
\qquad
C_I(\phi)=\sum_i\sin2\phi_i=0
$$

と同値である。したがって、未制約の関係位相力

$$
F_i(\phi)
=
\sum_{j\sim i}\sin(\phi_j-\phi_i)
$$

を零閉包多様体の接空間へ射影し、

$$
\dot\phi
=
P_\phi F(\phi)
$$

とすれば、連続な写像パラメータに沿って零閉包は厳密に保存される。有限反復では、この接ベクトルを零閉包多様体上へ戻す retraction を用い、

$$
\boxed{
\phi^{(n+1)}
=
R_{\phi^{(n)}}
\!\left(
\eta P_{\phi^{(n)}}F(\phi^{(n)})
\right)
}
$$

と定義する。これにより、

$$
\phi^{(n)}\in\mathcal Z_N
\Longrightarrow
\phi^{(n+1)}\in\mathcal Z_N
$$

が各反復で成立する。

読出し分解能 $N$ を高くし、隣接間隔 $h_N\to0$ とすると、離散 current は

$$
J_N^\mu
\longrightarrow
A^2\partial^\mu\phi
$$

へ移り、離散 continuity equation は

$$
\boxed{
\partial_\mu J^\mu=0
}
$$

へ収束する。

さらに離散作用自身も連続極限で

$$
S_{\mathrm{cont}}
=
\frac{A^2}{2}
\int
g^{\mu\nu}
\partial_\mu\phi
\partial_\nu\phi
\,d^dq
$$

へ移るため、Euler–Lagrange 方程式は

$$
\boxed{
\partial_\mu
\left(
A^2g^{\mu\nu}\partial_\nu\phi
\right)
=
0
}
$$

となる。$A$ と $g^{\mu\nu}$ が一定なら、

$$
\boxed{
\Box_g\phi=0
}
$$

である。

したがって本論文では、離散関係位相系から、零閉包を厳密保存する自己写像、Noether 型局所保存則、および標準的な偏微分場方程式への連続極限を一つの導出鎖として与える。

---

## 1. 研究の動機

前論文 [1] では、多数の対称性を少数の閉包条件から導出した。

しかし、二つの問題が残った。

第一に、

$$
\boxed{
\text{対称性は導出したが、
Noether 保存則は導出していなかった}
}
$$

ことである。

第二に、

$$
\boxed{
\text{自己無撞着構造から、
次状態を決める動力学へ接続していなかった}
}
$$

ことである。

本研究はこの二点を検討した。

検討の結果、Noether 型保存則は有限 $N$ の離散 current 保存則の $N\to\infty$ 連続極限として現れた。さらに同じ関係位相 current が、次位相を決める自己写像の生成量になることが分かった。

ただし、単純な更新

$$
\phi_i^{(n+1)}
=
\phi_i^{(n)}
+
\eta F_i
$$

では、零閉包の厳密保存は自動ではない。

したがって本論文では、零閉包を保存する接空間射影と有限 retraction を導入し、

$$
\boxed{
\mathcal F:
\mathcal Z_N\to\mathcal Z_N
}
$$

を満たす動力学として閉じる。

---

## 2. 基礎状態と零閉包

### 2.1 等振幅状態

$$
\boxed{
X_i=Ae^{i\phi_i},
\qquad
A=\mathrm{const.}
}
\tag{1}
$$

である。

### 2.2 零閉包

第一公理は、

$$
\sum_iX_i^2=0.
$$

式 (1) を代入すると、

$$
A^2\sum_i e^{2i\phi_i}=0.
$$

$A\neq0$ なので、

$$
\boxed{
C(\phi)
:=
\sum_i e^{2i\phi_i}
=
0
}
\tag{2}
$$

である。

実部・虚部に分けると、

$$
\boxed{
C_R(\phi)
=
\sum_i\cos2\phi_i
=
0
}
\tag{3}
$$

$$
\boxed{
C_I(\phi)
=
\sum_i\sin2\phi_i
=
0
}
\tag{4}
$$

である。

従って許容位相集合は、

$$
\boxed{
\mathcal Z_N
=
\left\{
\phi
\in
(\mathbb R/2\pi\mathbb Z)^M
\ \middle|\
C_R(\phi)=0,
\ C_I(\phi)=0,
\ U^N=I
\right\}.
}
\tag{5}
$$

---

## 3. 関係位相 current

隣接状態 $i,j$ について、

$$
X_i^*X_j
=
A^2e^{i(\phi_j-\phi_i)}.
\tag{6}
$$

向き付き辺 current を、

$$
\boxed{
J_{ij}
=
\operatorname{Im}(X_i^*X_j)
=
A^2\sin(\phi_j-\phi_i)
}
\tag{7}
$$

と置く。

辺反転で、

$$
\boxed{
J_{ji}=-J_{ij}.
}
\tag{8}
$$

また、

$$
\phi_i\mapsto\phi_i+\alpha
$$

に対して $J_{ij}$ は不変である。

---

## 4. 離散作用と離散 Noether 保存則

各辺の実部から、

$$
s_{ij}
=
-\operatorname{Re}(X_i^*X_j)
=
-A^2\cos(\phi_j-\phi_i)
\tag{9}
$$

を得る。

全作用は、

$$
\boxed{
S_N[\phi]
=
-A^2
\sum_{\langle ij\rangle}
\cos(\phi_j-\phi_i).
}
\tag{10}
$$

全体位相変換

$$
\phi_i\mapsto\phi_i+\alpha
$$

では位相差が変わらないため、

$$
\boxed{
S_N[\phi+\alpha]=S_N[\phi].
}
\tag{11}
$$

頂点位相 $\phi_i$ で変分すると、

$$
\frac{\partial S_N}{\partial\phi_i}
=
-A^2
\sum_{j\sim i}
\sin(\phi_j-\phi_i).
$$

従って、

$$
\boxed{
\frac{\partial S_N}{\partial\phi_i}
=
-\sum_{j\sim i}J_{ij}.
}
\tag{12}
$$

停留条件

$$
\frac{\partial S_N}{\partial\phi_i}=0
$$

は、

$$
\boxed{
\sum_{j\sim i}J_{ij}=0
}
\tag{13}
$$

となる。

離散 divergence を、

$$
\boxed{
(\operatorname{div}_dJ)_i
=
\sum_{j\sim i}J_{ij}
}
\tag{14}
$$

と定義すれば、

$$
\boxed{
\operatorname{div}_dJ=0.
}
\tag{15}
$$

---

## 5. 未制約の最小関係位相力

作用の負勾配を、

$$
F_i(\phi)
:=
-\frac{1}{A^2}
\frac{\partial S_N}{\partial\phi_i}
$$

と定義する。

式 (12) より、

$$
\boxed{
F_i(\phi)
=
\sum_{j\sim i}
\sin(\phi_j-\phi_i).
}
\tag{16}
$$

これは、

$$
\boxed{
\text{現在の関係位相差が次の位相変化を決める}
}
$$

という最小局所生成量である。

しかし、

$$
\phi^{(n+1)}
=
\phi^{(n)}
+
\eta F(\phi^{(n)})
$$

だけでは、一般には、

$$
C(\phi^{(n+1)})=0
$$

は保証されない。

ここを次節で閉じる。

---

## 6. 零閉包多様体の接空間

零閉包の実二条件は式 (3)、式 (4) である。

その勾配は、

$$
\frac{\partial C_R}{\partial\phi_i}
=
-2\sin2\phi_i,
$$

$$
\frac{\partial C_I}{\partial\phi_i}
=
2\cos2\phi_i.
$$

従って、

$$
\boxed{
g_R(\phi)
=
(-2\sin2\phi_1,\ldots,-2\sin2\phi_M)^T
}
\tag{17}
$$

$$
\boxed{
g_I(\phi)
=
(2\cos2\phi_1,\ldots,2\cos2\phi_M)^T.
}
\tag{18}
$$

位相変化ベクトル $v\in\mathbb R^M$ が零閉包多様体の接ベクトルである条件は、

$$
g_R^Tv=0,
\qquad
g_I^Tv=0.
\tag{19}
$$

これは直接、

$$
\frac{dC_R}{ds}=0,
\qquad
\frac{dC_I}{ds}=0
$$

を意味する。

---

## 7. 未制約力を接空間へ射影する

$M\times2$ 行列

$$
\boxed{
G(\phi)
=
\begin{pmatrix}
|&|\\
g_R&g_I\\
|&|
\end{pmatrix}
}
\tag{20}
$$

を定義する。

$g_R,g_I$ が独立な通常点では、

$$
G^TG
$$

は可逆である。

接空間への直交射影は、

$$
\boxed{
P_\phi
=
I
-
G(G^TG)^{-1}G^T.
}
\tag{21}
$$

実際、

$$
G^TP_\phi
=
G^T
-
G^TG(G^TG)^{-1}G^T
=
0.
$$

従って、

$$
\boxed{
G^TP_\phi=0.
}
\tag{22}
$$

未制約力 $F$ から、零閉包を壊さない接方向力を、

$$
\boxed{
F_T(\phi)
=
P_\phi F(\phi)
}
\tag{23}
$$

と定義する。

すると、

$$
g_R^TF_T=0,
\qquad
g_I^TF_T=0.
\tag{24}
$$

---

## 8. 連続自己写像パラメータで零閉包が厳密保存される

ここではまだ物理時間を導入しない。

単なる自己写像パラメータ $s$ を用い、

$$
\boxed{
\frac{d\phi}{ds}
=
F_T(\phi)
=
P_\phi F(\phi)
}
\tag{25}
$$

とする。

すると、

$$
\frac{dC_R}{ds}
=
\nabla C_R^T
\frac{d\phi}{ds}
=
g_R^TP_\phi F
=
0.
$$

同様に、

$$
\frac{dC_I}{ds}
=
g_I^TP_\phi F
=
0.
$$

従って、

$$
\boxed{
\frac{dC_R}{ds}
=
\frac{dC_I}{ds}
=
0.
}
\tag{26}
$$

初期状態が零閉包を満たせば、

$$
C_R(\phi(0))=C_I(\phi(0))=0
$$

なので、任意の $s$ で、

$$
\boxed{
C_R(\phi(s))
=
C_I(\phi(s))
=
0.
}
\tag{27}
$$

すなわち、

$$
\boxed{
\sum_iX_i(s)^2=0
}
\tag{28}
$$

が厳密に保存される。

これで連続な自己写像パラメータにおける零閉包保存は閉じた。

---

## 9. 有限反復で零閉包を厳密保存する

Euler 更新、

$$
\phi+\eta F_T
$$

は接方向には正しいが、有限 $\eta$ では多様体から $O(\eta^2)$ だけ外れる。

従って有限反復では retraction を用いる。

接方向仮更新を、

$$
\widetilde\phi
=
\phi^{(n)}
+
\eta
P_{\phi^{(n)}}F(\phi^{(n)})
\tag{29}
$$

とする。

次に、二つの法線方向 $g_R,g_I$ を用いて、

$$
\phi^{(n+1)}
=
\widetilde\phi
+
\lambda_Rg_R(\widetilde\phi)
+
\lambda_Ig_I(\widetilde\phi)
\tag{30}
$$

と置く。

$\lambda_R,\lambda_I$ を、

$$
C_R(\phi^{(n+1)})=0,
$$

$$
C_I(\phi^{(n+1)})=0
\tag{31}
$$

の二式から決める。

通常点では implicit function theorem により、十分小さい $\eta$ に対して局所一意な

$$
(\lambda_R,\lambda_I)
$$

が存在する。

この操作を、

$$
R_\phi(v)
$$

と書けば、

$$
\boxed{
\phi^{(n+1)}
=
R_{\phi^{(n)}}
\left(
\eta
P_{\phi^{(n)}}F(\phi^{(n)})
\right)
}
\tag{32}
$$

である。

定義から、

$$
\boxed{
\phi^{(n)}\in\mathcal Z_N
\Longrightarrow
\phi^{(n+1)}\in\mathcal Z_N.
}
\tag{33}
$$

従って、

$$
\boxed{
\mathcal F_N:
\mathcal Z_N
\to
\mathcal Z_N
}
\tag{34}
$$

が有限反復ごとに厳密に成立する。

これが零閉包を保存する離散動力学である。

---

## 10. 自己無撞着条件

自己無撞着固定点は、

$$
\boxed{
\phi_*
=
\mathcal F_N(\phi_*)
\quad
(\bmod 2\pi)
}
\tag{35}
$$

である。

$\phi_*$ は既に $\mathcal Z_N$ 上にあるので零閉包を満たす。

固定点では接方向更新が消えるため、

$$
\boxed{
P_{\phi_*}F(\phi_*)=0.
}
\tag{36}
$$

すなわち、

$$
F(\phi_*)
\in
\operatorname{span}
\{g_R(\phi_*),g_I(\phi_*)\}.
\tag{37}
$$

従って、零閉包を保った許容変分に対して作用は停留する。

制約付き停留条件は Lagrange multiplier を用いて、

$$
\boxed{
\frac{\partial S_N}{\partial\phi_i}
+
\lambda_R
\frac{\partial C_R}{\partial\phi_i}
+
\lambda_I
\frac{\partial C_I}{\partial\phi_i}
=
0.
}
\tag{38}
$$

これが自己無撞着と零閉包を同時に満たす有限 $N$ の離散 Euler–Lagrange 方程式である。

---

## 11. 制約なし sector と局所 current 保存

零閉包に対する法線反力が局所 current sector に寄与しない、すなわち、

$$
\lambda_R=\lambda_I=0
$$

の sector では、式 (38) は、

$$
\frac{\partial S_N}{\partial\phi_i}=0
$$

となる。

従って、

$$
\boxed{
\sum_{j\sim i}J_{ij}=0.
}
\tag{39}
$$

これは離散 continuity equation である。

一般の制約付き sector では、

$$
\sum_{j\sim i}J_{ij}
=
\lambda_R
\frac{\partial C_R}{\partial\phi_i}
+
\lambda_I
\frac{\partial C_I}{\partial\phi_i}.
\tag{40}
$$

従って、零閉包制約力を明示的 source として分離できる。

---

## 12. 読出し分解能 $N$ と格子間隔

正規化された読出し長を $L$ とする。

分解能 $N$ では、

$$
\boxed{
h_N=\frac{L}{N}.
}
\tag{41}
$$

従って、

$$
\boxed{
N\to\infty
\quad\Longleftrightarrow\quad
h_N\to0.
}
\tag{42}
$$

$N$ は存在する要素の絶対個数ではなく、同じ構造を読む分解能である。

---

## 13. 離散 current density の連続極限

方向 $\mu$ の辺 current density を、

$$
\boxed{
J_N^\mu
\left(
q+\frac{h_N}{2}e_\mu
\right)
=
\frac{A^2}{h_N}
\sin
\left[
\phi(q+h_Ne_\mu)-\phi(q)
\right].
}
\tag{43}
$$

とする。

Taylor 展開で、

$$
\phi(q+h_Ne_\mu)
-
\phi(q)
=
h_N\partial_\mu\phi
+
\frac{h_N^2}{2}
\partial_\mu^2\phi
+
O(h_N^3).
\tag{44}
$$

また、

$$
\sin u
=
u
-
\frac{u^3}{6}
+
O(u^5).
\tag{45}
$$

従って、

$$
J_N^\mu
=
A^2\partial_\mu\phi
+
O(h_N).
\tag{46}
$$

よって、

$$
\boxed{
N\to\infty
\quad\Longrightarrow\quad
J_\mu
=
A^2\partial_\mu\phi.
}
\tag{47}
$$

読出し計量 $g_{\mu\nu}$ を選べば、

$$
\boxed{
J^\mu
=
A^2g^{\mu\nu}\partial_\nu\phi.
}
\tag{48}
$$

---

## 14. 離散 divergence から連続 continuity equation

頂点 $q$ における正方向流出と負方向流入を組にすると、

$$
\sum_\mu
\left[
J_N^\mu
\left(
q+\frac{h_N}{2}e_\mu
\right)
-
J_N^\mu
\left(
q-\frac{h_N}{2}e_\mu
\right)
\right]
=0.
\tag{49}
$$

$h_N$ で割ると、

$$
\boxed{
\sum_\mu
\frac{
J_N^\mu(q+\tfrac12h_Ne_\mu)
-
J_N^\mu(q-\tfrac12h_Ne_\mu)
}{
h_N
}
=
0.
}
\tag{50}
$$

中心差分の Taylor 展開は、

$$
\frac{
J_N^\mu(q+\tfrac12h_Ne_\mu)
-
J_N^\mu(q-\tfrac12h_Ne_\mu)
}{
h_N
}
=
\partial_\mu J^\mu
+
O(h_N^2).
\tag{51}
$$

従って、

$$
\boxed{
N\to\infty
\quad\Longrightarrow\quad
\partial_\mu J^\mu=0.
}
\tag{52}
$$

これが連続 Noether 型保存則である。

---

## 15. 離散作用から連続作用へ

離散作用は、

$$
S_N
=
-A^2
\sum_{\langle ij\rangle}
\cos(\phi_j-\phi_i).
$$

小位相差で、

$$
\cos u
=
1-\frac{u^2}{2}
+
O(u^4).
\tag{53}
$$

定数項を除けば、

$$
S_N
\sim
\frac{A^2}{2}
\sum_{\langle ij\rangle}
(\phi_j-\phi_i)^2.
\tag{54}
$$

方向 $\mu$ について、

$$
\phi(q+h_Ne_\mu)-\phi(q)
=
h_N\partial_\mu\phi
+
O(h_N^2).
$$

従って、cell volume $h_N^d$ と計量重みを含めて Riemann 和を取ると、

$$
\boxed{
S_N
\longrightarrow
S_{\mathrm{cont}}
=
\frac{A^2}{2}
\int
g^{\mu\nu}
\partial_\mu\phi
\partial_\nu\phi
\,d^dq.
}
\tag{55}
$$

Lagrangian density は、

$$
\boxed{
\mathcal L
=
\frac{A^2}{2}
g^{\mu\nu}
\partial_\mu\phi
\partial_\nu\phi.
}
\tag{56}
$$

である。

---

## 16. 標準的な偏微分場方程式への接続

Euler–Lagrange 方程式は、

$$
\frac{\partial\mathcal L}{\partial\phi}
-
\partial_\mu
\left(
\frac{\partial\mathcal L}
{\partial(\partial_\mu\phi)}
\right)
=
0.
\tag{57}
$$

$\mathcal L$ は $\phi$ 自身に依存しないので、

$$
\frac{\partial\mathcal L}{\partial\phi}=0.
$$

また、

$$
\frac{\partial\mathcal L}
{\partial(\partial_\mu\phi)}
=
A^2g^{\mu\nu}\partial_\nu\phi
=
J^\mu.
\tag{58}
$$

従って、

$$
\boxed{
\partial_\mu J^\mu
=
\partial_\mu
\left(
A^2g^{\mu\nu}\partial_\nu\phi
\right)
=
0.
}
\tag{59}
$$

$A$ と $g^{\mu\nu}$ が一定なら、

$$
A^2
g^{\mu\nu}
\partial_\mu\partial_\nu\phi
=
0.
$$

従って、

$$
\boxed{
\Box_g\phi=0.
}
\tag{60}
$$

Euclid 読出しなら、

$$
\boxed{
\nabla^2\phi=0.
}
\tag{61}
$$

Lorentz 読出しなら、

$$
\boxed{
\Box\phi
=
-\partial_t^2\phi
+
\nabla^2\phi
=
0.
}
\tag{62}
$$

従って、$N\to\infty$ の連続近似では、離散関係位相系は標準的な massless scalar / phase field の偏微分方程式へ接続する。

ここで Lorentz 符号は基礎に時間を特権化して置いたのではない。読出し時に選ばれた実形式の符号から現れる。

---

## 17. 制約付き連続場方程式

零閉包を連続極限でも明示的に保持する場合は、作用に Lagrange multiplier を加える。

$$
\boxed{
S_{\mathrm{c}}
=
S_{\mathrm{cont}}
+
\int
\left[
\Lambda_R C_R
+
\Lambda_I C_I
\right]
d^dq.
}
\tag{63}
$$

変分すると、

$$
\boxed{
\partial_\mu
\left(
A^2g^{\mu\nu}\partial_\nu\phi_i
\right)
=
\Lambda_R
\frac{\partial C_R}{\partial\phi_i}
+
\Lambda_I
\frac{\partial C_I}{\partial\phi_i}.
}
\tag{64}
$$

すなわち、

$$
\boxed{
A^2\Box_g\phi_i
=
-2\Lambda_R\sin2\phi_i
+
2\Lambda_I\cos2\phi_i.
}
\tag{65}
$$

制約反力が消える sector では、

$$
\boxed{
\Box_g\phi_i=0
}
\tag{66}
$$

へ戻る。

従って一般形は制約付き wave / Laplace 型方程式であり、標準的自由場方程式はその零制約反力 sector として現れる。

---


## 18. 局所ゲージ構造の導出

本節では、前節までの global $U(1)$ 位相構造を局所化し、辺 connection、面 curvature、共変微分、および Maxwell / Yang–Mills 型連続極限までを導出する。

### 18.1 局所位相原点の無名性

global 位相変換

$$
\phi_i\mapsto\phi_i+\alpha
$$

では、関係位相差

$$
\phi_j-\phi_i
$$

は不変である。

しかし、各頂点の位相原点を独立に選び直す

$$
\boxed{
\phi_i\mapsto\phi_i+\alpha_i
}
\tag{67}
$$

を許すと、

$$
\phi_j-\phi_i
\mapsto
\phi_j-\phi_i+\alpha_j-\alpha_i
$$

となり、裸の位相差は不変ではない。

従って、異なる頂点間の位相を比較するためには、辺に補償量 $\theta_{ij}$ が必要である。

局所的に観測可能な関係位相を、

$$
\boxed{
\Delta_{ij}^{(g)}
=
\phi_j-\phi_i-\theta_{ij}
}
\tag{68}
$$

と定義する。

式 (67) の局所位相変換に対して $\Delta_{ij}^{(g)}$ を不変に保つには、

$$
\boxed{
\theta_{ij}
\mapsto
\theta_{ij}
+
\alpha_j-\alpha_i
}
\tag{69}
$$

でなければならない。

これは $U(1)$ connection の離散ゲージ変換則である。

従って、辺 connection は外部から追加された変数ではなく、

$$
\boxed{
\text{局所位相原点の無名性}
\Longrightarrow
\text{異なる頂点を比較する辺 connection}
}
\tag{70}
$$

として必要になる。

### 18.2 link variable

辺 connection を指数写像して、

$$
\boxed{
U_{ij}
=
e^{-i\theta_{ij}}
}
\tag{71}
$$

とする。

局所変換

$$
X_i
\mapsto
e^{i\alpha_i}X_i
$$

に対して、

$$
U_{ij}
\mapsto
e^{i\alpha_i}
U_{ij}
e^{-i\alpha_j}
$$

であれば、

$$
X_i^*U_{ij}X_j
$$

は不変である。

実際、

$$
X_i^*U_{ij}X_j
\mapsto
e^{-i\alpha_i}X_i^*
e^{i\alpha_i}U_{ij}e^{-i\alpha_j}
e^{i\alpha_j}X_j
=
X_i^*U_{ij}X_j.
$$

従って gauge-covariant 離散 current は、

$$
\boxed{
J_{ij}^{(g)}
=
\operatorname{Im}
\left(
X_i^*U_{ij}X_j
\right)
}
\tag{72}
$$

である。

等振幅

$$
X_i=Ae^{i\phi_i}
$$

を代入すると、

$$
\boxed{
J_{ij}^{(g)}
=
A^2
\sin
\left(
\phi_j-\phi_i-\theta_{ij}
\right)
}
\tag{73}
$$

となる。

### 18.3 gauge-covariant 離散作用

同様に実部から、

$$
\boxed{
S_N^{(g)}
=
-A^2
\sum_{\langle ij\rangle}
\cos
\left(
\phi_j-\phi_i-\theta_{ij}
\right)
}
\tag{74}
$$

を得る。

式 (68)、(69) より各辺の引数は gauge invariant だから、

$$
\boxed{
S_N^{(g)}
\text{ は局所 }U(1)\text{ 変換に厳密に不変}
}
\tag{75}
$$

である。

$\phi_i$ で変分すると、

$$
\frac{\partial S_N^{(g)}}{\partial\phi_i}
=
-\sum_{j\sim i}
J_{ij}^{(g)}.
$$

従って停留条件は、

$$
\boxed{
\sum_{j\sim i}
J_{ij}^{(g)}
=
0
}
\tag{76}
$$

である。

これは局所 gauge-covariant continuity equation の離散形である。

### 18.4 simplex 面から curvature が出る

三角形 $i\to j\to k\to i$ を考える。

閉路に沿う connection の和を、

$$
\boxed{
\Theta_{ijk}
=
\theta_{ij}
+
\theta_{jk}
+
\theta_{ki}
}
\tag{77}
$$

とする。

局所変換では、

$$
\theta_{ij}
\mapsto
\theta_{ij}+\alpha_j-\alpha_i,
$$

$$
\theta_{jk}
\mapsto
\theta_{jk}+\alpha_k-\alpha_j,
$$

$$
\theta_{ki}
\mapsto
\theta_{ki}+\alpha_i-\alpha_k.
$$

従って、

$$
\Theta_{ijk}
\mapsto
\Theta_{ijk}
+
(\alpha_j-\alpha_i)
+
(\alpha_k-\alpha_j)
+
(\alpha_i-\alpha_k)
=
\Theta_{ijk}.
$$

ゆえに、

$$
\boxed{
\Theta_{ijk}
\text{ は gauge invariant}
}
\tag{78}
$$

である。

従って simplex の階層は、

$$
\boxed{
\text{頂点}
\to
\text{位相場}
}
$$

$$
\boxed{
\text{辺}
\to
\text{connection}
}
$$

$$
\boxed{
\text{面}
\to
\text{curvature}
}
\tag{79}
$$

と対応する。

### 18.5 plaquette action

面 curvature に対する最小 $2\pi$ 周期作用を、

$$
\boxed{
S_F
=
\beta
\sum_p
\left(
1-\cos\Theta_p
\right)
}
\tag{80}
$$

とする。

$\Theta_p$ は gauge invariant なので、$S_F$ も gauge invariant である。

小さい plaquette では、

$$
1-\cos\Theta_p
=
\frac12\Theta_p^2
+
O(\Theta_p^4).
\tag{81}
$$

従って連続極限では curvature 二乗作用へ移る。

### 18.6 連続極限で共変微分が出る

読出し間隔を $h_N$ とする。

方向 $\mu$ の辺について、

$$
\theta_{ij}
=
gh_NA_\mu(q)
+
O(h_N^2)
\tag{82}
$$

と置く。

また、

$$
\phi(q+h_Ne_\mu)-\phi(q)
=
h_N\partial_\mu\phi
+
O(h_N^2).
$$

従って、

$$
\frac{
\phi_j-\phi_i-\theta_{ij}
}{
h_N
}
\longrightarrow
\partial_\mu\phi-gA_\mu.
$$

よって gauge current は、

$$
\boxed{
J^\mu_{(g)}
=
A^2
\left(
\partial^\mu\phi
-
gA^\mu
\right)
}
\tag{83}
$$

となる。

複素場

$$
\Psi=Ae^{i\phi}
$$

について、

$$
D_\mu
=
\partial_\mu
-
igA_\mu
$$

とすれば、

$$
D_\mu\Psi
=
iAe^{i\phi}
\left(
\partial_\mu\phi-gA_\mu
\right).
$$

従って、

$$
\boxed{
J^\mu_{(g)}
\propto
\operatorname{Im}
\left(
\Psi^*D^\mu\Psi
\right)
}
\tag{84}
$$

であり、標準的 covariant derivative へ一致する。

### 18.7 面 curvature の連続極限

$\mu\nu$ 面の小矩形 plaquette を考える。

閉路和は、

$$
\Theta_{\mu\nu}
=
\theta_\mu(q)
+
\theta_\nu(q+h_Ne_\mu)
-
\theta_\mu(q+h_Ne_\nu)
-
\theta_\nu(q).
\tag{85}
$$

式 (82) を Taylor 展開すると、

$$
\Theta_{\mu\nu}
=
gh_N^2
\left(
\partial_\mu A_\nu
-
\partial_\nu A_\mu
\right)
+
O(h_N^3).
$$

従って、

$$
\boxed{
\frac{\Theta_{\mu\nu}}{gh_N^2}
\longrightarrow
F_{\mu\nu}
=
\partial_\mu A_\nu
-
\partial_\nu A_\mu
}
\tag{86}
$$

である。

### 18.8 Maxwell action の極限

式 (80) と式 (81) より、

$$
S_F
\sim
\frac{\beta}{2}
\sum_p
\Theta_p^2.
$$

式 (86) を用いると、

$$
\Theta_p^2
=
g^2h_N^4
F_{\mu\nu}F^{\mu\nu}
+
O(h_N^5).
$$

cell volume と結合定数規格化を含めて Riemann 和を取れば、

$$
\boxed{
S_F
\longrightarrow
-\frac14
\int
F_{\mu\nu}F^{\mu\nu}
\,d^dq
}
\tag{87}
$$

へ接続する。

従って、

$$
\boxed{
\text{辺 connection}
+
\text{面 curvature}
\Longrightarrow
\text{Maxwell gauge dynamics}
}
\tag{88}
$$

である。

### 18.9 非可換化

頂点状態を、

$$
\Psi_i\in\mathbb C^r
$$

とし、局所内部基底変換を、

$$
\Psi_i
\mapsto
G_i\Psi_i,
\qquad
G_i\in SU(r)
\tag{89}
$$

とする。

異なる頂点の状態を比較するには、辺 transporter

$$
U_{ij}\in SU(r)
$$

が必要であり、

$$
\boxed{
U_{ij}
\mapsto
G_iU_{ij}G_j^{-1}
}
\tag{90}
$$

でなければならない。

閉路積

$$
\boxed{
W_{ijk}
=
U_{ij}U_{jk}U_{ki}
}
\tag{91}
$$

は、

$$
W_{ijk}
\mapsto
G_iW_{ijk}G_i^{-1}
$$

だから、

$$
\boxed{
\operatorname{Tr}W_{ijk}
\text{ は gauge invariant}
}
\tag{92}
$$

である。

連続極限で、

$$
U_{ij}
=
\exp
\left(
igh_NA_\mu^aT^a
\right)
$$

と置くと、plaquette 積の Baker–Campbell–Hausdorff 展開から、

$$
\boxed{
F_{\mu\nu}
=
\partial_\mu A_\nu
-
\partial_\nu A_\mu
-
ig[A_\mu,A_\nu]
}
\tag{93}
$$

を得る。

従って Yang–Mills の非可換項

$$
[A_\mu,A_\nu]
$$

は、有限 transporter の閉路積から現れる。

### 18.10 本論文でここまで導出できたこと

本節で得たのは、

$$
\boxed{
\text{局所位相無名性}
\to
\text{辺 connection}
\to
\text{面 curvature}
\to
\text{covariant derivative}
\to
\text{Maxwell / Yang--Mills 型連続極限}
}
\tag{94}
$$

である。

従って、標準模型との接続は前版より一段強化される。


## 19. 標準模型との整合性と接続範囲

本節では、本論文で得た連続極限が標準模型と矛盾しないか、また「標準模型へ接続した」と言える範囲を明確に検証する。

### 18.1 標準模型が要求する最小構造

標準模型は、

$$
\boxed{
G_{\mathrm{SM}}
=
SU(3)_C
\times
SU(2)_L
\times
U(1)_Y
}
\tag{67}
$$

を局所ゲージ群とする量子場理論であり、三世代のカイラル・フェルミオン、Higgs の複素 $SU(2)$ doublet、三つのゲージ結合、および Yukawa 結合を含む [6,7]。

従って、本論文の連続極限が標準模型全体と同一であるためには、少なくとも以下が必要である。

1. 局所 $U(1)$ ゲージ共変性
2. 非可換 $SU(2)_L$ および $SU(3)_C$ ゲージ場
3. ゲージ場自身の Yang–Mills 動力学
4. 左右非対称なカイラル・フェルミオン表現
5. Higgs doublet の動径自由度と対称性の自発的破れ
6. Yukawa 結合とフェルミオン質量
7. 量子化、経路積分、演算子構造、またはそれと等価な量子場理論の構成
8. anomaly cancellation を含む表現整合性

本論文で厳密に導出したのは、この全てではない。

本論文が直接導出したのは、

$$
\boxed{
\text{固定振幅複素位相場の }
U(1)
\text{ current とその連続保存則}
}
\tag{68}
$$

である。

従って、「標準模型に接続した」という表現は、以下の限定された意味で用いる。

$$
\boxed{
\text{標準模型に含まれる複素スカラー／位相 }
U(1)
\text{ sector と同型の連続場構造へ接続した}
}
\tag{69}
$$

標準模型全体を導出した、という意味ではない。

### 18.2 本論文の連続作用と標準的複素スカラー場

本論文の連続作用は、

$$
S_{\mathrm{cont}}
=
\frac{A^2}{2}
\int
g^{\mu\nu}
\partial_\mu\phi
\partial_\nu\phi
\,d^dq.
\tag{70}
$$

固定振幅複素場

$$
\Psi
=
Ae^{i\phi}
\tag{71}
$$

を考えると、

$$
\partial_\mu\Psi
=
iAe^{i\phi}\partial_\mu\phi,
$$

したがって、

$$
(\partial_\mu\Psi)^*
(\partial^\mu\Psi)
=
A^2
\partial_\mu\phi
\partial^\mu\phi.
\tag{72}
$$

従って、本論文の連続作用は固定振幅複素スカラー場の kinetic term と一致する。

Noether current も、

$$
j^\mu
=
\frac{1}{2i}
\left(
\Psi^*\partial^\mu\Psi
-
\Psi\partial^\mu\Psi^*
\right)
$$

より、

$$
\boxed{
j^\mu
=
A^2\partial^\mu\phi
}
\tag{73}
$$

となり、本論文の式と一致する。

従って、

$$
\boxed{
\partial_\mu j^\mu=0
}
\tag{74}
$$

は標準的な global $U(1)$ 位相対称性の保存 current と矛盾しない。

### 18.3 固定振幅という条件の意味

標準模型の Higgs 場は固定振幅の単一複素 scalar ではなく、複素 $SU(2)_L$ doublet であり、動径自由度を持つ [6,7]。

本論文は、

$$
|\Psi|=A
$$

を固定している。

従って標準模型 Higgs sector 全体ではなく、振幅モードを固定した phase-only sector に対応する。

これは矛盾ではない。

ただし、

$$
\boxed{
\text{Higgs boson の動径励起そのものは本論文からまだ導出していない}
}
\tag{75}
$$

ことを意味する。

### 18.4 global $U(1)$ から local $U(1)$ への拡張

標準模型の $U(1)_Y$ は global symmetry ではなく local gauge symmetry である。

従って、

$$
\phi(q)
\mapsto
\phi(q)+\alpha(q)
\tag{76}
$$

という位置依存変換に対して不変である必要がある。

本論文の裸の位相差、

$$
\phi_j-\phi_i
$$

は local 変換では、

$$
\phi_j-\phi_i
\mapsto
\phi_j-\phi_i
+
\alpha_j-\alpha_i
$$

となるため、そのままでは local gauge invariant ではない。

この点は標準模型との最初の明確な差である。

しかし離散系では、辺に gauge link

$$
\boxed{
U_{ij}^{(g)}
=
e^{-ig a A_{ij}}
}
\tag{77}
$$

を置けば解消できる。

局所変換を、

$$
X_i\mapsto e^{ig\alpha_i}X_i
$$

とし、

$$
U_{ij}^{(g)}
\mapsto
e^{ig\alpha_i}
U_{ij}^{(g)}
e^{-ig\alpha_j}
$$

とすれば、

$$
X_i^*
U_{ij}^{(g)}
X_j
$$

は gauge invariant である。

従って gauge-covariant current は、

$$
\boxed{
J_{ij}^{(g)}
=
\operatorname{Im}
\left(
X_i^*U_{ij}^{(g)}X_j
\right)
}
\tag{78}
$$

となる。

$U(1)$ では、

$$
J_{ij}^{(g)}
=
A^2
\sin
\left[
\phi_j-\phi_i
-gaA_{ij}
\right].
\tag{79}
$$

高分解能極限

$$
a\to0
$$

では、

$$
\phi_j-\phi_i
=
a\partial_\mu\phi+O(a^2),
$$

$$
aA_{ij}
=
aA_\mu+O(a^2)
$$

なので、

$$
\frac{1}{a}
\left(
\phi_j-\phi_i
-gaA_{ij}
\right)
\longrightarrow
\partial_\mu\phi-gA_\mu.
$$

従って、

$$
\boxed{
J^\mu_{(g)}
=
A^2
\left(
\partial^\mu\phi
-
gA^\mu
\right)
}
\tag{80}
$$

を得る。

これは、

$$
D_\mu
=
\partial_\mu
-
igA_\mu
$$

による標準的 gauge-covariant derivative の位相表示と一致する。

従って、本論文の離散 current は link variable を導入することで、局所 $U(1)$ gauge current へ自然に拡張できる。

### 18.5 gauge 場自身の動力学

標準模型では gauge field は単なる補助変数ではなく、自身の kinetic term を持つ [6,7]。

$U(1)$ なら、

$$
-\frac14
F_{\mu\nu}F^{\mu\nu},
$$

非可換群なら、

$$
-\frac14
F_{\mu\nu}^aF^{a\mu\nu}
$$

が必要である。

離散 gauge theory では、link variable の plaquette 積から Wilson gauge action を構成し、link spacing を 0 にする極限で Yang–Mills action が回復することが標準的に知られている [8]。

従って標準模型の gauge dynamics まで接続するための離散拡張は、

$$
\boxed{
\text{頂点位相}
+
\text{辺 link variable}
+
\text{plaquette action}
}
\tag{81}
$$

である。

本論文は頂点位相と辺 current までを導出した。

plaquette gauge action はまだ本公理から導出していない。

### 18.6 非可換 $SU(2)$ と $SU(3)$

前論文では $SU(2)$、$SU(3)$ へ接続する対称構造を扱った。

しかし本論文の動力学は現状、

$$
\phi_i\in U(1)
$$

という可換位相に対して明示的に閉じている。

非可換化では、scalar phase を group-valued state に置き換え、

$$
G_i\in SU(n),
$$

辺 transporter を、

$$
U_{ij}\in SU(n)
$$

とし、

$$
G_i^\dagger U_{ij}G_j
$$

の Lie algebra 成分を current として読む必要がある。

連続極限では、

$$
U_{ij}
=
\exp
\left(
igaA_\mu^aT^a
\right)
$$

から、

$$
\boxed{
D_\mu
=
\partial_\mu
-
igA_\mu^aT^a
}
\tag{82}
$$

へ接続する。

この構造は標準的 lattice gauge theory と整合する [8]。

従って、

$$
\boxed{
U(1)
\to
SU(2),SU(3)
}
$$

への拡張に数学的障害は見えないが、非可換 current と plaquette dynamics の公理からの一意導出は本論文の範囲外である。

### 18.7 フェルミオンとの比較

標準模型の matter field はカイラル fermion であり、kinetic term は、

$$
\bar\psi i\gamma^\mu D_\mu\psi
$$

型である [6,7]。

本論文の連続方程式、

$$
\Box_g\phi=0
$$

は scalar / phase field の方程式であり、Dirac 方程式ではない。

従って、

$$
\boxed{
\text{本論文の PDE を fermion 方程式と同一視することはできない}
}
\tag{83}
$$

。

ただし、前論文または別論文で得られた Bose/Fermi 分類と本 current 構造を接続する余地は残る。

標準模型への完全接続には、spinor representation と discrete Dirac operator の導出が必要である。

### 18.8 質量項

本論文の phase-only action は、

$$
\mathcal L
=
\frac{A^2}{2}
\partial_\mu\phi
\partial^\mu\phi
$$

であり、

$$
m^2\phi^2
$$

型の mass term を持たない。

従って連続極限は、

$$
\Box\phi=0
$$

という massless scalar equation である。

標準模型では Higgs potential と Yukawa coupling によって質量構造が生じる [6,7]。

従って、

$$
\boxed{
\text{質量生成は本論文の連続極限からはまだ出ていない}
}
\tag{84}
$$

。

これは矛盾ではなく、phase-only sector の範囲を示す。

### 18.9 量子場理論としての検証

本論文で導出したのは古典的な離散作用、current、自己写像、およびその連続極限である。

標準模型は量子場理論である。

従って完全な接続には、

$$
\boxed{
\text{離散状態空間上の量子化}
}
$$

または、

$$
\boxed{
Z
=
\int
\mathcal D\phi\,
e^{iS}
}
$$

に対応する測度・振幅構造を導出する必要がある。

本論文の Noether 保存則との一致は、古典場方程式レベルでは成立する。

量子補正、renormalization、anomaly までの一致はまだ検証していない。

### 18.10 標準模型との矛盾検査

以上から、本論文の結果と標準模型との関係は次のように整理できる。

| 項目 | 本論文 | 標準模型 | 判定 |
|---|---|---|---|
| 複素位相場 | あり | あり | 整合 |
| global $U(1)$ current | 導出 | 標準構造 | 整合 |
| $\partial_\mu J^\mu=0$ | 導出 | Noether 保存則 | 整合 |
| Lorentz 型 PDE | 読出し計量で $\Box\phi=0$ | 相対論的 QFT | 整合 |
| 固定振幅 | 仮定・既導出 | Higgs 等では一般に可変 | phase-only sector で整合 |
| local $U(1)$ | link 変数で拡張可能 | 必須 | 未導出だが接続式を構成 |
| $SU(2)_L$ | 対称構造は前論文 | 局所 gauge dynamics 必須 | 動力学未完 |
| $SU(3)_C$ | 対称構造は前論文 | Yang–Mills dynamics 必須 | 動力学未完 |
| gauge kinetic term | 未導出 | 必須 | 未完 |
| fermion / chirality | 本論文では未導出 | 必須 | 未完 |
| Higgs radial mode | なし | 必須 | 未完 |
| Yukawa / masses | なし | 必須 | 未完 |
| anomaly cancellation | 未検証 | 必須 | 未完 |
| 量子化 | 未完 | 必須 | 未完 |

従って、現時点で標準模型との直接矛盾は見つからない。

一方で、現時点で確立した接続を正確に書けば、

$$
\boxed{
\text{本論文は標準模型全体ではなく、
その複素位相／global }U(1)
\text{ 保存構造と固定振幅 scalar sector に接続した}
}
\tag{85}
$$

である。

### 18.11 標準模型へ接続する最短の次段階

本論文の構造を標準模型へ拡張する最短経路は、

$$
\boxed{
\phi_i
\to
(\phi_i,U_{ij})
}
$$

である。

すなわち、

1. 頂点に matter phase
2. 辺に gauge transporter
3. 面に plaquette curvature

を置く。

これにより、

$$
\boxed{
\text{頂点}
\to
\text{matter field}
}
$$

$$
\boxed{
\text{辺}
\to
\text{connection}
}
$$

$$
\boxed{
\text{面}
\to
\text{curvature / field strength}
}
$$

という simplex 階層がそのまま標準 gauge theory の幾何へ接続する。

この点は本公理系の simplex 閉包と特に整合的である。

したがって次論文の明確な課題は、

$$
\boxed{
\text{辺 link variable と plaquette curvature を
零閉包・有限回帰・自己無撞着から導出し、
}
SU(3)_C\times SU(2)_L\times U(1)_Y
\text{ の局所 gauge dynamics を閉じる}
}
\tag{86}
$$

ことである。

## 20. 離散動力学から連続動力学へ

離散自己写像の微小反復極限を取る。

反復間隔を $\Delta s$ とし、

$$
\eta=\Delta s.
$$

式 (32) の接方向部分は、

$$
\frac{
\phi^{(n+1)}-\phi^{(n)}
}{
\Delta s
}
=
P_\phi F
+
O(\Delta s).
$$

従って、

$$
\boxed{
\Delta s\to0
\quad\Longrightarrow\quad
\frac{d\phi}{ds}
=
P_\phi F.
}
\tag{67}
$$

さらに空間読出しも $N\to\infty$ とすると、

$$
F_i
=
\sum_{j\sim i}
\sin(\phi_j-\phi_i)
$$

は二階差分へ移る。

規則格子で、

$$
F_N(q)
=
\sum_\mu
\left[
\sin(\phi(q+h_Ne_\mu)-\phi(q))
+
\sin(\phi(q-h_Ne_\mu)-\phi(q))
\right].
\tag{68}
$$

Taylor 展開により、

$$
\boxed{
F_N(q)
=
h_N^2
\sum_\mu
\partial_\mu^2\phi
+
O(h_N^4).
}
\tag{69}
$$

従って、反復パラメータを

$$
\tau
=
s\,h_N^2
$$

ではなく、非自明極限を保つため

$$
d\tau
=
h_N^2\,ds
$$

として再規格化すると、

$$
\boxed{
\frac{\partial\phi}{\partial\tau}
=
P_\phi
\left(
\sum_\mu
\partial_\mu^2\phi
\right).
}
\tag{70}
$$

制約反力を Lagrange multiplier で書けば、

$$
\boxed{
\frac{\partial\phi_i}{\partial\tau}
=
\Box_g\phi_i
-
\lambda_R\frac{\partial C_R}{\partial\phi_i}
-
\lambda_I\frac{\partial C_I}{\partial\phi_i}.
}
\tag{71}
$$

これは自己無撞着解へ向かう連続 relaxational dynamics である。

一方、物理的な Lorentz 読出しを選び、作用の Euler–Lagrange 方程式として読む場合は、式 (65) の hyperbolic PDE が物理場方程式になる。

従って、

$$
\boxed{
\text{反復パラメータ }s
\neq
\text{物理時間 }t
}
\tag{72}
$$

である。

この区別により、時間を基礎から特権化せずに、離散自己写像と標準的連続場方程式を両立できる。

---

## 21. 導出鎖の完全形

本論文の全導出は、

$$
\boxed{
\sum_iX_i^2=0
}
$$

$$
\Downarrow
$$

$$
\boxed{
X_i=Ae^{i\phi_i},
\qquad
\sum_i e^{2i\phi_i}=0
}
$$

$$
\Downarrow
$$

$$
\boxed{
J_{ij}
=
A^2\sin(\phi_j-\phi_i)
}
$$

$$
\boxed{
S_N
=
-A^2
\sum_{\langle ij\rangle}
\cos(\phi_j-\phi_i)
}
$$

$$
\Downarrow
$$

$$
\boxed{
F
=
-\nabla_\phi S_N/A^2
}
$$

$$
\Downarrow
$$

零閉包接空間へ射影：

$$
\boxed{
F_T=P_\phi F
}
$$

$$
\Downarrow
$$

有限 retraction：

$$
\boxed{
\phi^{(n+1)}
=
R_{\phi^{(n)}}(\eta F_T)
}
$$

$$
\Downarrow
$$

$$
\boxed{
\mathcal F_N:
\mathcal Z_N
\to
\mathcal Z_N
}
$$

$$
\Downarrow
$$

自己無撞着固定点：

$$
\boxed{
P_{\phi_*}F(\phi_*)=0
}
$$

$$
\Downarrow
$$

制約なし current sector：

$$
\boxed{
\operatorname{div}_dJ=0
}
$$

$$
\Downarrow
\quad
N\to\infty
$$

$$
\boxed{
J^\mu
=
A^2g^{\mu\nu}\partial_\nu\phi
}
$$

$$
\boxed{
\partial_\mu J^\mu=0
}
$$

$$
\boxed{
\Box_g\phi=0
}
$$

である。

---

## 22. 結論

前論文では対称性を導出したが、Noether 保存則と動力学への接続が残った。

本論文では、等振幅複素状態の関係位相から離散 current と離散作用を構成し、さらに零閉包を厳密に保存する接空間射影と retraction を用いて、

$$
\boxed{
\mathcal F_N:
\mathcal Z_N
\to
\mathcal Z_N
}
$$

という離散自己写像を構成した。

従って「次の位相」が決まるだけではなく、その次状態が同じ零閉包許容集合に残ることまで導出された。

さらに、

$$
N\to\infty
$$

で、

$$
J^\mu=A^2g^{\mu\nu}\partial_\nu\phi
$$

および、

$$
\partial_\mu J^\mu=0
$$

を得る。

離散作用の連続極限は、

$$
S_{\mathrm{cont}}
=
\frac{A^2}{2}
\int
g^{\mu\nu}
\partial_\mu\phi
\partial_\nu\phi
\,d^dq
$$

であり、その Euler–Lagrange 方程式は、

$$
\boxed{
\partial_\mu
\left(
A^2g^{\mu\nu}\partial_\nu\phi
\right)=0.
}
$$

一定振幅・一定計量では、

$$
\boxed{
\Box_g\phi=0
}
$$

となる。


さらに、局所位相原点の無名性を要求すると、辺 connection、面 curvature、covariant derivative が必要になり、連続極限では Maxwell / Yang–Mills 型 gauge dynamics へ接続する。

従って、標準模型との接続は fixed-amplitude $U(1)$ phase sector にとどまらず、

$$
\boxed{
\text{matter phase}
\to
\text{connection}
\to
\text{curvature}
\to
\text{gauge-covariant dynamics}
}
$$

まで具体化された。

従って本論文で得られた接続は、

$$
\boxed{
\text{零閉包}
\to
\text{離散関係位相}
\to
\text{零閉包保存自己写像}
\to
\text{離散 Noether 保存則}
\to
N\to\infty
\to
\text{標準的偏微分場方程式}
}
$$

である。

---

## 参考文献

1. 木原範昭，「物理学の対称性は、本当に最初から与える必要があるのか」，前論文.
2. E. Noether, “Invariante Variationsprobleme,” *Nachrichten von der Gesellschaft der Wissenschaften zu Göttingen, Mathematisch-Physikalische Klasse*, 1918, 235–257.
3. J. E. Marsden and M. West, “Discrete Mechanics and Variational Integrators,” *Acta Numerica*, 10 (2001), 357–514. DOI: 10.1017/S096249290100006X.
4. V. A. Dorodnitsyn, “Noether-type theorems for difference equations,” *Applied Numerical Mathematics*, 39 (2001), 307–321. DOI: 10.1016/S0168-9274(00)00041-6.
5. M. Skopenkov, “Discrete Field Theory: Symmetries and Conservation Laws,” *Mathematical Physics, Analysis and Geometry*, 26 (2023), Article 19. DOI: 10.1007/s11040-023-09459-4.
6. F. Takahashi et al. (Particle Data Group), *Review of Particle Physics*, Int. J. Mod. Phys. A 41, 2630011 (2026), Standard Model reviews.
7. CERN, *The Standard Model*, CERN-TH-2024-106, Chapter 8, 2024.
8. Standard lattice gauge theory review: link variables $U_\mu(n)$, Wilson plaquette action, and continuum Yang–Mills limit; see e.g. CERN/INSPIRE lattice gauge theory lecture material cited in the text.
