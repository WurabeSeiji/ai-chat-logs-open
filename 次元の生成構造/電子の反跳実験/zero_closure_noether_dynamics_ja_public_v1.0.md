# 離散零閉包からの Noether 保存則と関係位相動力学――零閉包を厳密保存する離散自己写像、$N\to\infty$ 連続場方程式、局所ゲージ幾何、および標準模型一世代表現と chirality 選択

**著者:** 木原 範昭<br>
**ORCID:** 0009-0004-6753-4020<br>
**日付:** 2026年8月21日<br>
**Version DOI:** 10.5281/zenodo.22040736<br>
**Concept DOI:** 10.5281/zenodo.22040735<br>
**位置づけ:** 「次元の生成構造」シリーズ・閉包公理からの対称性導出 続編（動力学・保存則・標準模型表現）公開版 v1.0<br>
**前論文:** 零閉包・有限位数・自己無撞着幾何からの対称性生成 v1.0（Concept DOI 10.5281/zenodo.22028072）<br>
**ライセンス:** CC BY 4.0

> **主題**：前論文で残った二つの課題――Noether 型保存則と、次状態を決める動力学への接続――を扱う。動力学に関する中心主張は、状態の書換えを暗黙に許さず、動力学を許容状態空間の自己写像として構成することである。零閉包を各有限反復で厳密保存する離散自己写像を構成し、その $N\to\infty$ 極限で標準的な連続場方程式を得る。さらに局所ゲージ幾何、標準模型一世代表現、hypercharge、anomaly cancellation、chirality 選択までの理論接続を与え、数値検証仕様を付す。

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

ここで本論文の動力学に関する中心的な主張を明記する。**状態を書き換えること自体を動力学として暗黙に許さない。** 本公理系では自己無撞着が基礎条件であるため、次状態への写像は、許容状態空間を破壊しない自己写像でなければならない。とくに本論文で明示的に扱う零閉包については、位相 $\phi$ のみを更新しながら

$$
\boxed{
\mathcal F_N:\mathcal Z_N\to\mathcal Z_N
}
$$

を構成し、各有限反復で

$$
\sum_iX_i^2=0
$$

が厳密に保存されることを示す。したがって、**離散系の動力学は、公理を一度破って後から拘束し直すのではなく、零閉包を保つ許容位相写像として実現できる。** これは本論文における動力学導出の前提ではなく、検証すべき結果である。有限回帰、simplex 閉包、自己無撞着を含む他の既導出条件についても、後続の数値検証では「更新後にも保存されること」を独立の監査項目とし、保存されない更新則は本公理系の動力学として棄却する。

したがって本論文では、離散関係位相系から、零閉包を厳密保存する自己写像、Noether 型局所保存則、および標準的な偏微分場方程式への連続極限を一つの導出鎖として与える。

さらに前論文の5複素自由度と $S(U(3)\times U(2))$ を本論文の局所 gauge dynamics に接続すると、$U(1)$ generator の trace-zero 条件から hypercharge 比が固定され、5自由度と simplex 二体関係から

$$
V^*\oplus\Lambda^2V
$$

が現れる。この表現は標準模型一世代の

$$
d^c,\ L,\ u^c,\ Q,\ e^c
$$

の15左手 Weyl 成分に分解され、全 perturbative gauge anomaly および $SU(2)$ global anomaly の相殺を直接確認できる。標準模型の内部表現の導出鎖で最後に残っていた共役二 Weyl sector の選択についても、過去論文の A/B 二状態選択系を Weyl 共役 sector と同定し、mirror-odd 内部相関
$$
J=\operatorname{Im}(B^{*2}C)
$$
と既存の非線形選択項から
$$
\dot S_\chi=\lambda J+gS_\chi(1-S_\chi^2)
$$
という最小 normal form を得る。これにより、本論文が対象とする gauge group・局所 gauge geometry・一世代内部表現・hypercharge・anomaly cancellation・chirality selection までの理論接続は閉じる。残るのは、この同定と選択則が既存の数値 dynamics 上で成立するかの数値検証である。

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

ここで第二の問題には、本公理系に固有の厳しい条件がある。通常の場の理論では、状態 $\Phi$ から次状態 $\Phi'$ への更新則をまず置き、その更新を「時間発展」として扱うことができる。しかし本公理系では、**自己無撞着そのものが基礎条件である。したがって状態の書換えを無条件に仮定することはできない。** 更新後の状態が零閉包、有限回帰、simplex 閉包その他の許容条件から外れるなら、その更新則は本理論の動力学ではない。

従って本論文で問うべき動力学の問題は、単に

$$
\phi^{(n)}\mapsto\phi^{(n+1)}
$$

を与えることではなく、

$$
\boxed{
\phi^{(n)}\in\mathcal Z_N
\Longrightarrow
\phi^{(n+1)}\in\mathcal Z_N
}
$$

を満たす**許容状態空間内部の自己写像**を構成できるか、である。

検討の結果、Noether 型保存則は有限 $N$ の離散 current 保存則の $N\to\infty$ 連続極限として現れた。さらに同じ関係位相 current が、次位相を決める自己写像の生成量になることが分かった。しかも、零閉包多様体の接空間射影と有限 retraction を用いることで、位相 $\phi$ を更新しながら零閉包を有限ステップごとに厳密保存できる。**「次状態を作れる」ことと「公理を壊さず次状態を作れる」ことは別問題であり、本論文が閉じるのは後者である。**

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

この式の重要性は、単に数値計算を安定化する拘束処理にあるのではない。本公理系では自己無撞着が出発条件であるため、状態更新は許容集合の外へ出ることを許されない。式 (32)--(34) は、**位相を実際に更新しても零閉包を破壊しない離散自己写像が構成可能であること**を示している。すなわち、動力学を導入するために零閉包を一時的に捨てる必要はない。

$$
\boxed{
\text{許容状態}
\xrightarrow{\;\phi\text{ の更新}\;}
\text{許容状態}
}
$$

が有限 $N$ で成立する。これは自己無撞着を基礎公理に置く本系では必須の検証であり、暗黙の状態書換えでは代替できない。

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


## 19. 前論文の対称性導出との接続と標準模型表現の閉包

本節では、前論文で得られた対称性構造を、本論文で導出した局所 gauge connection・curvature・covariant derivative・Yang--Mills 型動力学へ接続する。

前論文では、複素6軸零閉包

$$
\sum_{n=1}^{6}X_n^2=0
$$

から複素拘束1本を差し引き、

$$
\boxed{
\dim_{\mathbb C}=5
}
\tag{95}
$$

という独立複素自由度を得た。

その5自由度を自己無撞着な Hermitian 分解

$$
\boxed{
V=V_3\oplus V_2,
\qquad
\dim_{\mathbb C}V_3=3,
\qquad
\dim_{\mathbb C}V_2=2
}
\tag{96}
$$

として読むと、分解保存群は

$$
U(3)\times U(2)
$$

である。

さらに全体位相冗長性を除くと、

$$
\boxed{
S(U(3)\times U(2))
\cong
\frac{
SU(3)\times SU(2)\times U(1)
}{
\mathbb Z_6
}
}
\tag{97}
$$

を得る。

前論文ではここまでを標準模型 global gauge group への条件付き厳密な接続としていた。

本論文では、局所位相原点の無名性から辺 connection が必要となり、simplex 面閉路から curvature が得られ、連続極限で

$$
D_\mu
=
\partial_\mu-igA_\mu
$$

および

$$
F_{\mu\nu}
=
\partial_\mu A_\nu-\partial_\nu A_\mu-ig[A_\mu,A_\nu]
$$

が現れることを前節で導出した。

従って、前論文の global stabilizer と本論文の local transporter を接続すると、

$$
\boxed{
S(U(3)\times U(2))
\quad\Longrightarrow\quad
\text{局所 }
\frac{SU(3)\times SU(2)\times U(1)}{\mathbb Z_6}
\text{ connection}
}
\tag{98}
$$

という局所 gauge dynamics への橋が得られる。

以下では、この5複素自由度構造から標準模型一世代の内部量子数がどこまで固定されるかを順に検算する。

### 19.1 $U(1)$ generator は trace-zero 条件で比が固定される

$V=V_3\oplus V_2$ 上で、$SU(3)$ と $SU(2)$ に可換な $U(1)$ generator は block diagonal に、

$$
Y
=
\begin{pmatrix}
y_3 I_3 & 0\\
0 & y_2 I_2
\end{pmatrix}
\tag{99}
$$

と書ける。

ここで $S(U(3)\times U(2))$ の Lie algebra 条件は全 trace が 0 であることなので、

$$
\operatorname{Tr}Y
=
3y_3+2y_2
=
0.
\tag{100}
$$

従って、

$$
\boxed{
3y_3+2y_2=0
}
\tag{101}
$$

であり、

$$
\boxed{
y_3:y_2=-2:3
}
\tag{102}
$$

と比が一意に決まる。

全体規格化は $U(1)$ charge の単位選択に対応する。標準模型の慣習に合わせて、

$$
\boxed{
y_3=-\frac13,
\qquad
y_2=\frac12
}
\tag{103}
$$

と取れば、

$$
\boxed{
Y
=
\operatorname{diag}
\left(
-\frac13,-\frac13,-\frac13,
\frac12,\frac12
\right)
}
\tag{104}
$$

である。

ここで重要なのは、$-1/3$ と $1/2$ を独立に入力していないことである。比 $-2:3$ は $3+2$ 分解と trace-zero 条件だけで固定され、残る共通倍数は charge の単位規格化である。

従って、

$$
\boxed{
\text{hypercharge 比}
\quad
3\left(-\frac13\right)
+
2\left(\frac12\right)
=0
}
\tag{105}
$$

は $S(U(3)\times U(2))$ の構造から得られる。

### 19.2 5自由度と simplex 二体関係

前論文の A3 は、$N$ 個の頂点と全二体関係を simplex として扱う。

5複素自由度を一頂点レジスタ空間

$$
V\simeq\mathbb C^5
$$

として読むと、向き付き二体関係は反対称二階テンソル

$$
\boxed{
\Lambda^2V
}
\tag{106}
$$

で表される。

その複素次元は、

$$
\boxed{
\dim_{\mathbb C}\Lambda^2V
=
\binom52
=
10
}
\tag{107}
$$

である。

一方、一体側の双対表現は

$$
V^*
$$

であり、

$$
\dim_{\mathbb C}V^*=5.
$$

従って、

$$
\boxed{
V^*\oplus\Lambda^2V
}
\tag{108}
$$

は、

$$
\boxed{
5+10=15
}
\tag{109}
$$

複素成分を持つ。

これは単なる次元一致として終わらない。次節で $S(U(3)\times U(2))$ 表現へ分解すると、標準模型一世代の左手 Weyl 表現と一致する。

### 19.3 $V^*$ の分解

式 (96)、(103) より、

$$
V
=
(\mathbf3,\mathbf1)_{-1/3}
\oplus
(\mathbf1,\mathbf2)_{1/2}.
\tag{110}
$$

双対表現は、

$$
\boxed{
V^*
=
(\overline{\mathbf3},\mathbf1)_{1/3}
\oplus
(\mathbf1,\mathbf2)_{-1/2}
}
\tag{111}
$$

である。

これは標準模型の左手 Weyl 記法で、

$$
\boxed{
d^c
:
(\overline{\mathbf3},\mathbf1)_{1/3}
}
\tag{112}
$$

および

$$
\boxed{
L
:
(\mathbf1,\mathbf2)_{-1/2}
}
\tag{113}
$$

の内部量子数と一致する。

### 19.4 $\Lambda^2V$ の分解

直和の exterior square は、

$$
\Lambda^2(V_3\oplus V_2)
=
\Lambda^2V_3
\oplus
(V_3\otimes V_2)
\oplus
\Lambda^2V_2.
\tag{114}
$$

これを各項ごとに計算する。

#### 19.4.1 $\Lambda^2V_3$

$SU(3)$ の基本表現について、

$$
\Lambda^2\mathbf3
\cong
\overline{\mathbf3}.
$$

hypercharge は加法的なので、

$$
Y
=
-\frac13-\frac13
=
-\frac23.
$$

従って、

$$
\boxed{
\Lambda^2V_3
=
(\overline{\mathbf3},\mathbf1)_{-2/3}
}
\tag{115}
$$

である。

これは、

$$
\boxed{
u^c
:
(\overline{\mathbf3},\mathbf1)_{-2/3}
}
\tag{116}
$$

に一致する。

#### 19.4.2 $V_3\otimes V_2$

表現は、

$$
(\mathbf3,\mathbf2)
$$

であり、hypercharge は、

$$
-\frac13+\frac12
=
\frac16.
$$

従って、

$$
\boxed{
V_3\otimes V_2
=
(\mathbf3,\mathbf2)_{1/6}
}
\tag{117}
$$

である。

これは、

$$
\boxed{
Q
:
(\mathbf3,\mathbf2)_{1/6}
}
\tag{118}
$$

に一致する。

#### 19.4.3 $\Lambda^2V_2$

$SU(2)$ の基本二重項では、

$$
\Lambda^2\mathbf2
\cong
\mathbf1.
$$

hypercharge は、

$$
\frac12+\frac12=1.
$$

従って、

$$
\boxed{
\Lambda^2V_2
=
(\mathbf1,\mathbf1)_1
}
\tag{119}
$$

である。

これは、

$$
\boxed{
e^c
:
(\mathbf1,\mathbf1)_1
}
\tag{120}
$$

に一致する。

### 19.5 一世代表現

式 (111)、(115)、(117)、(119) をまとめると、

$$
\boxed{
V^*\oplus\Lambda^2V
=
(\overline3,1)_{1/3}
\oplus
(1,2)_{-1/2}
\oplus
(\overline3,1)_{-2/3}
\oplus
(3,2)_{1/6}
\oplus
(1,1)_1
}
\tag{121}
$$

である。

すなわち、

$$
\boxed{
V^*\oplus\Lambda^2V
=
d^c\oplus L\oplus u^c\oplus Q\oplus e^c
}
\tag{122}
$$

となる。

これは右手粒子を左手 charge-conjugate field で表す標準的な左手 Weyl 記法での、ニュートリノ右手成分を含まない標準模型一世代15成分と一致する。

群論的には、

$$
\boxed{
V^*\oplus\Lambda^2V
\cong
\overline{\mathbf5}
\oplus
\mathbf{10}
}
\tag{123}
$$

という既知の分解と同型である。

ただし本論文では $SU(5)$ 大統一群を出発点に置いていない。

導出順序は、

$$
\boxed{
\text{複素零閉包}
\to
5\text{独立複素自由度}
\to
3\oplus2
\to
S(U(3)\times U(2))
\to
V^*\oplus\Lambda^2V
}
\tag{124}
$$

である。

### 19.6 anomaly cancellation の直接検算

式 (121) の表現内容について、perturbative gauge anomaly を直接計算する。

以下ではすべて左手 Weyl field で数える。

#### 19.6.1 $SU(3)^3$

$Q=(3,2)$ は $SU(2)$ の二成分を持つため、$SU(3)$ fundamental が2個ある。

$u^c,d^c$ は anti-fundamental が各1個である。

fundamental と anti-fundamental の cubic anomaly coefficient の符号をそれぞれ $+1,-1$ とすれば、

$$
2-1-1=0.
$$

従って、

$$
\boxed{
\mathcal A_{SU(3)^3}=0
}
\tag{125}
$$

である。

#### 19.6.2 $SU(3)^2U(1)_Y$

$T(\mathbf3)=T(\overline{\mathbf3})=1/2$ を用いる。

$Q$ は $SU(2)$ 二重項なので multiplicity 2 を持つ。

$$
\mathcal A_{SU(3)^2U(1)}
=
2\left(\frac16\right)\frac12
+
\left(-\frac23\right)\frac12
+
\left(\frac13\right)\frac12.
$$

従って、

$$
\mathcal A_{SU(3)^2U(1)}
=
\frac16-\frac13+\frac16
=
0.
$$

すなわち、

$$
\boxed{
\mathcal A_{SU(3)^2U(1)}=0
}
\tag{126}
$$

である。

#### 19.6.3 $SU(2)^2U(1)_Y$

$T(\mathbf2)=1/2$ を用いる。

$Q$ は color multiplicity 3 を持つので、

$$
\mathcal A_{SU(2)^2U(1)}
=
3\left(\frac16\right)\frac12
+
\left(-\frac12\right)\frac12.
$$

従って、

$$
\frac14-\frac14=0.
$$

よって、

$$
\boxed{
\mathcal A_{SU(2)^2U(1)}=0
}
\tag{127}
$$

である。

#### 19.6.4 $U(1)_Y^3$

全左手 Weyl 成分を multiplicity 付きで数えると、

$$
\mathcal A_{U(1)^3}
=
6\left(\frac16\right)^3
+
3\left(-\frac23\right)^3
+
3\left(\frac13\right)^3
+
2\left(-\frac12\right)^3
+
1^3.
\tag{128}
$$

各項は、

$$
\frac1{36},
\qquad
-\frac89,
\qquad
\frac19,
\qquad
-\frac14,
\qquad
1
$$

なので、

$$
\frac1{36}
-\frac{32}{36}
+\frac4{36}
-\frac9{36}
+\frac{36}{36}
=
0.
$$

従って、

$$
\boxed{
\mathcal A_{U(1)^3}=0
}
\tag{129}
$$

である。

#### 19.6.5 gravitational--$U(1)_Y$

hypercharge の総和は、

$$
6\left(\frac16\right)
+
3\left(-\frac23\right)
+
3\left(\frac13\right)
+
2\left(-\frac12\right)
+
1.
$$

従って、

$$
1-2+1-1+1=0.
$$

ゆえに、

$$
\boxed{
\mathcal A_{\mathrm{grav}^2U(1)}=0
}
\tag{130}
$$

である。

#### 19.6.6 $SU(2)$ global anomaly

左手 $SU(2)$ doublet の総数は、color multiplicity を含めて、

$$
3\quad(Q)
+
1\quad(L)
=
4.
$$

従って偶数であり、

$$
\boxed{
N_{\mathrm{doublet}}=4
}
\tag{131}
$$

なので $SU(2)$ global anomaly の条件も満たす。

以上から、

$$
\boxed{
V^*\oplus\Lambda^2V
\text{ の一世代表現は標準模型の anomaly cancellation を満たす}
}
\tag{132}
$$

ことが直接確認できる。

### 19.7 hypercharge と anomaly cancellation は独立の追加条件ではない

ここまでの結果をまとめる。

hypercharge 比は、

$$
\boxed{
3y_3+2y_2=0
}
$$

という $S(U(3)\times U(2))$ の trace-zero 条件から出る。

一世代表現は、

$$
\boxed{
V^*\oplus\Lambda^2V
}
$$

という5自由度と simplex 二体関係から構成される。

その表現に対する anomaly cancellation は、式 (125)--(131) の直接計算で成立する。

従って、

$$
\boxed{
\text{hypercharge}
+
\text{anomaly cancellation}
}
\tag{133}
$$

を、標準模型に合わせて独立に入力する必要はない。

導出鎖は、

$$
\boxed{
5\text{複素自由度}
\to
3\oplus2
\to
S(U(3)\times U(2))
\to
Y
\to
V^*\oplus\Lambda^2V
\to
\text{anomaly cancellation}
}
\tag{134}
$$

で閉じる。

### 19.8 Lorentz spinor と chirality

前論文では、$(x,y,z,t)$ の Lorentz 部分読出しについて、

$$
SO^+(3,1)
$$

およびその二重被覆

$$
\boxed{
Spin^+(3,1)
\cong
SL(2,\mathbb C)
}
\tag{135}
$$

への接続が得られている。

$SL(2,\mathbb C)$ の最小 Weyl spinor 表現は、

$$
\boxed{
\left(\frac12,0\right)
}
\tag{136}
$$

と、

$$
\boxed{
\left(0,\frac12\right)
}
\tag{137}
$$

の二つであり、互いに複素共役である。

内部表現側にも、

$$
\boxed{
V^*\oplus\Lambda^2V
}
\tag{138}
$$

と、その複素共役

$$
\boxed{
V\oplus\Lambda^2V^*
}
\tag{139}
$$

がある。

従って mirror を含む完全な候補は、

$$
\boxed{
\left(\frac12,0\right)
\otimes
\left(
V^*\oplus\Lambda^2V
\right)
}
\tag{140}
$$

と、

$$
\boxed{
\left(0,\frac12\right)
\otimes
\left(
V\oplus\Lambda^2V^*
\right)
}
\tag{141}
$$

の共役二 sector として書ける。

標準模型を左手 Weyl field だけで表記すると、観測される一世代は式 (140) 側の内部量子数

$$
d^c,\ L,\ u^c,\ Q,\ e^c
$$

としてまとめられる。

しかし、ここで重要な論理境界がある。

前論文の Fermi 型二重被覆・奇数倍音 sector は fermionic sector の存在を与えるが、それだけでは式 (140) と式 (141) のどちらが自己無撞着な物理 sector として選択されるかを一意に決めない。

従って、標準模型表現の閉包について最後に残る問題は、

$$
\boxed{
\text{chirality sector の自己無撞着選択}
}
\tag{142}
$$

である。

これは「chirality をゼロから発生させる」問題ではない。

既に存在する共役二 sector、

$$
\left(\frac12,0\right)
\otimes
(V^*\oplus\Lambda^2V)
$$

と、

$$
\left(0,\frac12\right)
\otimes
(V\oplus\Lambda^2V^*)
$$

のうち、なぜ一方が物理的に選択されるかを確定する問題である。

### 19.9 最後の課題を現在の公理系の言葉で書く

前論文には既に、

1. 向き付き simplex
2. $\partial^2=0$ と $d^2=0$
3. Lorentz spin double cover
4. Fermi 型の一重／二重被覆分類
5. A4 の自己無撞着固定点と stabilizer selection

が存在する。

従って、最後の課題は抽象的な「左右非対称性の起源」ではなく、

$$
\boxed{
\text{simplex orientation}
+
\text{Lorentz Weyl sector}
+
\text{Fermi 型二重被覆}
+
\text{A4 自己無撞着性}
}
\tag{143}
$$

が、

$$
\boxed{
V^*\oplus\Lambda^2V
}
\quad\text{と}\quad
\boxed{
V\oplus\Lambda^2V^*
}
\tag{144}
$$

のどちらを安定固定点 sector として選ぶかを解析する問題である。

すなわち、

$$
\boxed{
\mathcal F_N
\text{ が orientation reversal / complex conjugation に対して
どの固定点を安定化するか}
}
\tag{145}
$$

を求めればよい。

この選択が一方に固定されれば、

$$
\boxed{
\text{chirality}
\to
\text{hypercharge}
\to
\text{一世代表現}
\to
\text{anomaly cancellation}
}
$$

までの標準模型内部表現の導出鎖が閉じる。

### 19.10 現時点で標準模型と一致する範囲

本論文と前論文を合わせた現時点の到達点を整理する。

| 構造 | 導出経路 | 現在の状態 |
|---|---|---|
| Lorentz 型符号 | 複素虚軸 $it$ | 導出済み |
| $Spin^+(3,1)$ | Lorentz group の二重被覆 | 接続済み |
| 5複素自由度 | $\mathbb C^6$ の零閉包1本 | 導出済み |
| $3\oplus2$ 読出し | 5自由度の Hermitian 分解 | 既導出読出し、選択則の一般化は別課題 |
| $S(U(3)\times U(2))$ | 分解保存 + 全体位相除去 | 条件付き厳密 |
| local connection | 局所位相／基底原点の無名性 | 本論文で導出 |
| curvature | simplex 面閉路 | 本論文で導出 |
| Maxwell / Yang--Mills 型作用 | plaquette 極限 | 本論文で接続 |
| hypercharge 比 | $3y_3+2y_2=0$ | 導出 |
| 一世代15成分 | $V^*\oplus\Lambda^2V$ | 導出 |
| $Q,u^c,d^c,L,e^c$ の量子数 | $3\oplus2$ 表現分解 | 導出 |
| perturbative anomaly cancellation | 表現の直接和 | 検算成立 |
| $SU(2)$ global anomaly | doublet 数4 | 検算成立 |
| Bose/Fermi/混合 sector | 奇偶倍音・二重被覆 | 自己論文で既導出・数値確認 |
| chirality sector 選択 | orientation + Weyl + A4 | **最後の未閉包点** |

従って、標準模型の内部 gauge 表現については、

$$
\boxed{
\text{最後に残るのは chirality sector の自己無撞着選択}
}
\tag{146}
$$

である。

ただし、これは標準模型の全現象をすべて導出したという意味ではない。

Higgs の動径モード、Yukawa coupling、世代数、具体的質量階層、量子補正・renormalization などは、内部 gauge 表現の閉包とは別の動力学・量子化課題である。

本節で閉じたのは、

$$
\boxed{
\text{標準模型の gauge group}
+
\text{local gauge geometry}
+
\text{一世代内部表現}
+
\text{hypercharge}
+
\text{anomaly cancellation}
}
\tag{147}
$$

までである。


## 20. chirality sector の理論閉包

前節までで、内部 gauge 表現については、

$$
5\text{複素自由度}
\to
3\oplus2
\to
S(U(3)\times U(2))
\to
Y
\to
V^*\oplus\Lambda^2V
$$

までが得られた。残っていたのは、共役な二つの Weyl sector のうち一方がどのように自己無撞着に選択されるかである。

ここでは、過去論文で既に構成されていた A/B 二チャネル選択系を Weyl 共役二 sector と同定し、その選択則を chirality の秩序変数として書き直す。

### 20.1 共役二 Weyl sector

Lorentz 側の最小 Weyl 表現を、

$$
\left(\frac12,0\right),
\qquad
\left(0,\frac12\right)
$$

とする。

内部表現側には、

$$
V^*\oplus\Lambda^2V
$$

と、その複素共役

$$
V\oplus\Lambda^2V^*
$$

がある。

従って共役二 sector を、

$$
\boxed{
\mathcal H_L
=
\left(\frac12,0\right)
\otimes
\left(
V^*\oplus\Lambda^2V
\right)
}
\tag{148}
$$

$$
\boxed{
\mathcal H_R
=
\left(0,\frac12\right)
\otimes
\left(
V\oplus\Lambda^2V^*
\right)
}
\tag{149}
$$

と定義する。

complex conjugation / orientation reversal を $\mathcal P_\chi$ と書けば、

$$
\boxed{
\mathcal P_\chi:
\mathcal H_L
\leftrightarrow
\mathcal H_R
}
\tag{150}
$$

である。

### 20.2 過去の A/B 二チャネル系との同定

過去論文では、二つの複素振幅 $a,b$ から、

$$
p_A=|a|^2,
\qquad
p_B=|b|^2
$$

を定義し、

$$
\boxed{
S
=
\frac{p_A-p_B}{p_A+p_B}
}
\tag{151}
$$

を二状態選択の秩序変数として用いた。

$A\leftrightarrow B$ で、

$$
\boxed{
S\mapsto-S
}
\tag{152}
$$

である。

ここで、

$$
\boxed{
A\equiv\mathcal H_L,
\qquad
B\equiv\mathcal H_R
}
\tag{153}
$$

と同定する。

すると、

$$
p_L:=p_A,
\qquad
p_R:=p_B
$$

であり、

$$
\boxed{
S_\chi
:=
\frac{p_L-p_R}{p_L+p_R}
}
\tag{154}
$$

は chirality order parameter となる。

mirror 変換では、

$$
\boxed{
\mathcal P_\chi:
S_\chi\mapsto-S_\chi
}
\tag{155}
$$

である。

従って、過去の A/B 選択問題は、数学的には Weyl 共役二 sector の選択問題と同型である。

### 20.3 過去論文の内部相関量は mirror-odd である

過去論文では、交差項総和を、

$$
C
=
\sum_{m<n}x_mx_n
$$

とし、A/B 間に生じる位相差を読む量として、

$$
\boxed{
J
=
\operatorname{Im}
\left(
B^{*2}C
\right)
}
\tag{156}
$$

を定義した。

全複素量を共役すると、

$$
B^{*2}C
\mapsto
\left(B^{*2}C\right)^*.
$$

従って、

$$
\boxed{
J\mapsto-J
}
\tag{157}
$$

である。

また、

$$
\Delta\phi_{AB}
=
\arg
\left(
1+\frac{2C}{B^2}
\right)
$$

も複素共役で、

$$
\boxed{
\Delta\phi_{AB}
\mapsto
-\Delta\phi_{AB}
}
\tag{158}
$$

となる。

従って、

$$
J,
\qquad
\Delta\phi_{AB},
\qquad
S_\chi
$$

はすべて同じ mirror transformation に対して odd である。

すなわち、

$$
\boxed{
\mathcal P_\chi:
(J,\Delta\phi,S_\chi)
\mapsto
(-J,-\Delta\phi,-S_\chi)
}
\tag{159}
$$

である。

### 20.4 chirality 選択の最小 normal form

$S_\chi$ の発展を、

$$
\dot S_\chi
=
F(S_\chi,J)
$$

とする。

基礎方程式が mirror 対称であるため、

$$
\boxed{
F(-S_\chi,-J)
=
-F(S_\chi,J)
}
\tag{160}
$$

を満たさなければならない。

対称点

$$
S_\chi=0,
\qquad
J=0
$$

近傍で解析的に展開すると、mirror-even な定数項、$J^2$、$S_\chi^2$ などは許されない。

最低次数では、

$$
F(S_\chi,J)
=
\lambda J
+
aS_\chi
+
bS_\chi^3
+
O(J^3,J^2S_\chi,JS_\chi^2,S_\chi^5).
\tag{161}
$$

過去論文で既に使われていた非線形選択項は、

$$
gS(1-S^2)
=
gS-gS^3
$$

である。

従って最小 normal form は、

$$
\boxed{
\dot S_\chi
=
\lambda J
+
gS_\chi
\left(
1-S_\chi^2
\right)
}
\tag{162}
$$

となる。

これは新しい任意選択則を加えた式ではない。

- $\lambda J$ は過去論文で既に定義された mirror-odd 内部位相相関が選択変数へ写る最低次数結合、
- $gS(1-S^2)$ は過去論文で既に使われた非線形選択項、

である。

### 20.5 $J=0$ における自発的 chirality 選択

まず、

$$
J=0
$$

とする。

すると、

$$
\dot S_\chi
=
gS_\chi(1-S_\chi^2).
\tag{163}
$$

固定点は、

$$
\boxed{
S_\chi^*=0,\quad +1,\quad -1
}
\tag{164}
$$

である。

右辺を、

$$
f(S)=gS(1-S^2)
$$

とすれば、

$$
f'(S)
=
g(1-3S^2).
$$

従って、

$$
f'(0)=g,
$$

$$
f'(\pm1)=-2g.
$$

$g>0$ なら、

$$
\boxed{
S_\chi=0
\text{ は不安定}
}
\tag{165}
$$

であり、

$$
\boxed{
S_\chi=\pm1
\text{ は安定}
}
\tag{166}
$$

である。

従って基礎方程式は左右対称のまま、

$$
\boxed{
S_\chi=0
\longrightarrow
S_\chi=+1
\quad\text{または}\quad
S_\chi=-1
}
\tag{167}
$$

という自発的 chirality selection が起こる。

### 20.6 内部位相相関 $J$ が選択符号を決める

$J$ が小さいが非零であるとき、

$$
\dot S_\chi
=
\lambda J
+
gS_\chi(1-S_\chi^2)
$$

である。

対称点 $S_\chi=0$ では、

$$
\left.
\dot S_\chi
\right|_{S_\chi=0}
=
\lambda J.
\tag{168}
$$

従って、

$$
\lambda J>0
$$

なら $S_\chi$ は正方向へ、

$$
\lambda J<0
$$

なら負方向へ動き始める。

その後、$g>0$ の非線形項が差を増幅して、

$$
\boxed{
\lambda J>0
\Longrightarrow
S_\chi\to+1
}
\tag{169}
$$

$$
\boxed{
\lambda J<0
\Longrightarrow
S_\chi\to-1
}
\tag{170}
$$

となる。

mirror 変換では、

$$
J\to-J,
\qquad
S_\chi\to-S_\chi
$$

なので、方程式全体は不変である。

従って、一方の chirality を基礎公理として優遇していない。

選択される符号は、自己無撞着な内部位相相関 $J$ の符号によって決まる。

### 20.7 過去論文の因果鎖との一致

過去論文では、整数倍音・位相閉包によって交差項相関が長時間再閉鎖し、

$$
J\neq0
$$

が生じ、

$$
\Delta\phi_{AB}\neq0
$$

が現れ、

$$
S\neq0
$$

へ写り、非線形選択によって一状態へ増幅されるという因果鎖を検討していた。

chirality の言葉に書き直すと、

$$
\boxed{
\text{Fermi 型整数倍音・位相閉包}
}
$$

$$
\Downarrow
$$

$$
\boxed{
\text{交差相関 }C
}
$$

$$
\Downarrow
$$

$$
\boxed{
J
=
\operatorname{Im}(B^{*2}C)
\neq0
}
$$

$$
\Downarrow
$$

$$
\boxed{
\Delta\phi_{LR}\neq0
}
$$

$$
\Downarrow
$$

$$
\boxed{
S_\chi\neq0
}
$$

$$
\Downarrow
$$

$$
\boxed{
\dot S_\chi
=
\lambda J
+
gS_\chi(1-S_\chi^2)
}
$$

$$
\Downarrow
$$

$$
\boxed{
S_\chi\to\pm1
}
\tag{171}
$$

となる。

従って、chirality selection の力学は新規に追加する必要がなく、過去の Fermi 型非線形選択機構を Weyl 共役二 sector に同定することで得られる。

### 20.8 標準模型内部表現との理論接続の完了

前節までに、

$$
S(U(3)\times U(2))
$$

から hypercharge 比が固定され、

$$
V^*\oplus\Lambda^2V
$$

から一世代15 Weyl 成分が得られ、anomaly cancellation が成立した。

本節でさらに、

$$
\mathcal H_L
\leftrightarrow
\mathcal H_R
$$

という mirror pair に対し、過去の A/B 非線形選択機構を接続した。

従って、本論文で扱う**標準模型の内部 gauge 表現と chirality selection**について、理論導出鎖は、

$$
\boxed{
\text{零閉包}
}
$$

$$
\Downarrow
$$

$$
\boxed{
5\text{複素自由度}
}
$$

$$
\Downarrow
$$

$$
\boxed{
3\oplus2
}
$$

$$
\Downarrow
$$

$$
\boxed{
S(U(3)\times U(2))
}
$$

$$
\Downarrow
$$

$$
\boxed{
SU(3)\times SU(2)\times U(1)
\text{ の局所 gauge geometry}
}
$$

$$
\Downarrow
$$

$$
\boxed{
Y
}
$$

$$
\Downarrow
$$

$$
\boxed{
V^*\oplus\Lambda^2V
=
d^c\oplus L\oplus u^c\oplus Q\oplus e^c
}
$$

$$
\Downarrow
$$

$$
\boxed{
\text{anomaly cancellation}
}
$$

$$
\Downarrow
$$

$$
\boxed{
\text{mirror-odd }J
}
$$

$$
\Downarrow
$$

$$
\boxed{
\text{chirality selection }S_\chi\to\pm1
}
\tag{172}
$$

として閉じる。

この意味で、**標準模型の gauge group、局所 gauge geometry、一世代内部表現、hypercharge、anomaly cancellation、および chirality selection までの理論的接続は完了した**。

ここでの「完了」は、本論文が対象とする gauge-representation/chirality chain に対するものである。

Higgs 動径モード、Yukawa coupling、三世代の起源、具体的質量階層、量子補正・renormalization は別の問題であり、本論文の理論閉包範囲には含めない。

---

## 21. 数値検証仕様――既存実験系を最小変更して chirality selection を検証する

理論導出が閉じた後に残るのは、既存の Fermi 型 A/B 実験系が、Weyl 共役二 sector の同定の下で理論式 (162) を再現するかの数値検証である。

新しい基礎 dynamics は追加しない。

既存の更新則・散乱則を変更せず、主として**読出し量の再定義と追加計測**を行う。

### 21.1 既存 A/B 変数の再解釈

既存コードの、

$$
A,\qquad B
$$

を、

$$
\boxed{
A\equiv\mathcal H_L,
\qquad
B\equiv\mathcal H_R
}
\tag{173}
$$

として読む。

既存の、

$$
p_A=|a|^2,
\qquad
p_B=|b|^2
$$

を、

$$
p_L:=p_A,
\qquad
p_R:=p_B
$$

とする。

既存の選択秩序変数は、

$$
\boxed{
S_\chi
=
\frac{p_L-p_R}{p_L+p_R}
}
\tag{174}
$$

としてそのまま利用する。

### 21.2 追加すべき観測量

既存の運動方程式は変更せず、各反復で以下を記録する。

$$
C
=
\sum_{m<n}x_mx_n
\tag{175}
$$

$$
\boxed{
J_\chi
=
\operatorname{Im}
\left(
B^{*2}C
\right)
}
\tag{176}
$$

$$
\boxed{
\Delta\phi_\chi
=
\arg
\left(
1+\frac{2C}{B^2}
\right)
}
\tag{177}
$$

$$
S_\chi
$$

および、

$$
\boxed{
G_\chi
=
gS_\chi
\left(
1-S_\chi^2
\right)
}
\tag{178}
$$

である。

さらに有限差分から、

$$
\boxed{
\dot S_{\chi,\mathrm{num}}
=
\frac{
S_\chi(n+1)-S_\chi(n)
}{
\Delta s
}
}
\tag{179}
$$

を記録する。

### 21.3 最小 normal form の直接検証

理論式、

$$
\dot S_\chi
=
\lambda J_\chi
+
gS_\chi(1-S_\chi^2)
$$

を数値時系列へ fit する。

すなわち、

$$
\boxed{
\dot S_{\chi,\mathrm{num}}
=
\lambda_{\mathrm{fit}}J_\chi
+
g_{\mathrm{fit}}
S_\chi(1-S_\chi^2)
+
\epsilon
}
\tag{180}
$$

として、

$$
\lambda_{\mathrm{fit}},
\qquad
g_{\mathrm{fit}}
$$

を推定する。

理論支持条件は、

$$
\boxed{
\lambda_{\mathrm{fit}}\neq0,
\qquad
g_{\mathrm{fit}}>0
}
\tag{181}
$$

である。

さらに residual $\epsilon$ が主要二項より十分小さいことを確認する。

### 21.4 mirror test

初期状態と全内部位相を複素共役した mirror run を作る。

元 run を、

$$
\mathcal R
$$

とし、mirror run を、

$$
\mathcal R^*
$$

とする。

理論上、

$$
\boxed{
J_\chi^{*}(n)
=
-J_\chi(n)
}
\tag{182}
$$

$$
\boxed{
\Delta\phi_\chi^{*}(n)
=
-\Delta\phi_\chi(n)
}
\tag{183}
$$

$$
\boxed{
S_\chi^{*}(n)
=
-S_\chi(n)
}
\tag{184}
$$

でなければならない。

一方、mirror-even な量、

$$
Q=p_L+p_R,
$$

総ノルム、

$$
|C|,
$$

局在度 $L$、

$$
N_{\mathrm{eff}}
$$

などは一致しなければならない。

従って mirror test の判定条件は、

$$
\boxed{
\max_n
|S_\chi^*(n)+S_\chi(n)|
<
\varepsilon_{\mathrm{mirror}}
}
\tag{185}
$$

および同様の $J_\chi,\Delta\phi_\chi$ 条件である。

### 21.5 対称初期条件 test

理論では、

$$
S_\chi=0,
\qquad
J_\chi=0
$$

は対称固定点である。

従って機械的に完全対称な初期条件では、浮動小数誤差を除けば一方が恣意的に選ばれてはならない。

次に、符号のみ異なる微小 seed、

$$
J_\chi(0)=+\epsilon
$$

と、

$$
J_\chi(0)=-\epsilon
$$

を用いる。

理論予測は、

$$
\boxed{
+\epsilon
\Longrightarrow
S_\chi\to+1
}
\tag{186}
$$

$$
\boxed{
-\epsilon
\Longrightarrow
S_\chi\to-1
}
\tag{187}
$$

である。

seed の絶対値を、

$$
10^{-k},
\qquad
k=2,\ldots,12
$$

程度まで下げ、符号だけが最終 chirality を決める領域を測定する。

### 21.6 位相閉包破壊 test

過去論文の因果鎖では、

$$
\text{整数倍音・位相閉包}
\to
J_\chi
$$

が必要である。

従って次の4条件を同一振幅分布で比較する。

1. 整数倍音＋位相閉包
2. 整数倍音＋位相ランダム化
3. 非整数周波数＋位相閉包
4. 非整数周波数＋位相ランダム化

理論支持条件は、条件1だけが長時間持続する、

$$
|J_\chi|>0
$$

と chirality selection を同時に示すことである。

### 21.7 因果順序 test

相互相関を、

$$
R_{J\Delta\phi}(\ell)
=
\operatorname{corr}
\left(
J_\chi(n),
\Delta\phi_\chi(n+\ell)
\right)
$$

$$
R_{\Delta\phi S}(\ell)
=
\operatorname{corr}
\left(
\Delta\phi_\chi(n),
S_\chi(n+\ell)
\right)
$$

とする。

理論上必要な順序は、

$$
\boxed{
J_\chi
\to
\Delta\phi_\chi
\to
S_\chi
\to
|S_\chi|\approx1
}
\tag{188}
$$

である。

正の lag で相関最大値が現れることを確認する。

### 21.8 $g$ の役割の分離

過去資料が既に強調していたように、

$$
\boxed{
J_\chi\text{ の生成}
\neq
S_\chi\text{ の増幅}
}
\tag{189}
$$

である。

従って、

$$
g<0,\quad g=0,\quad g>0
$$

を比較する。

- $g=0$：内部相関だけで符号付き小 $S_\chi$ が生成されるか。
- $g>0$：その seed が $\pm1$ へ増幅されるか。
- $g<0$：$S_\chi$ は0へ復元されるが $J_\chi$ 自体は残るか。

この分離が成立すれば、

$$
\boxed{
J_\chi
=
\text{chirality seed},
\qquad
gS_\chi(1-S_\chi^2)
=
\text{selection amplifier}
}
\tag{190}
$$

という理論解釈が直接検証される。

### 21.9 Bose 型対照

過去論文で Bose 型線形写像では一状態選択が成立しないことを対照とした。

同じ比較を chirality 読出しでも維持する。

Fermi 型 sector では、

$$
S_\chi\to\pm1
$$

が起こる一方、Bose 型対照では、

$$
S_\chi\approx0
$$

または周期的混合に留まることを確認する。

これにより chirality selection が単なる二チャネル交換一般の現象ではなく、Fermi 型内部閉包と非線形応答に結び付いていることを検証する。

### 21.10 公理保存監査――動力学として採用するための必須条件

本検証で最重要なのは、$S_\chi$ の選択が観測できることだけではない。**位相更新そのものが基礎公理を破壊していないことを同時に確認する。** 自己無撞着を基礎条件とする本系では、出力がもっともらしくても公理保存に失敗する更新則は動力学として採用しない。

最低限、各反復 $n$ について、

$$
\epsilon_C^{(n)}
=
\left|\sum_i X_i^{(n)2}\right|
$$

を記録し、零閉包誤差が数値精度内に留まることを確認する。さらに、有限回帰 $U^N=I$、simplex 閉包、ノルム、current 保存、mirror 対称性についても、既存コードで定義済みの監査量を更新前後で比較する。

判定原則は明瞭である。

$$
\boxed{
\text{観測上の選択が成立}
\;\land\;
\text{公理保存が成立}
}
$$

の場合だけを理論支持とする。選択が起きても零閉包その他の許容条件を破壊する run は棄却する。逆に、公理を保存したまま $J_\chi\to\Delta\phi_\chi\to S_\chi$ が再現されれば、**離散位相更新が許容状態空間内部の動力学として機能すること**の数値確認になる。

各 run で、

$$
\left|
\sum_iX_i^2
\right|
$$

を記録する。

さらに、

$$
p_L+p_R,
$$

全ノルム、

$$
\operatorname{div}_dJ,
$$

および本論文で導出した gauge current の保存誤差を記録する。

chirality selection が起きても、

$$
\boxed{
\text{零閉包}
}
$$

$$
\boxed{
\text{全体保存則}
}
$$

を破ってはならない。

### 21.11 既存コードの変更方針

既存の dynamics を変更しないことを原則とする。

変更は次の三種類に限定する。

1. **名称の再解釈**
   - `A -> L_sector`
   - `B -> R_sector`
   - `S -> S_chi`

2. **観測量の追加**
   - `C_cross`
   - `J_chi`
   - `delta_phi_chi`
   - `dS_chi`
   - `normal_form_residual`

3. **mirror run の追加**
   - 全複素状態の共役初期条件を生成
   - 同一 dynamics を適用
   - odd/even observable の変換則を比較

従って、理論を検証するために新しい相互作用項をコードへ追加しない。

### 21.12 判定基準

理論接続を数値的に支持するための必要条件を次とする。

1. mirror run で
   $$
   J_\chi\to-J_\chi,\quad
   \Delta\phi_\chi\to-\Delta\phi_\chi,\quad
   S_\chi\to-S_\chi
   $$
   が成立する。

2. mirror-even 保存量は一致する。

3. 整数倍音＋位相閉包条件でのみ持続的 $J_\chi$ が生じる。

4. $J_\chi$ が $\Delta\phi_\chi$ に先行または同時発生する。

5. $\Delta\phi_\chi$ が $S_\chi$ の成長に先行する。

6. normal form fit で
   $$
   \lambda_{\mathrm{fit}}\neq0,
   \qquad
   g_{\mathrm{fit}}>0
   $$
   が得られる。

7. seed の符号反転で最終 $S_\chi$ の符号が反転する。

8. Bose 型対照では同じ安定 chirality selection が成立しない。

9. 零閉包・ノルム・current 保存誤差が許容範囲に留まる。

これらが全て成立すれば、

$$
\boxed{
\text{過去の Fermi 型 A/B 選択系}
=
\text{Weyl chirality selection の数値実現}
}
\tag{191}
$$

として支持される。

### 21.13 後続論文への接続

本論文で理論導出は閉じた。

後続論文の主題は、新しい理論を追加することではなく、

$$
\boxed{
\text{既存数値実験系による chirality selection 導出の検証}
}
\tag{192}
$$

である。

後続論文では、

- mirror symmetry test,
- $J_\chi\to\Delta\phi_\chi\to S_\chi$ の因果順序,
- normal form fit,
- 位相閉包破壊,
- Bose/Fermi 対照,
- 零閉包保存監査

を実行し、本論文の理論接続が既存 dynamics 上で再現されるかを検証する。

## 22. 離散動力学から連続動力学へ

### 22.1 自己写像パラメータは物理時間ではない

本節の $s$ は物理時間を代用するパラメータではない。これは、零閉包を破壊せず自己無撞着解へ到達するための**構成・選択パラメータ**である。物理時間 $t$ は前論文で、他の軸と同じ複素軸から Lorentz 読出しによって現れる観測方向として扱った。従って本論文は「$s$ に沿う緩和流を物理時間発展と同一視する」構成ではない。

この区別は動力学の主張を弱めるものではなく、むしろ本公理系の核心である。求めるべきものは、外部時間に沿って状態を無条件に書き換える規則ではなく、**時間軸を含む全配置について公理を保つ自己無撞着解を選ぶ写像**である。その固定点を Lorentz 読出しで表現したとき、物理時間を含む偏微分場方程式として読まれる。

### 22.2 微小反復極限

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

この区別により、時間を基礎から特権化せずに、離散自己写像と標準的連続場方程式を両立できる。従って、$s$ に沿う relaxational dynamics が物理時間発展そのものだと解釈してはならない。$s$ は許容自己無撞着配置を構成する写像パラメータであり、物理時間 $t$ は固定点配置の Lorentz 読出しの内部に含まれる。

---

## 23. 導出鎖の完全形

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

## 24. 結論

前論文では、零閉包、有限回帰、simplex 閉包、自己無撞着から多数の対称構造を導出したが、Noether 保存則、局所 gauge dynamics、および chirality selection への接続が残っていた。

本論文ではまず、関係位相から離散 current と離散作用を構成し、零閉包接空間への射影と retraction により、

$$
\mathcal F_N:
\mathcal Z_N\to\mathcal Z_N
$$

という零閉包保存自己写像を構成した。

この結果は本論文の動力学に関する中心的結論である。**離散系では、位相 $\phi$ を更新することで状態を変化させながら、零閉包を破壊しない動力学を構成できる。** 多くの連続理論では状態更新が形式上先に与えられるが、本系では自己無撞着が基礎条件であるため、そのような暗黙の更新は許されない。更新則自身が

$$
\boxed{
\mathcal F_N(\mathcal Z_N)\subseteq\mathcal Z_N
}
$$

を満たすことを示して初めて、本公理系の動力学として採用できる。本論文では零閉包についてこれを解析的に閉じ、他の既導出条件についても数値検証時の必須監査条件とした。

有限 $N$ の離散 continuity equation は、

$$
N\to\infty
$$

の高分解能極限で、

$$
\partial_\mu J^\mu=0
$$

へ移り、連続作用は、

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

へ接続した。

さらに局所位相原点の無名性から辺 connection、simplex 面から curvature を得て、

$$
D_\mu
=
\partial_\mu-igA_\mu
$$

および、

$$
F_{\mu\nu}
=
\partial_\mu A_\nu-\partial_\nu A_\mu-ig[A_\mu,A_\nu]
$$

へ接続した。

前論文の5複素自由度と $3\oplus2$ 分解を合わせると、

$$
S(U(3)\times U(2))
$$

が現れ、trace-zero 条件

$$
3y_3+2y_2=0
$$

から hypercharge 比が固定される。

さらに、

$$
V^*\oplus\Lambda^2V
$$

は、

$$
d^c,\quad L,\quad u^c,\quad Q,\quad e^c
$$

という標準模型一世代15左手 Weyl 成分へ分解され、その表現について、

$$
SU(3)^3,
\quad
SU(3)^2U(1),
\quad
SU(2)^2U(1),
\quad
U(1)^3,
\quad
\mathrm{grav}^2U(1)
$$

の anomaly cancellation および $SU(2)$ global anomaly の不存在を直接検算した。

最後に、共役二 Weyl sector、

$$
\mathcal H_L
=
\left(\frac12,0\right)
\otimes
(V^*\oplus\Lambda^2V)
$$

$$
\mathcal H_R
=
\left(0,\frac12\right)
\otimes
(V\oplus\Lambda^2V^*)
$$

を、過去論文の A/B 二チャネルと同定した。

mirror 変換で、

$$
J\to-J,
\qquad
\Delta\phi\to-\Delta\phi,
\qquad
S_\chi\to-S_\chi
$$

となるため、対称点近傍の最小選択方程式は、

$$
\boxed{
\dot S_\chi
=
\lambda J
+
gS_\chi(1-S_\chi^2)
}
$$

となる。

$g>0$ では、

$$
S_\chi=0
$$

は不安定で、

$$
S_\chi=\pm1
$$

が安定である。

従って、基礎方程式に左右優先を入れずに、内部位相相関 $J$ の符号を seed として chirality が自発選択される。

以上から、本論文が対象とする範囲では、

$$
\boxed{
\text{零閉包}
\to
\text{Noether 保存則}
\to
\text{局所 gauge geometry}
\to
S(U(3)\times U(2))
\to
\text{hypercharge}
\to
\text{一世代表現}
\to
\text{anomaly cancellation}
\to
\text{chirality selection}
}
$$

という**標準模型の gauge-representation/chirality 構造への理論接続は閉じた**。

従って、この導出鎖について残る主課題は新しい理論項の追加ではなく、既存の数値実験系で、

$$
\boxed{
J_\chi
\to
\Delta\phi_\chi
\to
S_\chi
\to
\pm1
}
$$

が mirror symmetry、零閉包、保存則を保ったまま再現されることを確認する**数値検証**である。

後続論文では、既存 A/B Fermi 型実験系を変更せず、観測量の追加と mirror run を行い、

$$
\dot S_\chi
=
\lambda J_\chi
+
gS_\chi(1-S_\chi^2)
$$

の直接 fit、位相閉包破壊、Bose/Fermi 対照、因果順序、保存則監査を実施する。

なお、本論文でいう「標準模型との理論接続」は gauge group、local gauge geometry、一世代内部表現、hypercharge、anomaly cancellation、chirality selection の鎖を指す。Higgs 動径モード、Yukawa coupling、世代数、具体的質量階層、量子補正・renormalization は別の理論問題として残る。

## 参考文献

1. 木原範昭，「物理学の対称性は、本当に最初から与える必要があるのか」，前論文.
2. E. Noether, “Invariante Variationsprobleme,” *Nachrichten von der Gesellschaft der Wissenschaften zu Göttingen, Mathematisch-Physikalische Klasse*, 1918, 235–257.
3. J. E. Marsden and M. West, “Discrete Mechanics and Variational Integrators,” *Acta Numerica*, 10 (2001), 357–514. DOI: 10.1017/S096249290100006X.
4. V. A. Dorodnitsyn, “Noether-type theorems for difference equations,” *Applied Numerical Mathematics*, 39 (2001), 307–321. DOI: 10.1016/S0168-9274(00)00041-6.
5. M. Skopenkov, “Discrete Field Theory: Symmetries and Conservation Laws,” *Mathematical Physics, Analysis and Geometry*, 26 (2023), Article 19. DOI: 10.1007/s11040-023-09459-4.
6. S. Navas et al. (Particle Data Group), "Review of Particle Physics," *Phys. Rev. D* 110, 030001 (2024). DOI: 10.1103/PhysRevD.110.030001.
7. M. E. Peskin and D. V. Schroeder, *An Introduction to Quantum Field Theory*, Addison-Wesley (1995), Chapter 20 (gauge theories with spontaneous symmetry breaking; Standard Model representation content and anomaly cancellation).
8. K. G. Wilson, "Confinement of quarks," *Phys. Rev. D* 10, 2445 (1974). DOI: 10.1103/PhysRevD.10.2445 (link variables, plaquette action, continuum Yang–Mills limit).
9. H. Georgi and S. L. Glashow, "Unity of All Elementary-Particle Forces," *Phys. Rev. Lett.* 32, 438 (1974). DOI: 10.1103/PhysRevLett.32.438 (the $\overline{\mathbf 5}\oplus\mathbf{10}$ one-generation content).
10. 木原範昭，「零閉包・有限位数・自己無撞着幾何からの対称性生成――唯一の外部指定パラメータ $N$ と、一般化・動力学導出という残された課題」，Zenodo, Concept DOI: 10.5281/zenodo.22028072, Version DOI: 10.5281/zenodo.22028073 (2026)．（本文中の [1]・「前論文」）
