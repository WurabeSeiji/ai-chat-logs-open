# 万能相互作用 Einstein–Maxwell–GH–CF4 相互作用完全仕様 v2.0
## — admissible state ψ から K(ψ), S(ψ), ψ′ を一意に生成するための実装完全仕様 —

**著者:** 木原範昭 (Noriaki Kihara)  
**仕様整理:** ChatGPT  
**版:** v2.0  
**日付:** 2026-09-30  
**位置づけ:** `universal_interaction_einstein_maxwell_cf4_full_spec_ja_v1.0_20260930.md` のうち、相互作用生成に関する節を監査し、実装者判断が残らないように完全化する。  
**重要:** 本文書は **初期化仕様ではない**。初期化器が admissible state `ψ_n` を生成した後の、1 step の状態遷移だけを完全定義する。

---

# 0. 完成条件

本仕様で「相互作用仕様が完全」とは、次を意味する。

任意の admissible state

\[
\psi_n
\]

が与えられたとき、ケース名、質量、電荷、スピン、軌道情報などを外部引数として与えず、また実装者が追加の物理判断・数値判断を行わずに、

\[
\boxed{
\psi_n
\longrightarrow
K(\psi_n)
\longrightarrow
S(\psi_n)
\longrightarrow
\psi_{n+1}
}
\]

が一意に決まることをいう。

相互作用は

\[
\boxed{
\psi_{n+1}=S(\psi_n)\psi_n
}
\]

であり、

\[
\boxed{
S(\psi_n)=e^B e^A
}
\]

とする。

この仕様で未定義の外部係数・外部履歴・外部コントローラ・case switch を参照した時点で実装監査 FAIL とする。

---

# 1. v1.0 からの重要修正

v1.0 の相互作用部には以下の未完成点があった。

1. `∇_(a H_b)` の全 spacetime 成分の構成法が未定義だった。
2. `g,Pi,Phi` から 4D Christoffel、3D Christoffel、extrinsic curvature、`K` を作る順序と式が未固定だった。
3. Maxwell evolution の `E^i,B^i` を contravariant state とした場合の metric-derivative 項が named K へ完全には展開されていなかった。
4. Lie derivative の partial derivative と covariant derivative が混在し得る書き方だった。
5. canonical affine lift で、どの項を derivative block に入れ、どの項を `ONE` 列へ入れるかが完全には機械的に固定されていなかった。
6. `D_i`, modal multiplication, projection に必要な離散演算子の実体が `ψ` に完全には含まれていなかった。
7. matrix exponential action の実装方式が複数候補のままだった。
8. CF4 の式そのものは正しかったが、`K` と `k_i=hK_i` の区別を実装契約として明記していなかった。
9. `rac12` という転記誤りが GH local source に存在した。

v2.0 はこれらを修正する。

---

# 2. 物理系と符号規約

唯一の物理系は Einstein–Maxwell electrovac とする。

計量符号：

\[
(-,+,+,+).
\]

Einstein 方程式：

\[
G_{ab}=8\pi G T^{\rm EM}_{ab}.
\]

単位規約は interaction 内部で

\[
\boxed{c=1,\qquad 4\pi\epsilon_0=1}
\]

とする。`c` が state に保持されている場合、interaction 開始時に `c=1` を監査し、異なる値を dynamics へ silently 換算しない。`G` は state 値を Einstein coupling に使用する。

Maxwell 方程式：

\[
\nabla_aF^{ab}=0,
\qquad
\nabla_{[a}F_{bc]}=0.
\]

GH first-order variables：

\[
\boxed{
\Pi_{ab}=-n^c\partial_cg_{ab}
}
\]

\[
\boxed{
\Phi_{iab}=\partial_i g_{ab}
}
\]

extrinsic curvature の符号は

\[
\boxed{
K_{ij}=-\frac{1}{2}\mathcal L_n\gamma_{ij}
}
\]

とする。

---

# 3. 相互作用が要求する state contract

相互作用関数は次の state block だけを読む。

\[
\psi=
(\psi_{unit},\psi_{time},\psi_{phys0},\psi_{numerics},\psi_{grid},\psi_{field})^T.
\]

`ψ_phys0` は相互作用物理 block から読んではならない。

## 3.1 unit

\[
\boxed{ONE=1}
\]

## 3.2 dynamic fields

各 retained spatial coefficient index `I` について、

- `g_ab[I]` : 10
- `Pi_ab[I]` : 10
- `Phi_iab[I]` : 30
- `H_a[I]` : 4
- `Theta_a[I]` : 4
- `E^i[I]` : 3
- `B^i[I]` : 3

合計 64 field components / retained spatial coefficient。

## 3.3 相互作用に必要な数値恒等状態

最低限、次を `ψ_numerics` に恒等保持する。

\[
\gamma_0,\gamma_1,\gamma_2,\gamma_3,\gamma_4,\gamma_5,\epsilon_5,
\]

\[
\mu_H,\eta_H,\mu_L,\mu_S,p_H,
\]

\[
\kappa_\tau,D_\tau,
\]

\[
\epsilon_K,\epsilon_{exp},m_{exp,max},exp\_method\_id.
\]

v2.0 の必須値：

```text
exp_method_id = AL_MOHY_HIGHAM_2011_EXPMV
```

別方式へ切り替える場合は別 run で `exp_method_id` state を変更する。

---

# 4. 空間離散演算子を隠れ外部情報にしない

v1.0 では `basis_id`, `quadrature_id`, map coefficient から演算子を作る余地が実装者に残っていた。

v2.0 では、相互作用が直接使用する **完成済み semi-discrete operator** を `ψ_grid` の恒等状態に含める。

必須 operator state：

\[
\boxed{\mathbf D_i}\qquad i=1,2,3
\]

retained coefficient vectorに対する global physical-coordinate partial derivative operator。

また非線形積の aliasing 規則を一意にするため、

\[
\boxed{\mathbf V_*}
\]

retained coefficient → overintegration/evaluation point values、

\[
\boxed{\mathbf P_*}
\]

overintegration/evaluation point values → retained coefficient projection

を恒等 state として保持する。

要件：

\[
\mathbf P_*\mathbf V_*=I
\]

が retained space 上で設定 tolerance 内。

これらの生成法は別の grid/initialization 仕様で定めるが、**相互作用 step は生成法を知らない**。

`D_i,V_*,P_*` は domain interface、excision、outer treatment を含めて **semi-discrete spatial problem が閉じた完成済み global operator** でなければならない。interaction 中に外部 boundary value、外部 penalty、外部 ghost state を参照してはならない。外部値を必要とする operator state は admissible state ではない。

相互作用は `ψ_grid` に入っている具体的行列だけを読む。

したがって同じ `ψ_n` からは同じ `K(ψ_n)` が一意に出る。

---

# 5. 離散微分・乗算の唯一の定義

任意の retained coefficient vector `u` に対して、

\[
\boxed{
\partial_i u := \mathbf D_i u
}
\]

とする。

stage state `ψ_*` から評価した scalar/tensor coefficient field `f(x;ψ_*)` による frozen multiplication operator は、

\[
\boxed{
\mathscr M_*[f]
=\mathbf P_*\,\mathrm{diag}(f_*)\,\mathbf V_*
}
\]

とする。

ここで `f_*` は overintegration/evaluation points 上で `ψ_*` から計算した値。

すべての nonlinear product は **必ず evaluation space で形成してから `P_*` で retained space に戻す**。

retained coefficients 上で直接 component-wise product してはならない。

local source projection は

\[
\boxed{
\mathcal P_*[R]
:=\mathbf P_*R_*
}
\]

とする。

---

# 6. stage state から導出する 3+1 幾何量

以下は毎 `K(ψ_*)` 評価時に stage state から作る一時量であり、永続 state ではない。

spatial metric：

\[
\boxed{\gamma_{ij}=g_{ij}}
\]

inverse spatial metric：

\[
\gamma^{ik}\gamma_{kj}=\delta^i{}_j.
\]

shift covector：

\[
\boxed{\beta_i=g_{0i}}
\]

shift vector：

\[
\boxed{\beta^i=\gamma^{ij}\beta_j}
\]

lapse：

\[
\boxed{
\alpha=\sqrt{-g_{00}+\beta_i\beta^i}
}
\]

正の branch のみを採用する。

\[
\alpha\le0
\]

または平方根中身が非正なら run-time FAIL。

normal covector/vector：

\[
\boxed{n_a=(-\alpha,0,0,0)}
\]

\[
\boxed{
n^a=\left(\alpha^{-1},-\alpha^{-1}\beta^i\right)
}
\]

spatial determinant：

\[
\gamma=\det(\gamma_{ij})>0.
\]

---

# 7. metric の全 spacetime partial derivative

4D Christoffel を作るため、metric の coordinate partial derivativeを次で一意に定義する。

spatial derivative：

\[
\boxed{
\partial_i g_{ab}=\Phi_{iab}
}
\]

time derivative は `Pi` の定義から直接、

\[
\Pi_{ab}
=-\frac1\alpha
(\partial_tg_{ab}-\beta^i\partial_i g_{ab})
\]

ゆえに

\[
\boxed{
\partial_0g_{ab}
=\beta^i\Phi_{iab}-\alpha\Pi_{ab}
}
\]

とする。

ここでは GH evolution equation を代入してはならない。

`∂_0 g` は `Pi` の定義から作る。

---

# 8. 4D Christoffel と GH constraint

first-kind Christoffel：

\[
\boxed{
\Gamma_{abc}
=\frac12(
\partial_b g_{ac}
+\partial_c g_{ab}
-\partial_a g_{bc})
}
\]

second-kind：

\[
\boxed{
\Gamma^a{}_{bc}=g^{ad}\Gamma_{dbc}
}
\]

contracted lower-index Christoffel：

\[
\boxed{
\Gamma_a=g^{bc}\Gamma_{abc}
}
\]

GH gauge constraint：

\[
\boxed{
\mathcal C_a=H_a+\Gamma_a
}
\]

three-index first-order constraint：

\[
\boxed{
\mathcal C_{iab}=\partial_i g_{ab}-\Phi_{iab}
}
\]

semi-discreteでは `∂_i g` に `D_i g` を用いる。

---

# 9. 3D Christoffel、extrinsic curvature、trace K

spatial metric derivative：

\[
\boxed{
\partial_k\gamma_{ij}=\Phi_{kij}
}
\]

3D Christoffel：

\[
\boxed{
{}^{(3)}\Gamma^i{}_{jk}
=\frac12\gamma^{il}
(\Phi_{jlk}+\Phi_{klj}-\Phi_{ljk})
}
\]

extrinsic curvature は GH variables から

\[
\boxed{
K_{ij}
=\frac12\Pi_{ij}
+\frac12 n^a(\Phi_{ija}+\Phi_{jia})
}
\]

とする。

trace：

\[
\boxed{K=\gamma^{ij}K_{ij}}
\]

---

# 10. derived lapse/shift derivative

Maxwell local source に必要な metric derivative を一意にする。

\[
\partial_k\gamma^{ij}
=-\gamma^{im}\gamma^{jn}\Phi_{kmn}.
\]

\[
\partial_k\beta_i=\Phi_{k0i}.
\]

\[
\boxed{
\partial_k\beta^i
=(\partial_k\gamma^{ij})\beta_j
+\gamma^{ij}\Phi_{k0j}
}
\]

lapse derivative：

\[
2\alpha\partial_k\alpha
=-\Phi_{k00}
+(\partial_k\beta_i)\beta^i
+\beta_i\partial_k\beta^i.
\]

すなわち

\[
\boxed{
\partial_k\alpha
=\frac{-\Phi_{k00}+\Phi_{k0i}\beta^i+\beta_i\partial_k\beta^i}{2\alpha}
}
\]

とする。

---

# 11. gauge target の完全定義

Damped-wave target：

\[
\boxed{
F_a
=\mu_L\log\left(\frac{\gamma^{p_H}}{\alpha}\right)n_a
-\mu_S\alpha^{-1}g_{ai}\beta^i
}
\]

ここで `γ=det γ_ij`。

引数が非正で log が定義不能なら run-time FAIL。

v2.0 では gauge blend を使用しない。

`T_gauge-blend` は相互作用物理 blockから読まない。

将来 blend を使う場合は別仕様版とする。

---

# 12. gauge-driver RHS と ∇aHb の完全構成

唯一の gauge-driver は

\[
\boxed{
\partial_t H_a-\beta^k\partial_kH_a
=-\mu_H(H_a-F_a)+\Theta_a
}
\]

\[
\boxed{
\partial_t\Theta_a
+\eta_H\beta^k\partial_kH_a
=-\eta_H\Theta_a
}
\]

である。

したがって現在 stage における coordinate derivative は

\[
\boxed{
\partial_0H_a
=\beta^k\partial_kH_a
-\mu_H(H_a-F_a)+\Theta_a
}
\]

\[
\boxed{
\partial_iH_a=\mathbf D_iH_a
}
\]

と一意に決める。

これを使い、

\[
\boxed{
\nabla_aH_b
=\partial_aH_b-\Gamma^c{}_{ab}H_c
}
\]

\[
\boxed{
\nabla_{(a}H_{b)}
=\frac12(\nabla_aH_b+\nabla_bH_a)
}
\]

を構成する。

`∂_t H` を別の predictor から取ってはならない。

---

# 13. Maxwell state から F_ab を一意に再構成

state は contravariant spatial vectors `E^i,B^i`。

4D volume form の向きは座標 `(t,x,y,z)` に対して

\[
\boxed{\epsilon_{0123}=+\sqrt{-g}}
\]

と固定する。

spatial vector の4D contravariant extension は

\[
\boxed{E^a=(0,E^i),\qquad B^a=(0,B^i)}
\]

とする。これは `n_aE^a=n_aB^a=0` を満たす。

\[
E_i=\gamma_{ij}E^j,
\qquad
B_i=\gamma_{ij}B^j.
\]

4D spatial extensions は

\[
E_a=(\beta^iE_i,E_i),
\qquad
B_a=(\beta^iB_i,B_i)
\]

とし、`E_a n^a=B_a n^a=0` を満たす。

spatial Levi-Civita tensor：

\[
\epsilon_{ijk}=\sqrt\gamma\,[ijk],
\qquad
\epsilon^{ijk}=\frac{1}{\sqrt\gamma}[ijk].
\]

4D Maxwell tensor：

\[
\boxed{
F_{ab}
=n_aE_b-n_bE_a+\epsilon_{abc}B^c
}
\]

ここで

\[
\epsilon_{abc}=n^d\epsilon_{dabc}
\]

を使用する。

この定義により

\[
E_a=F_{ab}n^b
\]

を満たす符号規約を固定する。

---

# 14. Maxwell stress-energy の唯一の構成

primary definition は 4D tensor から直接作る。

\[
\boxed{
T^{\rm EM}_{ab}
=\frac1{4\pi}
\left(
F_{ac}F_b{}^c
-\frac14g_{ab}F_{cd}F^{cd}
\right)
}
\]

trace

\[
T^{EM}=g^{ab}T^{EM}_{ab}
\]

を毎 stage 計算し、

\[
|T^{EM}|<\epsilon_T
\]

を監査する。Einstein RHS では理論上の `T=0` を代入して省略せず、必ず `T_ab-1/2 g_ab T` を形成する。

3+1 quantities は監査用に

\[
\rho_{EM}=\frac{E_iE^i+B_iB^i}{8\pi},
\]

\[
S_i^{EM}=\frac1{4\pi}\epsilon_{ijk}E^jB^k,
\]

\[
S_{ij}^{EM}=\frac1{4\pi}
\left[
-E_iE_j-B_iB_j
+\frac12\gamma_{ij}(E^2+B^2)
\right]
\]

も計算してよいが、GH matter source の正本は 4D `T_ab` とする。

---

# 15. GH evolution — 唯一の continuous RHS

mesh velocity は v2.0 で 0。

## 15.1 g_ab

\[
\boxed{
\partial_tg_{ab}
=(1+\gamma_1)\beta^k\partial_kg_{ab}
-\alpha\Pi_{ab}
-\gamma_1\beta^i\Phi_{iab}
}
\]

## 15.2 Pi_ab

\[
\boxed{
\begin{aligned}
\partial_t\Pi_{ab}
={}&\beta^k\partial_k\Pi_{ab}
-\alpha\gamma^{ki}\partial_k\Phi_{iab}
+\gamma_1\gamma_2\beta^k\partial_kg_{ab}\\
&+R^{GH}_{\Pi,ab}
-16\pi G\alpha\left(T^{EM}_{ab}-\frac12g_{ab}T^{EM}\right).
\end{aligned}
}
\]

\[
\boxed{
\begin{aligned}
R^{GH}_{\Pi,ab}
={}&2\alpha g^{cd}
\left(
\gamma^{ij}\Phi_{ica}\Phi_{jdb}
-\Pi_{ca}\Pi_{db}
-g^{ef}\Gamma_{ace}\Gamma_{bdf}
\right)\\
&-2\alpha\nabla_{(a}H_{b)}
-\frac12\alpha n^cn^d\Pi_{cd}\Pi_{ab}
-\alpha n^c\Pi_{ci}\gamma^{ij}\Phi_{jab}\\
&+\alpha\gamma_0
\left(2\delta^c{}_{(a}n_{b)}-(1+\gamma_3)g_{ab}n^c\right)\mathcal C_c\\
&+2\gamma_4\alpha\Pi_{ab}n^c\mathcal C_c\\
&-\gamma_5\alpha n^c\mathcal C_c
\frac{\mathcal C_a\mathcal C_b-\frac12g_{ab}\mathcal C_d\mathcal C^d}
{\epsilon_5+2(n^d\mathcal C_d)^2+\mathcal C_d\mathcal C^d}\\
&-\gamma_1\gamma_2\beta^i\Phi_{iab}.
\end{aligned}
}
\]

baseline では

\[
\gamma_3=\gamma_4=\gamma_5=0.
\]

## 15.3 Phi_iab

\[
\boxed{
\begin{aligned}
\partial_t\Phi_{iab}
={}&\beta^k\partial_k\Phi_{iab}
-\alpha\partial_i\Pi_{ab}
+\alpha\gamma_2\partial_i g_{ab}\\
&+\frac12\alpha n^cn^d\Phi_{icd}\Pi_{ab}
+\alpha\gamma^{jk}n^c\Phi_{ijc}\Phi_{kab}
-\alpha\gamma_2\Phi_{iab}.
\end{aligned}
}
\]

---

# 16. Maxwell evolution — contravariant state に対する完全展開

採用する electrovac equations：

\[
(\partial_t-\mathcal L_\beta)E^i
=\epsilon^{ijk}D_j(\alpha B_k)+\alpha K E^i,
\]

\[
(\partial_t-\mathcal L_\beta)B^i
=-\epsilon^{ijk}D_j(\alpha E_k)+\alpha K B^i.
\]

Lie derivative は spatial vector に対して

\[
\boxed{
(\mathcal L_\beta E)^i
=\beta^j\partial_jE^i-E^j\partial_j\beta^i
}
\]

と **partial derivative で**評価する。

curl は connection cancellation を使い、

\[
\epsilon^{ijk}D_j(\alpha B_k)
=\epsilon^{ijk}\partial_j(\alpha\gamma_{kl}B^l)
\]

とする。

従って完全展開は

\[
\boxed{
\begin{aligned}
\partial_tE^i
={}&\beta^j\partial_jE^i
+\alpha\epsilon^{ijk}\gamma_{kl}\partial_jB^l\\
&-E^j\partial_j\beta^i
+\epsilon^{ijk}\partial_j(\alpha\gamma_{kl})B^l
+\alpha K E^i.
\end{aligned}
}
\]

\[
\boxed{
\begin{aligned}
\partial_tB^i
={}&\beta^j\partial_jB^i
-\alpha\epsilon^{ijk}\gamma_{kl}\partial_jE^l\\
&-B^j\partial_j\beta^i
-\epsilon^{ijk}\partial_j(\alpha\gamma_{kl})E^l
+\alpha K B^i.
\end{aligned}
}
\]

ここで

\[
\boxed{
\partial_j(\alpha\gamma_{kl})
=(\partial_j\alpha)\gamma_{kl}+\alpha\Phi_{jkl}
}
\]

であり、v1.0 で不足していた metric-derivative contribution を明示的に含める。

Maxwell constraints：

\[
\boxed{
\mathcal C_E=D_iE^i
=\partial_iE^i+{}^{(3)}\Gamma^i{}_{ik}E^k
}
\]

\[
\boxed{
\mathcal C_B=D_iB^i
=\partial_iB^i+{}^{(3)}\Gamma^i{}_{ik}B^k
}
\]

v2.0 では divergence-cleaning field は追加しない。

---

# 17. canonical affine lift の唯一の規則

stage state を `ψ_*` とする。

各 continuous RHS を

\[
F_A(\psi_*)
\]

とする。

`K(ψ_*)` は arbitrary augmented vector `v` に作用する frozen affine operator として定義する。

規則は以下だけ。

### Rule L

**retained dynamic field に spatial derivative が直接作用している項だけ**を corresponding dynamic column block に入れる。

### Rule ONE

それ以外の全 algebraic/local term は、線形・非線形を問わず stage state `ψ_*` で完全評価し、`ONE` 列に入れる。

従って

\[
\boxed{
(K(\psi_*)v)_A
=\sum_B L_A{}^B(\psi_*)v_B
+R_A(\psi_*)v_{ONE}
}
\]

である。

`v_ONE` は exponential action 中の augmented coordinate。

この規則により factorization の自由度を禁止する。

---

# 18. K の named derivative blocks

以下だけが dynamic-field derivative column として非零。

## 18.1 g row

\[
\boxed{
K[g,g]
=\mathscr M[(1+\gamma_1)\beta^k]\mathbf D_k
}
\]

`-αPi-γ1βPhi` は `ONE` 列。

## 18.2 Pi row

\[
\boxed{
K[Pi,Pi]=\mathscr M[\beta^k]\mathbf D_k
}
\]

\[
\boxed{
K[Pi,Phi_i]
=-\mathscr M[\alpha\gamma^{ki}]\mathbf D_k
}
\]

\[
\boxed{
K[Pi,g]
=\mathscr M[\gamma_1\gamma_2\beta^k]\mathbf D_k
}
\]

## 18.3 Phi_i row

\[
\boxed{
K[Phi_i,Phi_i]
=\mathscr M[\beta^k]\mathbf D_k
}
\]

\[
\boxed{
K[Phi_i,Pi]
=-\mathscr M[\alpha]\mathbf D_i
}
\]

\[
\boxed{
K[Phi_i,g]
=\mathscr M[\alpha\gamma_2]\mathbf D_i
}
\]

## 18.4 H row

\[
\boxed{
K[H,H]
=\mathscr M[\beta^k]\mathbf D_k
}
\]

## 18.5 Theta row

\[
\boxed{
K[Theta,H]
=-\mathscr M[\eta_H\beta^k]\mathbf D_k
}
\]

## 18.6 E row

\[
\boxed{
K[E^i,E^j]
=\delta^i{}_j\,\mathscr M[\beta^k]\mathbf D_k
}
\]

\[
\boxed{
K[E^i,B^l]
=\mathscr M[\alpha\epsilon^{ijk}\gamma_{kl}]\mathbf D_j
}
\]

## 18.7 B row

\[
\boxed{
K[B^i,B^j]
=\delta^i{}_j\,\mathscr M[\beta^k]\mathbf D_k
}
\]

\[
\boxed{
K[B^i,E^l]
=-\mathscr M[\alpha\epsilon^{ijk}\gamma_{kl}]\mathbf D_j
}
\]

## 18.8 time row

\[
\boxed{
K[T_\tau,T_\tau]=\kappa_\tau
}
\]

---

# 19. K の ONE 列 — 完全定義

各 field row の `ONE` column は continuous RHS から前節 derivative terms を差し引いた残りを **stage state で評価した retained coefficients** とする。

## g

\[
\boxed{
R_g=-\alpha\Pi_{ab}-\gamma_1\beta^i\Phi_{iab}
}
\]

## Pi

\[
\boxed{
R_\Pi
=R^{GH}_{\Pi,ab}
-16\pi G\alpha\left(T^{EM}_{ab}-\frac12g_{ab}T^{EM}\right)
}
\]

## Phi

\[
\boxed{
R_{\Phi,iab}
=\frac12\alpha n^cn^d\Phi_{icd}\Pi_{ab}
+\alpha\gamma^{jk}n^c\Phi_{ijc}\Phi_{kab}
-\alpha\gamma_2\Phi_{iab}
}
\]

## H

\[
\boxed{
R_{H,a}=-\mu_H(H_a-F_a)+\Theta_a
}
\]

## Theta

\[
\boxed{
R_{\Theta,a}=-\eta_H\Theta_a
}
\]

## E

\[
\boxed{
R_E^i
=-E^j\partial_j\beta^i
+\epsilon^{ijk}\partial_j(\alpha\gamma_{kl})B^l
+\alpha K E^i
}
\]

## B

\[
\boxed{
R_B^i
=-B^j\partial_j\beta^i
-\epsilon^{ijk}\partial_j(\alpha\gamma_{kl})E^l
+\alpha K B^i
}
\]

そして

\[
\boxed{
K[A,ONE]=\mathcal P_*[R_A]
}
\]

とする。

---

# 20. K の完全 block type 一覧

非零になり得る field block type は v2.0 でも以下の21種だけ。

1. `K[g,g]`
2. `K[g,ONE]`
3. `K[Pi,Pi]`
4. `K[Pi,Phi]`
5. `K[Pi,g]`
6. `K[Pi,ONE]`
7. `K[Phi,Phi]`
8. `K[Phi,Pi]`
9. `K[Phi,g]`
10. `K[Phi,ONE]`
11. `K[H,H]`
12. `K[H,ONE]`
13. `K[Theta,H]`
14. `K[Theta,ONE]`
15. `K[E,E]`
16. `K[E,B]`
17. `K[E,ONE]`
18. `K[B,B]`
19. `K[B,E]`
20. `K[B,ONE]`
21. `K[T_tau,T_tau]`

ただし v2.0 では **各 block の内容を前節までの式で一意化した**。

---

# 21. 恒等 state と禁止 coupling

任意の identity state `C`：

\[
\boxed{K[C,B]=0\quad\forall B}
\]

初期化専用 state `I0` から任意の dynamic field row `u` へ

\[
\boxed{K[u,I0]=0}
\]

を要求する。

従って dynamics は `M^(0),Q^(0),J^(0),C0,D0,case_id` を一切読まない。

---

# 22. direct RHS evaluator F_direct

`F_direct(ψ)` は第15,16,12節の continuous equations を semi-discrete operator `D_i,V_*,P_*` で直接評価する独立関数。

`F_direct` は `K` の block assembly routine を呼んではならない。

`K` 側も `F_direct` を呼んではならない。

毎 CF4 stage で

\[
\boxed{
\|K(\psi)\psi-F_{direct}(\psi)\|_\infty<\epsilon_K
}
\]

を要求する。

不一致ならその stage で停止。

---

# 23. CF4 — 唯一の time update

`K_j=K(ψ^[j])` とし、`k_j=hK_j` と定義する。

\[
h=\Delta\tau=\frac{\log D_\tau}{\kappa_\tau}.
\]

`κ_tau=0` は v2.0 では不許可。

stage 1：

\[
\psi^{[1]}=\psi_n,
\qquad K_1=K(\psi^{[1]}).
\]

stage 2：

\[
\boxed{
\psi^{[2]}=e^{\frac12hK_1}\psi_n
}
\]

\[
K_2=K(\psi^{[2]}).
\]

stage 3：

\[
\boxed{
\psi^{[3]}=e^{\frac12hK_2}\psi_n
}
\]

\[
K_3=K(\psi^{[3]}).
\]

stage 4：

\[
\boxed{
\psi^{[4]}
=e^{h(K_3-\frac12K_1)}\psi^{[2]}
}
\]

\[
K_4=K(\psi^{[4]}).
\]

final exponent 1：

\[
\boxed{
A=\frac{h}{12}(3K_1+2K_2+2K_3-K_4)
}
\]

\[
\psi_{n+1/2}=e^A\psi_n.
\]

final exponent 2：

\[
\boxed{
B=\frac{h}{12}(-K_1+2K_2+2K_3+3K_4)
}
\]

\[
\boxed{
\psi_{n+1}=e^B\psi_{n+1/2}
}
\]

したがって

\[
\boxed{
S(\psi_n)=e^Be^A
}
\]

である。

積順序は固定。

---

# 24. matrix exponential action の唯一の実装契約

v2.0 で許可する exponential action は

```text
AL_MOHY_HIGHAM_2011_EXPMV
```

のみ。

対象は

\[
y=e^{A}x
\]

であり、dense `e^A` を materialize する必要はない。

実装は Al-Mohy–Higham 2011 の scaling + truncated Taylor action algorithm に従う。

state から読む値：

- `epsilon_exp`
- `m_exp_max`
- method parameter table / theta table version

method parameter table も再現性 state として保持する。

各 action は

- selected Taylor degree `m`
- scaling count `s`
- estimated backward error
- final residual estimate

をログする。

\[
error_{exp}>\epsilon_{exp}
\]

なら step FAIL。

別 exponential algorithm を silently fallback してはならない。

---

# 25. CF4 stage ごとの評価順序

各 `K(ψ_*)` は必ず次の順序で評価する。

1. identity/numerical/grid state validity check
2. `V_*` で stage fields を evaluation points へ展開
3. inverse metric, lapse, shift, normal, determinant
4. metric spacetime derivatives
5. 4D Christoffel, 3D Christoffel
6. extrinsic curvature, trace K
7. GH constraints
8. gauge target `F_a`
9. gauge-driver `∂0H`
10. `∇_(aH_b)`
11. Maxwell `F_ab`
12. Maxwell stress-energy `T_ab`
13. continuous direct RHS
14. derivative blocks `L`
15. local residuals `R`
16. `K` operator object生成
17. `Kψ=F_direct` audit

この順序を変えて前 stage の一時量を再利用してはならない。

---

# 26. admissibility / run-time FAIL conditions

interaction step は次のいずれかで停止する。

1. `ONE != 1` beyond tolerance
2. `c != 1` beyond unit-convention tolerance
3. `D_tau<=0` または `T_tau<=0`
4. identity state drift
5. `det(g_ab)=0` または inverse failure
6. `γ<=0`
7. `α<=0`
8. damped-wave log argument `γ^p/α<=0`
9. NaN/Inf
10. `Kψ-F_direct` audit failure
11. exponential action tolerance failure
12. metric-signature audit failure
13. Maxwell stress tensor trace audit failure
14. spatial operator state inconsistency (`P_*V_* != I` beyond tolerance)
15. semi-discrete operator が外部 boundary / ghost / penalty 値を要求

FAIL 時に state を補正して続行してはならない。

---

# 27. interaction が絶対に読んではならないもの

- `case_id`
- current readout mass
- current readout charge
- current readout spin
- orbit radius readout
- GW readout
- EM radiation readout
- horizon finder history
- previous CF4 stage from prior step
- adaptive step history
- external configuration file value not copied into state
- initial mass/charge/spin labels as dynamics coefficient

---

# 28. interaction が生成してよい一時量

- inverse metrics
- lapse/shift/normal
- Christoffel symbols
- extrinsic curvature
- constraints
- `F_a`
- `∂aH_b`
- `T_ab`
- multiplication operators
- `K_1..K_4`
- CF4 stage states
- exponential work vectors

これらは step 終了後に破棄する。

---

# 29. interaction output contract

1 step 成功時に永続的に更新されるのは

- `T_tau`
- `g_ab`
- `Pi_ab`
- `Phi_iab`
- `H_a`
- `Theta_a`
- `E^i`
- `B^i`

だけ。

すべての identity state は bitwise または指定 tolerance 内で不変でなければならない。

---

# 30. 数学的固定点監査

もし stage state が連続方程式の exact stationary state で

\[
F_{direct}(\psi_*)=0
\]

なら

\[
K(\psi_*)\psi_*=0.
\]

従って理想的 exact exponential では

\[
S(\psi_*)\psi_*=\psi_*.
\]

有限演算ではその差を exponential tolerance と spatial truncation error で監査する。

この fixed-point test を interaction 単体の必須 validation とする。

---

# 31. v1.0 に対する決定的修正点

特に Maxwell block は v1.0 の

\[
R_E^i=-E^jD_j\beta^i+\epsilon^{ijk}(D_j\alpha)B_k+\alpha K E^i
\]

という省略形を廃止する。

contravariant state `B^i` を使用する以上、正本は

\[
\boxed{
R_E^i
=-E^j\partial_j\beta^i
+\epsilon^{ijk}[(\partial_j\alpha)\gamma_{kl}+\alpha\Phi_{jkl}]B^l
+\alpha K E^i
}
\]

である。

同様に

\[
\boxed{
R_B^i
=-B^j\partial_j\beta^i
-\epsilon^{ijk}[(\partial_j\alpha)\gamma_{kl}+\alpha\Phi_{jkl}]E^l
+\alpha K B^i
}
\]

とする。

これにより curved spatial metric における metric derivative contribution を欠落させない。

---

# 32. 相互作用仕様の自己完結性

本仕様だけで以下は一意に実装できる。

- stage state から幾何量生成
- full GH RHS
- gauge-driver RHS
- full electrovac Maxwell RHS
- Maxwell stress-energy backreaction
- canonical affine `K(ψ)`
- named K block content
- CF4 stage
- exponential action method selection
- final `S(ψ)` action
- interaction-level audit / fail conditions

一方、以下は本仕様の対象外であり、別仕様が必要。

1. **初期化:** superposed Kerr–Newman + XCTS + electromagnetic constraints から admissible `ψ_0` を作る方法
2. **grid/operator construction:** `D_i,V_*,P_*` をどの domain/basis から生成するか
3. **readout:** horizon/GW/EM flux 等の観測量生成

ただし 2 の生成法は対象外でも、生成済み operator 自体は interaction state の必須成分なので、interaction に隠れ自由度は残らない。

---

# 33. 実装開始許可条件

相互作用コードの実装開始前に、次を文書監査で全 PASS する。

- [ ] 64 dynamic components の全 RHS が本仕様に存在
- [ ] 全 derived geometry 式が存在
- [ ] `∇_(aH_b)` の time component が存在
- [ ] Maxwell metric-derivative term が存在
- [ ] `K` factorization rule が一意
- [ ] 21 named block type が一意
- [ ] `Kψ=F_direct` 独立監査が定義済み
- [ ] CF4 stage と積順序が固定
- [ ] exponential method が固定
- [ ] external case switch が不存在
- [ ] initialization labels から field rows への coupling がゼロ

一つでも未定義なら相互作用コードを書いてはならない。

---

# 34. 参照標準

相互作用の標準物理式の確認に用いた主な参照：

1. Lindblom, Scheel, Kidder, Owen, Rinne, first-order generalized harmonic Einstein system, Class. Quantum Grav. 23 (2006) S447–S462. 特に first-order GH equations (35)–(37).
2. Lindblom & Szilágyi, “An Improved Gauge Driver for the Generalized Harmonic Einstein System”, Phys. Rev. D 80, 084019 (2009). Gauge driver Eqs. (9),(11), damped-wave target Eq. (A15).
3. SpECTRE `gh::TimeDerivative` documentation: full GH RHS and matter source conventions.
4. SpECTRE generalized-harmonic extrinsic-curvature helper: `K_ij = 1/2 Pi_ij + Phi_(ij)a n^a`.
5. Standard 3+1 Maxwell electrovac equations with lapse/shift/extrinsic-curvature coupling.
6. Celledoni, Marthinsen, Owren, “Commutator-free Lie group methods”, Future Generation Computer Systems 19 (2003) 341–352. Fourth-order commutator-free scheme used here.
7. Al-Mohy & Higham, “Computing the Action of the Matrix Exponential, with an Application to Exponential Integrators”, SIAM J. Sci. Comput. 33 (2011) 488–511.

---

# 35. 最終定義

v2.0 における万能相互作用は、

\[
\boxed{
\mathcal U:
\psi_n
\mapsto
K_1,K_2,K_3,K_4
\mapsto
A,B
\mapsto
S(\psi_n)=e^Be^A
}
\]

である。

そして

\[
\boxed{
\psi_{n+1}=S(\psi_n)\psi_n
}
\]

以外の time evolution route を認めない。

**本仕様 v2.0 を満たさない相互作用実装は、数値結果が見かけ上もっともらしくても正式実験には使用しない。**
