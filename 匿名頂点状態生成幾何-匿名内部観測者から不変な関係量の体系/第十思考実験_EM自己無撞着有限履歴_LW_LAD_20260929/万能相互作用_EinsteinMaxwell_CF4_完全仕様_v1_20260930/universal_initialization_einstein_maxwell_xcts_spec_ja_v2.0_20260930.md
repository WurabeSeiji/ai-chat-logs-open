# 万能相互作用 Einstein–Maxwell 初期化仕様 v2.0
## — superposed Kerr–Newman + XCTS + electromagnetic constraints から interaction v2.0 admissible state ψ₀ を一意に生成するための改訂仕様 —

**著者:** 木原範昭 (Noriaki Kihara)  
**仕様整理:** ChatGPT  
**版:** v2.0  
**日付:** 2026-09-30  
**上位相互作用仕様:** `universal_interaction_operator_complete_spec_ja_v2.0_20260930.md`  
**旧仕様:** `universal_interaction_einstein_maxwell_cf4_full_spec_ja_v1.0_20260930.md` の初期化部を置換する。

---

# 0. 本文書の目的と完成条件

本仕様の目的は、初期物理条件と数値設定から、相互作用仕様 v2.0 が要求する admissible state

\[
\boxed{\psi_0}
\]

を一意に生成することである。

初期化写像を

\[
\boxed{
\mathcal I:
(\psi_{\rm phys0},\psi_{\rm initcfg},\psi_{\rm grid})
\longmapsto
\psi_0
}
\]

と定義する。

ここで `ψ_grid` は interaction v2.0 の state contract を満たす完成済み grid/operator state であり、少なくとも

\[
\mathbf D_i,\quad \mathbf V_*,\quad \mathbf P_*
\]

を含む。

**本仕様で「初期化仕様が完全」とは、承認済み `ψ_grid` が与えられたとき、実装者が追加の物理判断・数値判断・ケース別の代替初期化を行わずに、唯一の `ψ₀` または明示的 FAIL を返せることをいう。**

`ψ_grid` 自体の生成法は interaction v2.0 が別仕様に分離しているため、本書ではその生成アルゴリズムを勝手に再定義しない。ただし、本初期化器が受理可能な `ψ_grid` の契約と監査条件は本書で完全に定める。

---

# 1. 絶対ルール

1. 初期化法は全ケースで同じとする。
2. `case_id` による MP / RN / KN / vacuum / opposite-charge などの初期化アルゴリズム切替を禁止する。
3. 個々のブラックホール free data は常に Kerr–Newman in Kerr–Schild form を同一式で生成する。
4. spin=0 は Kerr–Newman の `a=0` 極限として扱う。
5. charge=0 は Kerr–Newman の `Q=0` 極限として扱う。
6. `|Q|/M=1` は同じ Kerr–Newman/XCTS 初期化器へ通す。Majumdar–Papapetrou 厳密解へ勝手に置換しない。
7. `|Q|/M>1` または Kerr–Newman bound を破る入力は救済せず reject する。
8. Einstein–Maxwell constraints を満たさない近似 superposition を `ψ₀` として採用しない。
9. XCTS / Maxwell constraint solve が収束しない場合、別の初期化器へ fallback しない。
10. 初期化後の dynamics は初期ラベル `M^(0),Q^(0),J^(0)` を読まない。全物理情報は `ψ_field` に写像済みでなければならない。
11. 初期化に使用した設定・solver parameter・attenuation parameter・target parameter はすべて恒等 state として保存する。
12. 初期化で使う spatial derivative / projection は interaction と同じ `D_i,V_*,P_*` を使用する。別の隠れ derivative operator を使用しない。

---

# 2. 単位・符号・添字規約

interaction v2.0 と完全一致させる。

計量符号：

\[
(-,+,+,+).
\]

初期化 v2.0 の単位系：

\[
\boxed{G=c=1,\qquad 4\pi\epsilon_0=1}
\]

と固定する。

`ψ` に保存された `G,c` はともに 1 でなければならない。異なる値を silent conversion してはならない。

4D volume form orientation は interaction v2.0 と同じく

\[
\epsilon_{0123}=+\sqrt{-g},
\qquad
\epsilon^{0123}=-\frac1{\sqrt{-g}}.
\]

spacetime index：`a,b,c=0..3`  
spatial index：`i,j,k=1..3`  
black-hole label：`A=1..N_BH`。

---

# 3. 初期化器が受け取る state

## 3.1 ψ_phys0

最低限：

\[
N_{BH},
M_A^{(0)},Q_A^{(0)},\mathbf J_A^{(0)},
\mathbf X_A^{(0)},
P_0,E_0,\phi_0,\chi_0,\tau_0,
q_0,\nu_0,n_q,s_A,
\lambda_A^{(0)},C_0,D_0.
\]

`N_BH` は 0,1,2 を本仕様で許可する。

- `N_BH=0`: Minkowski validation
- `N_BH=1`: single Kerr–Newman / RN validation
- `N_BH=2`: binary experiment

これらは**同じ sum-over-A 形式の initializer**で扱う。

## 3.2 ψ_initcfg

初期化の再現性に影響する全設定を恒等 state として保持する。

必須：

```text
initializer_method_id = SKN_XCTS_EM_V2
zeta_att
Omega0
dot_a0
v0_x,v0_y,v0_z
Omega_r_A_x,Omega_r_A_y,Omega_r_A_z
epsilon_init
epsilon_linear
epsilon_horizon_match
epsilon_mass_match
epsilon_charge_match
epsilon_spin_match
epsilon_ADM_P
I_Newton_max
I_GMRES_max
I_mass_match_max
I_line_search_max
armijo_c
line_search_factor
jacobian_method_id = FORWARD_AD
linear_solver_id = GMRES_BLOCK_JACOBI
```

baseline：

```text
zeta_att = 2
Omega0 = Omega_guess
dot_a0 = 0
v0 = (0,0,0)
armijo_c = 1e-4
line_search_factor = 0.5
```

`zeta_att=2` は charged-binary initial data で使用実績のある attenuation choice を本研究の共通 baseline として固定する。別値は別 run の明示 state 変更としてのみ許可する。

## 3.3 ψ_grid

interaction v2.0 の state contract をそのまま要求する。

最低限：

\[
\mathbf D_1,\mathbf D_2,\mathbf D_3,
\mathbf V_*,\mathbf P_*,
\]

及び domain/map/basis/quadrature coefficients。

受理条件：

\[
\|\mathbf P_*\mathbf V_*-I\|_\infty<\epsilon_{grid}
\]

かつ `D_i` が physical-coordinate derivative operator であること。

initializer は `ψ_grid` 以外の hidden grid、hidden finite-difference stencil、hidden ghost value を参照してはならない。

---

# 4. 初期物理入力の admissibility gate

各 black hole A について

\[
\chi_A^{(0)}=\frac{|\mathbf J_A^{(0)}|}{M_A^{(0)2}},
\qquad
\mathcal Q_A^{(0)}=\frac{Q_A^{(0)}}{M_A^{(0)}}.
\]

Kerr–Newman regularity：

\[
\boxed{
\chi_A^{(0)2}+\mathcal Q_A^{(0)2}\le1
}
\]

を必須とする。

違反時：

```text
INITIALIZATION_REJECTED_SUPEREXTREMAL
```

で終了。

等号 `=1` は許可する。extremal case を別解へ置換しない。

また

\[
M_A^{(0)}>0
\]

を必須とする。

---

# 5. baseline binary coordinate data

旧 `P_0` は v1.0 と同じく dimensionless separation scale とする。

\[
\boxed{d_0=R_0=P_0 M_{tot}^{(0)}}
\]

\[
M_{tot}^{(0)}=\sum_A M_A^{(0)}.
\]

binary の center of mass：

\[
\mathbf r_{CM}
=\frac{\sum_A M_A^{(0)}\mathbf X_A^{(0)}}{M_{tot}^{(0)}}.
\]

baseline equal-mass では

\[
\mathbf X_1^{(0)}=(+d_0/2,0,0),
\quad
\mathbf X_2^{(0)}=(-d_0/2,0,0).
\]

旧係数

\[
C_0=1-\lambda_1^{(0)}\lambda_2^{(0)}
\]

から quasi-circular initial angular-velocity guess：

\[
\boxed{
\Omega_{guess}^2
=\frac{M_{tot}^{(0)}C_0}{d_0^3}
}
\]

とする。

`Omega0` state の baseline 値は `Omega_guess`。

\[
\boxed{\dot a_0=0}
\]

を baseline とする。

個々の superposed Kerr–Newman seed の boost velocity は

\[
\boxed{
\mathbf v_A
=\mathbf\Omega_0\times(\mathbf X_A^{(0)}-\mathbf r_{CM})
+\dot a_0(\mathbf X_A^{(0)}-\mathbf r_{CM})
}
\]

で固定する。

ここで `Ω0=(0,0,Omega0)`。

`P_A^(0)` は dynamics coefficient ではなく、この velocity から

\[
\mathbf P_A^{(0)}
=\gamma_A M_A^{(0)}\mathbf v_A,
\qquad
\gamma_A=(1-|\mathbf v_A|^2)^{-1/2}
\]

として保存する audit identity とする。

`|v_A|>=1` は initialization FAIL。

---

# 6. attenuation width の唯一の定義

各 black hole の Gaussian attenuation width は

\[
\boxed{
w_A
=\zeta_{att}\frac{M_A^{seed}}{\sum_B M_B^{seed}}d_0
}
\]

とする。

weight：

\[
\boxed{
W_A(\mathbf x)=\exp[-r_A^2/w_A^2]
}
\]

\[
r_A=|\mathbf x-\mathbf X_A^{(0)}|.
\]

metric、K、electric background、magnetic vector potential で**同じ W_A**を使う。

---

# 7. individual Kerr–Newman seed — 唯一の式

各 A の seed parameter：

\[
M_A^{seed},\quad Q_A^{seed}=Q_A^{(0)},
\quad
\mathbf a_A=\frac{\mathbf J_A^{seed}}{M_A^{seed}}.
\]

initial outer-loop guess：

\[
M_A^{seed}=M_A^{(0)},
\qquad
\mathbf J_A^{seed}
=\boldsymbol\chi_A^{(0)} M_A^{seed2}.
\]

各 hole の rest-frame centered coordinate を

\[
\mathbf x_A=\mathbf x-\mathbf X_A^{(0)}
\]

とする。

\[
\rho_A^2=\mathbf x_A\cdot\mathbf x_A.
\]

Kerr–Newman radial coordinate：

\[
\boxed{
r_A^2
=\frac12(\rho_A^2-a_A^2)
+\sqrt{
\frac14(\rho_A^2-a_A^2)^2
+(\mathbf a_A\cdot\mathbf x_A)^2
}
}
\]

metric function：

\[
\boxed{
H_A
=\frac{M_A^{seed}r_A^3-\frac12(Q_A^{seed}r_A)^2}
{r_A^4+(\mathbf a_A\cdot\mathbf x_A)^2}
}
\]

null vector rest-frame components：

\[
l_0=1,
\]

\[
\boxed{
\mathbf l
=\frac{r_A[\mathbf x_A-(\hat{\mathbf a}_A\cdot\mathbf x_A)\hat{\mathbf a}_A]
-\mathbf a_A\times\mathbf x_A}{r_A^2+a_A^2}
+\frac{(\hat{\mathbf a}_A\cdot\mathbf x_A)\hat{\mathbf a}_A}{r_A}
}
\]

`a_A=0` では連続極限を使用し、任意方向 `ahat` を参照してはならない。

4-metric：

\[
\boxed{
g_{\mu\nu}^{A}=\eta_{\mu\nu}+2H_A l_\mu l_\nu
}
\]

3+1 fields：

\[
\gamma_{ij}^A=\delta_{ij}+2H_A l_i l_j,
\]

\[
\alpha_A=(1+2H_A l_0^2)^{-1/2},
\]

\[
\beta_i^A=2H_A l_0l_i.
\]

Kerr–Newman 4-potential：

\[
\boxed{
A_\mu^A
=-\frac{Q_A^{seed}r_A^3}
{r_A^4+(\mathbf a_A\cdot\mathbf x_A)^2}
l_\mu
}
\]

とする。

\[
F_{\mu\nu}^A=\partial_\mu A_\nu^A-\partial_\nu A_\mu^A.
\]

interaction v2.0 と同じ normal / dual convention で

\[
E_A^\mu=-n_\nu F_A^{\nu\mu},
\qquad
B_A^\mu=-n_\nu {}^*F_A^{\nu\mu}
\]

を作る。

---

# 8. boost の唯一の定義

各 A を constant Lorentz boost `v_A` で inertial coordinates へ写す。

Lorentz matrix：

\[
\Lambda^0{}_0=\gamma,
\quad
\Lambda^0{}_i=\gamma v_i,
\quad
\Lambda^i{}_0=\gamma v^i,
\]

\[
\Lambda^i{}_j
=\delta^i{}_j+\frac{\gamma-1}{v^2}v^iv_j.
\]

`v=0` は identity limit。

rest-frame coordinate は

\[
\bar x^\alpha=(\Lambda^{-1})^\alpha{}_\beta x^\beta
\]

で得る。

scalar `H` は composition、vector/covector は Lorentz tensor transformation を用いる。

**boosted seed の空間・時間微分は finite difference で近似してはならない。**

```text
seed_derivative_method_id = FORWARD_AUTODIFF_ANALYTIC_KN
```

を固定し、analytic Kerr–Newman formula + Lorentz map 全体を forward automatic differentiation して

\[
\partial_\alpha g^A_{\mu\nu},
\quad
\partial_\alpha A^A_\mu
\]

を生成する。

これより

\[
K^A_{ij}
=\frac1{2\alpha_A}
(D_i\beta^A_j+D_j\beta^A_i-\partial_t\gamma^A_{ij})
\]

を得る。

---

# 9. excision surface seed

各 Kerr–Newman seed の outer horizon radius：

\[
\boxed{
r_{+,A}
=M_A^{seed}
\left[1+\sqrt{1-\chi_A^2-\mathcal Q_A^2}\right]
}
\]

\[
\chi_A=\frac{|\mathbf J_A^{seed}|}{M_A^{seed2}},
\qquad
\mathcal Q_A=\frac{Q_A^{seed}}{M_A^{seed}}.
\]

excision surface `S_A` は rest-frame Kerr–Newman horizon を Lorentz boost した surface とする。

`ψ_grid` の inner boundary map はこの `S_A` と一致していなければならない。

一致誤差が `epsilon_surface_match` を超える場合

```text
GRID_EXCISION_MISMATCH
```

で FAIL。

**initializer は grid shape を勝手に変更しない。**

---

# 10. XCTS free data

source-derived constructionを固定する。

conformal metric：

\[
\boxed{
\tilde\gamma_{ij}
=\delta_{ij}
+\sum_{A=1}^{N_{BH}}
W_A(\gamma^A_{ij}-\delta_{ij})
}
\]

extrinsic-curvature trace free datum：

\[
\boxed{
K
=\sum_{A=1}^{N_{BH}}W_A K_A
}
\]

where

\[
K_A=\gamma_A^{ij}K^A_{ij}.
\]

free time derivatives：

\[
\boxed{\tilde u_{ij}=0}
\]

\[
\boxed{\partial_tK=0}
\]

と固定する。

---

# 11. XCTS unknowns

唯一の gravitational elliptic unknowns：

\[
\boxed{
\psi,\quad X:=\alpha\psi,\quad \beta^i
}
\]

とする。

\[
\bar\alpha=\alpha\psi^{-6}=X\psi^{-7}.
\]

conformal longitudinal operator：

\[
(\tilde L\beta)^{ij}
=\tilde D^i\beta^j+\tilde D^j\beta^i
-\frac23\tilde\gamma^{ij}\tilde D_k\beta^k.
\]

\[
\boxed{
\tilde A^{ij}
=\frac{\psi^7}{2X}
[(\tilde L\beta)^{ij}-\tilde u^{ij}]
}
\]

とする。

physical geometry：

\[
\gamma_{ij}=\psi^4\tilde\gamma_{ij}
\]

\[
A^{ij}=\psi^{-10}\tilde A^{ij}
\]

\[
K_{ij}=A_{ij}+\frac13\gamma_{ij}K.
\]

---

# 12. electromagnetic unknown and fields

唯一の EM elliptic unknown：

\[
\boxed{\varphi}
\]

とする。

background electric field：

\[
\boxed{
(E_{sp})^i
=\sum_A W_A E_A^i
}
\]

background spatial vector potential：

\[
\boxed{
(A_{sp})_i
=\sum_A W_A (A_A)_i
}
\]

solved electric field：

\[
\boxed{
E^i=(E_{sp})^i+D^i\varphi
}
\]

magnetic field：

\[
\boxed{
B^i
=\frac1{\sqrt\gamma}[ikl]\partial_k(A_{sp})_l
}
\]

とする。

`B` は独立 elliptic unknown にしない。

全 derivatives は current physical metric と approved `D_i` から構成する。

---

# 13. electromagnetic source terms for XCTS

\[
\rho
=\frac{E_iE^i+B_iB^i}{8\pi}
\]

\[
J^i
=\frac1{4\pi}(E\times B)^i
\]

\[
(E\times B)^i
=\frac1{\sqrt\gamma}[ikl]E_kB_l
\]

Maxwell では stress trace spatial quantity

\[
\boxed{\sigma=\rho}
\]

を使用する。

これらは Newton residual evaluation ごとに current unknowns `ψ,X,β,φ` から再計算する。

前 Newton iteration の source を frozen してはならない。

---

# 14. coupled XCTS + EM equations

未知量

\[
U=(\psi,X,\beta^1,\beta^2,\beta^3,\varphi)
\]

を**同時に**解く。

Hamiltonian：

\[
\boxed{
\tilde D^2\psi
-\frac18\tilde R\psi
-\frac1{12}K^2\psi^5
+\frac18\psi^{-7}\tilde A_{ij}\tilde A^{ij}
=-2\pi\psi^5\rho
}
\]

Momentum：

\[
\boxed{
\tilde D_j\left[
\frac{\psi^7}{2X}(\tilde L\beta)^{ij}
\right]
-\frac23\psi^6\tilde D^iK
-\tilde D_j\left[
\frac{\psi^7}{2X}\tilde u^{ij}
\right]
=8\pi\psi^{10}J^i
}
\]

Lapse：

\[
\boxed{
\tilde D^2X
-X\left[
\frac18\tilde R
+\frac5{12}K^2\psi^4
+\frac78\psi^{-8}\tilde A_{ij}\tilde A^{ij}
\right]
+\psi^5(\partial_tK-\beta^k\partial_kK)
=2\pi X\psi^4(\rho+2\sigma)
}
\]

Electric Gauss correction：

\[
\boxed{
D_iD^i\varphi
=-D_i(E_{sp})^i
}
\]

である。

computational domain は electrovac なので volume charge density は 0。

---

# 15. outer boundary conditions

spatial infinity `i0`：

\[
\boxed{\psi=1}
\]

\[
\boxed{X=\alpha\psi=1}
\]

\[
\boxed{\beta^i=v_0^i}
\]

\[
\boxed{\varphi=0}
\]

とする。

corotation `Omega0 x r` と expansion `dot_a0 r` は outer boundary へ入れない。

---

# 16. inner boundary — shift

各 excision surface `S_A` で

\[
\boxed{
\beta^i
=\alpha s^i
-\Omega_{r,A}^k\xi^i_{(k),A}
-[\mathbf\Omega_0\times(\mathbf r-\mathbf r_{CM})]^i
-\dot a_0(\mathbf r-\mathbf r_{CM})^i
}
\]

とする。

flat rotational vector：

\[
\boxed{
\xi^i_{(k),A}
=[ikl](r-X_A)_l
}
\]

を使用する。

---

# 17. inner boundary — conformal factor and lapse

conformal outward normal `tilde s^i`、surface conformal metric `tilde h_ij` を用い、

\[
\boxed{
\tilde s^k\partial_k\psi
=-\frac{\psi^3}{8\alpha}
\tilde s^i\tilde s^j(\tilde L\beta)_{ij}
-\frac\psi4\tilde h^{ij}\tilde D_i\tilde s_j
+\frac16K\psi^3
}
\]

を imposed apparent-horizon condition とする。

lapse gauge condition：

\[
\boxed{
X
=1+\sum_A W_A(\alpha_A-1)
\quad\text{on }S_A
}
\]

とする。

---

# 18. inner boundary — electric potential and charge control

各 `S_A` の desired charge：

\[
Q_{d,A}=Q_A^{(0)}.
\]

superposed-field charge：

\[
\boxed{
Q_{sp,A}
=\frac1{4\pi}\oint_{S_A}(E_{sp})^is_i\,dA
}
\]

とする。

## nonzero target charge

\[
|Q_{d,A}|>\epsilon_Q
\]

なら、

\[
|Q_{sp,A}|\le\epsilon_{Qsp}
\]

で FAIL：

```text
INITIALIZATION_FAIL_ZERO_QSP
```

それ以外は

\[
\boxed{
s^i\partial_i\varphi
=\left(\frac{Q_{d,A}}{Q_{sp,A}}-1\right)
s^i(E_{sp})_i
}
\]

とする。

## zero target charge

\[
|Q_{d,A}|\le\epsilon_Q
\]

なら 0/0 を避けるため、その連続的な zero-flux condition として

\[
\boxed{
s^i\partial_i\varphi
=-s^i(E_{sp})_i
}
\]

を用いる。

これにより final normal electric flux を 0 target にする。

---

# 19. nonlinear solve method — 唯一の実装

XCTS+EM residual を

\[
\mathcal R(U)=0
\]

とする。

唯一の nonlinear solver：

```text
NEWTON_RAPHSON
```

Jacobian action：

```text
FORWARD_AD
```

で residual code を自動微分して生成する。finite-difference Jacobianは禁止。

Newton linear system：

\[
J(U_n)\,\delta U_n=-\mathcal R(U_n)
\]

唯一の linear solver：

```text
GMRES_BLOCK_JACOBI
```

とする。

line search：Armijo backtracking。

\[
U_{n+1}=U_n+\lambda\delta U_n
\]

`lambda=1` から開始し、

\[
\|R(U_n+\lambda\delta U_n)\|_2
\le(1-c\lambda)\|R(U_n)\|_2
\]

を満たすまで

\[
\lambda\leftarrow0.5\lambda
\]

とする。

`I_line_search_max` を超えたら FAIL。別solverへfallback禁止。

---

# 20. Newton initial guess

一意に

\[
\boxed{\psi^{[0]}=1}
\]

\[
\boxed{
X^{[0]}
=1+\sum_AW_A(\alpha_A-1)
}
\]

\[
\boxed{
\beta^{i[0]}
=v_0^i+\sum_AW_A\beta_A^i
}
\]

\[
\boxed{\varphi^{[0]}=0}
\]

とする。

別 initial guess へ自動切替しない。

---

# 21. Newton iteration 内の v0 control

ADM linear momentum を current solution iterate から計算し、source construction に従い

\[
\boxed{
\mathbf v_{0,n+1}
=\mathbf v_{0,n}
-\frac{\mathbf P_{ADM,n}}{M_T}
}
\]

とする。

初期値：

\[
\boxed{\mathbf v_{0,0}=0}
\]

`M_T` は current sum of target/seed horizon masses とし、iteration中に定義を変えない。

---

# 22. Newton iteration 内の spin control

各 hole の rotational shift parameter 初期値：

\[
\boxed{
\mathbf\Omega_{r,A,0}
=\frac{\mathbf J_{KN,A}}
{4M_{Chr,KN,A}M_{irr,KN,A}^2}
}
\]

とする。

current horizon angular momentum `J_A,n` を測り、

\[
\boxed{
\mathbf\Omega_{r,A,n+1}
=\mathbf\Omega_{r,A,n}
+\frac{M_{Chr,A,n}^2\boldsymbol\chi_{KN,A}-\mathbf J_{A,n}}
{4M_{Chr,A,n}M_{irr,A,n}^2}
}
\]

と更新する。

spinless target はこの同じ式の `chi=0` 極限。

---

# 23. initialization-target horizon readouts

initial matching のためだけに以下を用いる。これらは dynamics state ではない。

charge：

\[
\boxed{
Q_{H,A}
=\frac1{4\pi}\oint_{S_A}E^is_i\,dA
}
\]

irreducible mass：

\[
\boxed{
M_{irr,A}=\sqrt{\frac{A_A}{16\pi}}
}
\]

charged horizon angular momentum は EM contribution を含め、

\[
\boxed{
J_{A,k}
=\frac1{8\pi}
\oint_{S_A}
\left[K_{ij}+2(A_{sp})_iE_j\right]
\xi^i_{(k),A}s^j\,dA
}
\]

とする。

Christodoulou–Ruffini mass：

\[
\boxed{
M_{H,A}^2
=\left(M_{irr,A}+\frac{Q_{H,A}^2}{4M_{irr,A}}\right)^2
+\frac{J_A^2}{4M_{irr,A}^2}
}
\]

を用いる。

---

# 24. seed mass outer matching loop

旧 v1.0 の「seedを反復補正する」を次で一意化する。

outer control variable：

\[
\mathbf m^{seed}=(M_1^{seed},...,M_{N_{BH}}^{seed}).
\]

target residual：

\[
r_A^M=M_{H,A}-M_A^{(0)}.
\]

初期値：

\[
M_A^{seed,(0)}=M_A^{(0)}.
\]

各 outer iteration で full coupled XCTS+EM solve を完了してから `r^M` を評価する。

Jacobian は symmetric central perturbation full-solve により

\[
J^M_{AB}
=\frac{r_A^M(M_B^{seed}+\delta_B)-r_A^M(M_B^{seed}-\delta_B)}{2\delta_B}
\]

\[
\delta_B
=\sqrt{\epsilon_{mach}}\max(1,|M_B^{seed}|)
\]

で作る。

\[
J^M\Delta\mathbf m=-\mathbf r^M
\]

を解き、Armijo line search を同じ `c=1e-4, factor=0.5` で適用する。

seed mass更新時、dimensionless target spin `chi_A^(0)` は固定し、

\[
\mathbf J_A^{seed}
=\boldsymbol\chi_A^{(0)}M_A^{seed2}
\]

を再構成する。

charge target `Q_A^(0)` は変更しない。

`I_mass_match_max` で収束しなければ FAIL。別 mass-renormalization 法へfallback禁止。

---

# 25. coupled solve convergence conditions

Newton内部収束：

\[
\|\mathcal R_{XCTS}\|_\infty<\epsilon_{init}
\]

\[
\|\mathcal R_{\varphi}\|_\infty<\epsilon_{init}
\]

\[
\|\mathcal R_{BC}\|_\infty<\epsilon_{init}
\]

を同時に満たす。

physical constraint re-evaluation：

\[
\boxed{
\|\mathcal H\|_\infty<\epsilon_{init}
}
\]

\[
\boxed{
\|\mathcal M_i\|_\infty<\epsilon_{init}
}
\]

\[
\boxed{
\|D_iE^i\|_\infty<\epsilon_{init}
}
\]

\[
\boxed{
\|D_iB^i\|_\infty<\epsilon_{init}
}
\]

を solver residualとは独立に再計算する。

horizon target：

\[
|M_{H,A}-M_A^{(0)}|<\epsilon_{mass-match}
\]

\[
|Q_{H,A}-Q_A^{(0)}|<\epsilon_{charge-match}
\]

\[
|\mathbf J_{H,A}-\mathbf J_A^{(0)}|<\epsilon_{spin-match}
\]

ADM momentum：

\[
|\mathbf P_{ADM}|<\epsilon_{ADM-P}.
\]

全条件を満たすまで `ψ₀` を生成したことにしてはならない。

---

# 26. physical ADM fields after solve

\[
\gamma_{ij}=\psi^4\tilde\gamma_{ij}
\]

\[
\alpha=X/\psi
\]

\[
A^{ij}=\psi^{-10}\tilde A^{ij}
\]

\[
K_{ij}=A_{ij}+\frac13\gamma_{ij}K.
\]

これと solved `E^i,B^i` が ADM initial data の唯一の正本。

superposed seed fieldsを dynamicsへ直接渡してはならない。

---

# 27. ADM -> spacetime metric

interaction v2.0 と一致する

\[
\boxed{
g_{ij}=\gamma_{ij}}
\]

\[
\boxed{
g_{0i}=\gamma_{ij}\beta^j}
\]

\[
\boxed{
g_{00}=-\alpha^2+\gamma_{ij}\beta^i\beta^j}
\]

を構成する。

---

# 28. Phi initialization

interaction と同じ derivative operator を使用し、

\[
\boxed{
\Phi_{iab}=\mathbf D_i g_{ab}
}
\]

とする。

別 analytic derivative、別 FD derivativeを使用しない。

これにより initial three-index GH constraint は semi-discrete level で

\[
\mathcal C_{iab}=D_i g_{ab}-\Phi_{iab}=0
\]

となる。

---

# 29. initial gauge source H

interaction v2.0 の唯一の damped-wave target

\[
F_a
=\mu_L\log(\gamma^{p_H}/\alpha)n_a
-\mu_S\alpha^{-1}g_{ai}\beta^i
\]

を ADM solved metric から計算する。

初期 gauge state：

\[
\boxed{H_a^{(0)}=F_a^{(0)}}
\]

と固定する。

これは initial gauge-driver fixed target choice である。

---

# 30. Pi_ij の一意な初期化

interaction v2.0 の identity

\[
K_{ij}
=\frac12\Pi_{ij}
+\frac12n^a(\Phi_{ija}+\Phi_{jia})
\]

を逆に解き、

\[
\boxed{
\Pi_{ij}
=2K_{ij}-n^a(\Phi_{ija}+\Phi_{jia})
}
\]

とする。

これにより `K_ij` と GH variables の符号を一意に一致させる。

---

# 31. Pi_00, Pi_0i の一意な初期化

残る4成分を lapse/shift time derivative の任意選択で決めてはならない。

未知 vector：

\[
p_A=(\Pi_{00},\Pi_{01},\Pi_{02},\Pi_{03}).
\]

現在の既知 `g,Phi,Pi_ij` と trial `p_A` から interaction v2.0 の式

\[
\partial_0g_{ab}
=\beta^i\Phi_{iab}-\alpha\Pi_{ab}
\]

を用いて 4D Christoffel と `Gamma_a(p)` を構成する。

初期 GH gauge constraint

\[
\boxed{
\mathcal C_a
=H_a+\Gamma_a=0
}
\]

を4本の線形方程式として

\[
\boxed{
\Gamma_a(p)=-H_a
}
\]

と解き、`Pi_00,Pi_0i` を一意に決定する。

4x4 coefficient matrix が singular / ill-conditioned (`cond>A_cond_max`) なら FAIL。

これにより hidden `partial_t alpha, partial_t beta` choice を完全に排除する。

---

# 32. Theta initialization

interaction v2.0 gauge-driver

\[
\partial_tH_a
=\beta^kD_kH_a-\mu_H(H_a-F_a)+\Theta_a
\]

で `H=F` を採用したため、initial `partial_t H=0` を固定するには

\[
\boxed{
\Theta_a^{(0)}=-\beta^kD_kH_a^{(0)}
}
\]

とする。

この選択は同時に second gauge-driver equation の initial fixed-point condition と整合する。

---

# 33. interaction state への projection

全 field は interaction と同じ evaluation/projection pair を用い、

\[
\boxed{
u_{coeff}=\mathbf P_*u_*}
\]

で retained state coefficients へ格納する。

格納対象：

- `g_ab`
- `Pi_ab`
- `Phi_iab`
- `H_a`
- `Theta_a`
- `E^i`
- `B^i`

である。

retained coefficient上で独自 interpolation を行わない。

---

# 34. time/unit/identity state initialization

\[
\boxed{ONE=1}
\]

\[
\boxed{T_\tau=e^{\kappa_\tau\tau_0}}
\]

\[
\boxed{D_\tau=e^{\kappa_\tau\Delta\tau}}
\]

を設定する。

`kappa_tau=0` は interaction v2.0 で不許可。

`D_tau>0`, `T_tau>0` を必須とする。

全 `ψ_phys0`,`ψ_initcfg`,`ψ_grid`,`ψ_numerics` identity state を同じ `ψ₀` に格納する。

---

# 35. interaction compatibility audit

`ψ₀` を返す前に interaction v2.0 の derived-geometry routinesだけを使って独立監査する。

必須：

1. `ONE=1`
2. `c=1`
3. metric signature `(-,+,+,+)`
4. `gamma>0`
5. `alpha>0`
6. `P_* V_* = I`
7. `C_iab=0` within `epsilon_GH_init`
8. `C_a=0` within `epsilon_GH_init`
9. interaction formula から再構成した `K_ij` が XCTS solved `K_ij` と一致
10. interaction Maxwell `F_ab` / `T_ab` reconstruction と initializer `rho,J,sigma` が一致
11. Maxwell trace audit
12. all identity state finite and fixed

不一致なら `ψ₀` を出力せず FAIL。

---

# 36. initialization fixed-point / limiting-case audits

同じ initializer の極限として以下を必須 validation にする。

## N_BH=0

flat free data, zero EM で

\[
\psi=1,\alpha=1,\beta=0,E=B=0
\]

を回復。

## N_BH=1, Q=0, J=0

single Schwarzschild/Kerr-Schild-XCTS limit。

## N_BH=1, Q!=0, J=0

single Reissner–Nordstrom / Kerr–Newman `a=0` limit。

## N_BH=1, generic admissible Q,J

single Kerr–Newman limit。

これらを別 initializer で作って比較してはならない。

---

# 37. extremal case rule

\[
\chi_A^2+\mathcal Q_A^2=1
\]

でも同じ SKN-XCTS-EM initializer を使用する。

許可されない行為：

- MP exact solutionへの置換
- `Q/M=0.999999` への変更
- excision radius の任意拡大
- mass/charge renormalization
- alternative puncture initializer fallback

数値的に収束しなければ

```text
INITIALIZATION_FAILED_EXTREMAL
```

として記録する。

---

# 38. superextremal rule

\[
\chi_A^2+\mathcal Q_A^2>1
\]

は regular black-hole initializer の domain 外。

```text
INITIALIZATION_REJECTED_SUPEREXTREMAL
```

を返し、state生成を開始しない。

---

# 39. 禁止事項

1. case-specific initializer switch
2. MP substitution
3. double-RN exact field substitution
4. point-particle initial field
5. unconstrained superposed KNをそのままψ0へ採用
6. solver failure時のFD/RBF/puncture法へのfallback
7. hidden derivative operator
8. hidden solver tolerance
9. `Q,M,J` target の無断変更
10. `zeta_att` のrun途中変更
11. Newton initial guess の自動変更
12. charge boundary条件の無断変更
13. `Pi_0a` を arbitrary `dt alpha,dt beta` で決めること
14. interactionと異なる Maxwell sign convention
15. interactionと異なる projection/dealiasing

---

# 40. FAIL code

最低限：

```text
INITIALIZATION_REJECTED_SUPEREXTREMAL
GRID_CONTRACT_FAIL
GRID_EXCISION_MISMATCH
INITIALIZATION_FAIL_ZERO_QSP
XCTS_NEWTON_DIVERGED
XCTS_LINEAR_SOLVE_FAILED
XCTS_LINE_SEARCH_FAILED
XCTS_CONSTRAINT_FAIL
HORIZON_MASS_MATCH_FAIL
HORIZON_CHARGE_MATCH_FAIL
HORIZON_SPIN_MATCH_FAIL
ADM_MOMENTUM_MATCH_FAIL
GH_PI_GAUGE_SOLVE_FAIL
GH_CONSTRAINT_INIT_FAIL
INTERACTION_COMPATIBILITY_FAIL
INITIALIZATION_FAILED_EXTREMAL
```

FAIL時に別ルールへ切り替えない。

---

# 41. ψ₀ output contract

成功時に返すもの：

\[
\boxed{
\psi_0
=(\psi_{unit},\psi_{time},\psi_{phys0},\psi_{numerics},\psi_{grid},\psi_{field})
}
\]

と、監査ログ：

```text
INITIALIZATION_INPUT.json
INITIALIZATION_METHOD.json
KN_SEED_PARAMETERS.json
XCTS_FREE_DATA.json
XCTS_NEWTON_HISTORY.csv
XCTS_LINEAR_SOLVE_HISTORY.csv
MASS_MATCH_HISTORY.csv
HORIZON_TARGET_HISTORY.csv
CONSTRAINT_AUDIT.json
GH_CONVERSION_AUDIT.json
INTERACTION_COMPATIBILITY_AUDIT.json
INITIAL_STATE.npy
STATE_SCHEMA.json
RUN_HASHES.sha256
```

を保存する。

---

# 42. v1.0 からの修正点

v1.0 では「superposed Kerr–Newman + XCTS + electromagnetic constraints」と方針だけを書き、以下が未定義だった。

v2.0 では次を固定した。

1. Kerr–Newman seed の具体式
2. boost rule
3. attenuation weight と baseline zeta
4. XCTS free data
5. XCTS unknowns
6. XCTS 3 elliptic equations
7. electric-potential elliptic equation
8. magnetic field construction
9. outer/inner geometry BC
10. charge-fixing EM BC
11. zero-charge 0/0 limit rule
12. Newton/GMRES/AD/line-search solver
13. Newton initial guess
14. ADM momentum control
15. spin control
16. mass-seed outer matching loop
17. full physical convergence gates
18. ADM→GH conversion
19. Pi_ij construction
20. Pi_0a gauge-constraint linear solve
21. H=F initialization
22. Theta initialization
23. same D,V,P use
24. interaction-v2 compatibility audit
25. extremal case no-substitution rule

---

# 43. source-derived部分と本研究固定選択の区別

## source-derived

Mukherjee et al. (2022) に基づく：

- superposed boosted Kerr–Newman conformal metric / K
- Gaussian attenuation
- `u_tilde=0`, `dt K=0`
- XCTS equations
- outer / inner boundary conditions
- electric scalar potential correction
- magnetic vector-potential curl construction
- charge-fixing Neumann BC
- Newton-Raphson + GMRES + block-Jacobi
- v0 ADM-momentum iteration
- Omega_r spin iteration
- attenuation relation `w_A=zeta M_A/M_T d`
- Kerr–Newman Kerr-Schild formula / boost formula

## project-fixed to close v1.0 gaps

本仕様 v2.0 で追加固定したもの：

- `zeta_att=2` を共通baseline
- Newton Jacobian actionを FORWARD_AD に固定
- Armijo parameters
- mass-seed outer Newton matching
- zero-charge phi Neumann limit
- H=F initial gauge
- Pi_ij from interaction K identity
- Pi_0a from `C_a=0` 4x4 solve
- Theta=-beta·grad H
- interactionと同じ D,V,P の使用
- extremalでMPへ切替禁止

これらはsourceに記載された唯一の方法であると主張せず、**本研究で実装自由度を消すために固定した仕様選択**である。

---

# 44. 参照標準

1. S. Mukherjee, N. K. Johnson-McDaniel, W. Tichy, S. L. Liebling, “Conformally curved initial data for charged, spinning black hole binaries on arbitrary orbits,” arXiv:2202.12133 / Phys. Rev. D.
2. G. Lovelace, R. Owen, H. P. Pfeiffer, T. Chu, superposed Kerr-Schild XCTS binary-black-hole initial data, Phys. Rev. D 78, 084017 (2008).
3. G. B. Cook, H. P. Pfeiffer, “Excision boundary conditions for black-hole initial data,” Phys. Rev. D 70, 104016 (2004).
4. J. W. York Jr.; H. P. Pfeiffer & J. W. York Jr., XCTS formulation.
5. 上位仕様 `universal_interaction_operator_complete_spec_ja_v2.0_20260930.md`。

---

# 45. 実装開始条件

初期化コードを実装してよいのは以下を全て満たした後のみ。

- [ ] interaction v2.0 が正本として承認済み
- [ ] 本 initialization v2.0 が承認済み
- [ ] `ψ_grid` / `D_i,V_*,P_*` が interaction state contract を満たす
- [ ] excision surface と grid inner boundary が一致
- [ ] all init/numerical parameters が state schema に登録済み
- [ ] no hidden fallback path
- [ ] no case-specific initializer branch other than deterministic admissibility/zero-charge boundary handling

一つでも未完成なら正式実験コードを走らせない。

---

# 46. 最終定義

v2.0 初期化写像は

\[
\boxed{
\mathcal I:
\{M_A,Q_A,J_A,X_A,P_0,\text{all init states},\psi_{grid}\}
\xrightarrow{\text{boosted KN free data}}
\{\tilde\gamma,K,\tilde u,\partial_tK,E_{sp},A_{sp}\}
}
\]

\[
\boxed{
\xrightarrow{\text{coupled XCTS + EM}}
\{\gamma,K_{ij},\alpha,\beta,E,B\}
\xrightarrow{\text{interaction-compatible ADM→GH}}
\psi_0
}
\]

である。

**この写像以外の初期化経路で生成された state は、interaction v2.0 が動作しても正式実験の `ψ₀` として認めない。**
