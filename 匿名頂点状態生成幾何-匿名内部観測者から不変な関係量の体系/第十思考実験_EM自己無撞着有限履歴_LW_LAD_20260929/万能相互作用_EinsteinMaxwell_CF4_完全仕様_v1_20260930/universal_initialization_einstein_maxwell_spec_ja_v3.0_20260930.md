# 万能相互作用 Einstein–Maxwell 初期化仕様 v3.0
## — 乗法正式系・加法検算系へ同一の完全初期状態 Ψ₀ を供給する共通初期化仕様 —

**著者:** 木原範昭 (Noriaki Kihara)  
**仕様整理:** ChatGPT  
**版:** v3.0  
**日付:** 2026-09-30  
**上位相互作用仕様:** `universal_interaction_einstein_maxwell_multiplicative_spec_ja_v3.0_20260930.md`  
**置換対象:** `universal_initialization_einstein_maxwell_xcts_spec_ja_v2.0_20260930.md`

---

# 0. 改定理由

v2.0 には以下の問題があった。

1. 上位仕様が interaction v2.0 / CF4 / K-block 前提のままだった。
2. `N_BH=0,1,2` を同じ initializer としながら、binary 専用量 `d0,C0,Omega_guess` を 0/1 体にも適用する自己矛盾があった。
3. XCTS Newton 未知量を `U=(psi,X,beta^i,varphi)` と定義しながら、Newton iteration 中に `v0` と `Omega_r,A` を変更し、解く方程式自体を iteration 中に変えていた。
4. seed mass を outer loop で変更すると seed Kerr–Newman horizon が変わる一方、grid inner boundary を seed horizon に一致させて固定するという矛盾があった。
5. interaction v2.0 固有の `D_i,V_*,P_*` を initializer にまで要求していたが、interaction v3.0 では物理 RHS と数値 solver の既存離散化をそのまま用いる方針へ改定された。
6. `ONE,T_tau,delta_tau,D_tau,kappa_tau` の v3.0 time-state contract を初期化側が満たしていなかった。
7. 乗法系用・加法系用に別初期化を作らないという実験上の最重要条件が明文化されていなかった。

v3.0 はこれらを修正する。

---

# 1. 最上位原則

初期化器の唯一の仕事は、正式乗法系と検算加法系の双方へ渡す **一つだけの完全状態**

\[
\boxed{\Psi_0}
\]

を生成することである。

初期化写像を

\[
\boxed{
\mathcal I:
(\Psi_{phys0},\Psi_{initcfg},\Psi_{numerics},\Psi_{grid})
\longmapsto \Psi_0
}
\]

とする。

成功後、実験 harness はこの一つの `Psi_0` を bitwise copy して

\[
\boxed{
\Psi_0^{(\times)}=\Psi_0^{(+)}=\Psi_0
}
\]

とする。

initializer 自身は「乗法用」「加法用」という二種類の state を作ってはならない。

---

# 2. 上位 interaction v3.0 との state contract

初期化成功時の state ordering は上位仕様と一致させる。

\[
\boxed{
\Psi=
(\Psi_{unit},\Psi_{time},\Psi_{phys0},\Psi_{numerics},\Psi_{grid},\Psi_{field})^T
}
\]

## 2.1 unit

\[
\boxed{ONE=1}
\]

## 2.2 time state

必須成分：

\[
\boxed{T_\tau,\;\delta\tau,\;D_\tau,\;\kappa_\tau}
\]

定義：

\[
\boxed{T_\tau=e^{\kappa_\tau\tau}}
\]

\[
\boxed{D_\tau=e^{\kappa_\tau\delta\tau}}
\]

初期時刻 `tau0` から

\[
\boxed{T_{\tau,0}=e^{\kappa_\tau\tau_0}}
\]

とする。

`delta_tau,D_tau,kappa_tau` は初期化後の恒等 state であり、initializer は

\[
\boxed{
\delta\tau-\frac{\log D_\tau}{\kappa_\tau}=0
}
\]

を `epsilon_time_state_init` 内で監査する。

外部 `dt` は state へコピーせずに dynamics へ渡してはならない。

## 2.3 field state

GH–Maxwell 表現を採る run では

- `g_ab`
- `Pi_ab`
- `Phi_iab`
- `H_a`
- `Theta_a`
- `E^i`
- `B^i`

を一つの `u_0` として格納する。

---

# 3. 初期化と数値発展 solver の責任分離

interaction v3.0 は、採用した数値 Einstein–Maxwell solver の空間離散化・境界処理を `F_h` の一部としてそのまま使用する。

したがって initializer v3.0 は、interaction v2.0 のように独自の `D_i,V_*,P_*` を上位物理原則として要求しない。

代わりに次を要求する。

1. 初期場は、正式発展 solver が受理する **同一の field representation と grid representation** へ変換して格納する。
2. interpolation / projection / prolongation が必要なら、採用 solver の公式変換 routine を使う。
3. initializer が独自に作った hidden FD stencil / hidden RBF stencil を正式発展 field の定義へ混入しない。
4. evolution grid と initialization elliptic grid が異なる場合、その変換 algorithm・tolerance・version/hash を `Psi_initcfg` に恒等記録する。
5. 変換後に physical constraints を evolution representation 上で再評価する。

---

# 4. 全ケース共通の物理初期化法

binary black-hole experiment の物理初期化法は

\[
\boxed{
\text{superposed boosted Kerr–Newman free data}
\rightarrow
\text{coupled XCTS + electromagnetic constraints}
}
\]

とする。

重要：これは `Q=0`, `J=0`, 同符号電荷, 逆符号電荷, 非対称電荷, near-extremal の全 binary case で同じである。

禁止：

- MP exact solution への case switch
- double-RN exact solution への case switch
- neutral case だけ別 puncture initializer へ切替
- opposite-charge case だけ別式へ切替
- solver failure を別 initializer で救済

Kerr–Newman regularity gate：

\[
\boxed{
\left(\frac{Q_A^{(0)}}{M_A^{(0)}}\right)^2
+
\left(\frac{|J_A^{(0)}|}{M_A^{(0)2}}\right)^2
\le 1
}
\]

違反は `INITIALIZATION_REJECTED_SUPEREXTREMAL`。

---

# 5. N_BH=0,1,2 の扱いを分離して矛盾を除く

`N_BH` は topology input であり、case-specific physics switch ではない。

## 5.1 N_BH=0

Minkowski validation 専用。

- black-hole sum は空和。
- `d0,C0,Omega0,Omega_r,A` は **未使用**。
- horizon inner boundary は存在しない。
- free data は flat / zero EM。
- `u_0` は flat-space fixed point を作る。

## 5.2 N_BH=1

single Kerr–Newman validation。

- `d0,C0,Omega_guess` は **未使用**。
- center は `X_1^(0)`。
- translational boost `v_1` と spin `J_1` は明示 input からのみ生成。
- 1 個の horizon/excision boundary を持つ。

## 5.3 N_BH=2

正式 binary experiment。

この場合のみ

\[
M_{tot}^{(0)}=M_1^{(0)}+M_2^{(0)}
\]

\[
\boxed{d_0=P_0M_{tot}^{(0)}}
\]

\[
C_0=1-\lambda_1^{(0)}\lambda_2^{(0)}
\]

\[
\boxed{
\Omega_{guess}^2=\frac{M_{tot}^{(0)}C_0}{d_0^3}
}
\]

を定義する。

`C0<0` で `Omega_guess` が実数にならない場合は、quasi-circular binary baseline の domain 外として明示 FAIL し、絶対値化等をしない。

---

# 6. binary baseline coordinate state

等質量 baseline：

\[
\mathbf X_1^{(0)}=(+d_0/2,0,0),
\qquad
\mathbf X_2^{(0)}=(-d_0/2,0,0)
\]

center of mass：

\[
\mathbf r_{CM}
=\frac{\sum_A M_A^{(0)}\mathbf X_A^{(0)}}{M_{tot}^{(0)}}
\]

baseline：

\[
\boxed{\Omega_0=\Omega_{guess}},\qquad
\boxed{\dot a_0=0},\qquad
\boxed{\mathbf v_0=0}
\]

seed boost：

\[
\boxed{
\mathbf v_A
=\mathbf\Omega_0\times(\mathbf X_A^{(0)}-\mathbf r_{CM})
+\dot a_0(\mathbf X_A^{(0)}-\mathbf r_{CM})
+\mathbf v_0
}
\]

`|v_A|>=1` は FAIL。

---

# 7. Kerr–Newman seed

各 hole A の seed parameters：

\[
M_A^{seed},\quad Q_A^{seed}=Q_A^{(0)},\quad
\mathbf J_A^{seed}
\]

初期 outer-control guess：

\[
\boxed{M_A^{seed}=M_A^{(0)}}
\]

\[
\boxed{
\mathbf J_A^{seed}
=\boldsymbol\chi_A^{(0)}(M_A^{seed})^2
}
\]

個々の seed は Kerr–Newman in Kerr–Schild form の同一公式から生成し、`Q=0` と `J=0` はその連続極限とする。

seed の Lorentz boost は `v_A` により同じ tensor transformation で行う。

seed derivative の生成法は `seed_derivative_method_id` として恒等記録し、同一 run 中に変更しない。

---

# 8. attenuation

superposition attenuation は全 tensor/EM seed sectors で同一 rule を使う。

\[
W_A(x)=\exp[-r_A^2/w_A^2]
\]

\[
\boxed{
w_A=\zeta_{att}\frac{M_A^{seed}}{\sum_BM_B^{seed}}d_0
}
\]

`zeta_att` は `Psi_initcfg` の恒等 state。

binary 以外では `d0` ruleを使わない。single-hole validation では attenuation 自体を不要とし、`W_1=1` とする。Minkowski では空和。

---

# 9. XCTS + EM inner solve の unknowns を固定

**inner nonlinear solve 中に変更してよい未知量は場だけ**とする。

\[
\boxed{
U=(\psi,X,\beta^1,\beta^2,\beta^3,\varphi)
}
\]

ここで

\[
X=\alpha\psi
\]

である。

inner solve 中は以下を固定する。

- `M_A^seed`
- `Q_A^seed`
- `J_A^seed`
- `v0`
- `Omega0`
- `dot_a0`
- `Omega_r,A`
- excision surface geometry
- grid/map
- all solver tolerances

**Newton iteration 内でこれら control parameter を更新してはならない。**

---

# 10. XCTS free data

固定 control vector に対して free dataを一度構成する。

\[
\boxed{
\tilde\gamma_{ij}
=\delta_{ij}+\sum_AW_A(\gamma^A_{ij}-\delta_{ij})
}
\]

\[
\boxed{K=\sum_AW_AK_A}
\]

\[
\boxed{\tilde u_{ij}=0},\qquad
\boxed{\partial_tK=0}
\]

physical fields：

\[
\gamma_{ij}=\psi^4\tilde\gamma_{ij}
\]

\[
\tilde A^{ij}
=\frac{\psi^7}{2X}(\tilde L\beta)^{ij}
\]

\[
A^{ij}=\psi^{-10}\tilde A^{ij}
\]

\[
K_{ij}=A_{ij}+\frac13\gamma_{ij}K
\]

---

# 11. electromagnetic constraint unknown

EM elliptic unknown は一つの scalar correction

\[
\boxed{\varphi}
\]

とする。

\[
(E_{sp})^i=\sum_AW_AE_A^i
\]

\[
(A_{sp})_i=\sum_AW_A(A_A)_i
\]

\[
\boxed{E^i=(E_{sp})^i+D^i\varphi}
\]

\[
\boxed{
B^i=\frac1{\sqrt\gamma}[ikl]\partial_k(A_{sp})_l
}
\]

Gauss correction：

\[
\boxed{D_iD^i\varphi=-D_i(E_{sp})^i}
\]

---

# 12. coupled inner equations

inner solve は Hamiltonian, Momentum, lapse XCTS, electric Gauss correction を同時に解く。

Hamiltonian：

\[
\tilde D^2\psi
-\frac18\tilde R\psi
-\frac1{12}K^2\psi^5
+\frac18\psi^{-7}\tilde A_{ij}\tilde A^{ij}
=-2\pi\psi^5\rho
\]

Momentum：

\[
\tilde D_j\left[\frac{\psi^7}{2X}(\tilde L\beta)^{ij}\right]
-\frac23\psi^6\tilde D^iK
=8\pi\psi^{10}J^i
\]

Lapse：

\[
\tilde D^2X
-X\left[
\frac18\tilde R
+\frac5{12}K^2\psi^4
+\frac78\psi^{-8}\tilde A_{ij}\tilde A^{ij}
\right]
-\psi^5\beta^k\partial_kK
=2\pi X\psi^4(\rho+2\sigma)
\]

with `partial_t K=0`。

EM source：

\[
\rho=\frac{E^2+B^2}{8\pi},\qquad
J^i=\frac1{4\pi}(E\times B)^i,\qquad
\sigma=\rho
\]

を current inner iterate から毎回再評価する。

---

# 13. boundary conditions

outer infinity：

\[
\boxed{\psi\to1,\quad X\to1,\quad\beta^i\to v_0^i,\quad\varphi\to0}
\]

各 black-hole inner boundaryでは quasi-equilibrium horizon boundary condition を用いる。

shift：

\[
\boxed{
\beta^i
=\alpha s^i
-\Omega_{r,A}^k\xi^i_{(k),A}
-[\mathbf\Omega_0\times(\mathbf r-\mathbf r_{CM})]^i
-\dot a_0(\mathbf r-\mathbf r_{CM})^i
}
\]

conformal factor は apparent-horizon condition、lapseは指定した inner lapse gauge、electric correction は target horizon chargeを満たす Neumann conditionを用いる。

neutral targetでは `Q_target/Q_sp` の 0/0 を使わず、final normal electric flux zeroを直接課す。

---

# 14. excision surface と seed mass control の矛盾を解消

v2.0 の「inner boundary = 現在 seed Kerr–Newman horizon surface」は撤回する。

v3.0 では各 excision surface geometry

\[
\boxed{\mathcal S_A^{exc}}
\]

を `Psi_initcfg/Psi_grid` の**固定初期化設定**として与える。

`S_A^exc` は target physical parameters から事前構成された候補 surface であり、outer-control loop中に seed mass を変更しても移動させない。

要求：inner solve後の apparent horizon `H_A` が存在し、

\[
\boxed{\mathcal S_A^{exc}\subset\mathrm{interior}(H_A)}
\]

であることを `epsilon_excision_margin` 付きで監査する。

満たさない場合は `EXCISION_NOT_INSIDE_HORIZON` として FAIL。

**自動で excision surface を動かして救済しない。**

---

# 15. inner nonlinear solver

inner nonlinear solver ID、linear solver ID、preconditioner、restart length、norm、全 tolerance は `Psi_initcfg` に恒等記録する。

baseline algorithm：

```text
nonlinear_solver_id = NEWTON_RAPHSON
jacobian_method_id = FORWARD_AD
linear_solver_id = GMRES
preconditioner_id = BLOCK_JACOBI
line_search_id = ARMIJO_BACKTRACKING
```

GMRESについて最低限、

- `gmres_restart`
- `I_GMRES_max`
- `epsilon_linear_abs`
- `epsilon_linear_rel`
- `linear_norm_id`

をstate化する。

block-Jacobi block partitionも `preconditioner_block_schema_id` として固定する。

inner Newton は固定 control vectorのもとで `R(U)=0` を解き切る。

---

# 16. outer control solve を inner Newton から完全分離

binary の target matching/control vectorを

\[
\boxed{
C=
(M_1^{seed},M_2^{seed},\mathbf v_0,
\mathbf\Omega_{r,1},\mathbf\Omega_{r,2})
}
\]

とする。

spinless runでも `Omega_r,A=0` が初期値であり、同じ control schemaを使う。

inner solveを完全収束させた後だけ、outer residualを評価する。

\[
\boxed{
R_C=
(
M_{H,1}-M_1^{(0)},
M_{H,2}-M_2^{(0)},
\mathbf P_{ADM},
\mathbf J_{H,1}-\mathbf J_1^{(0)},
\mathbf J_{H,2}-\mathbf J_2^{(0)}
)
}
\]

2-hole caseでは未知 control数と residual数はともに11。

outer control更新は **inner solveの外側**で行う。

各 trial control vectorごとに、free dataを再構成し、inner XCTS+EM solveを最初から行う。

これにより v2.0 の「Newton iteration 中に boundary parameter を変える」矛盾を除去する。

outer Jacobian生成法・perturbation幅・line-search法も `Psi_initcfg` へ恒等記録する。

---

# 17. charge targetは outer control で変更しない

\[
\boxed{Q_A^{seed}=Q_A^{(0)}}
\]

を全 outer-control iteration で固定する。

charge matchingは EM inner boundary condition と constraint solveの結果として監査する。

\[
|Q_{H,A}-Q_A^{(0)}|<\epsilon_{charge_match}
\]

を満たさない場合に target chargeを変更してはならない。

---

# 18. convergence gates

inner residual：

\[
\|R_{XCTS+EM}\|_\infty<\epsilon_{init}
\]

physical constraintsを独立再評価：

\[
\|\mathcal H\|_\infty<\epsilon_H
\]

\[
\|\mathcal M_i\|_\infty<\epsilon_M
\]

\[
\|D_iE^i\|_\infty<\epsilon_E
\]

\[
\|D_iB^i\|_\infty<\epsilon_B
\]

horizon targets：

\[
|M_{H,A}-M_A^{(0)}|<\epsilon_{mass_match}
\]

\[
|Q_{H,A}-Q_A^{(0)}|<\epsilon_{charge_match}
\]

\[
\|\mathbf J_{H,A}-\mathbf J_A^{(0)}\|<\epsilon_{spin_match}
\]

ADM momentum：

\[
\|\mathbf P_{ADM}\|<\epsilon_{ADM_P}
\]

excision inclusionもPASSするまで `Psi_0` を返さない。

---

# 19. ADM initial data -> GH field state

inner/outer solve成功後の正本は

\[
\gamma_{ij},K_{ij},\alpha,\beta^i,E^i,B^i
\]

である。

spacetime metric：

\[
\boxed{g_{ij}=\gamma_{ij}}
\]

\[
\boxed{g_{0i}=\gamma_{ij}\beta^j}
\]

\[
\boxed{g_{00}=-\alpha^2+\gamma_{ij}\beta^i\beta^j}
\]

`Phi_iab` は evolution representationへ変換後、その solver が採用する spatial derivative operator で

\[
\boxed{\Phi_{iab}=\partial_i g_{ab}}
\]

を構成する。

`Pi_ab` は

\[
\Pi_{ab}=-n^c\partial_cg_{ab}
\]

と GH gauge constraint を同時に満たすよう構成し、`K_ij`再構成値との一致を監査する。

`Pi_ij` は

\[
\boxed{
\Pi_{ij}=2K_{ij}-n^a(\Phi_{ija}+\Phi_{jia})
}
\]

から固定する。

残る `Pi_00,Pi_0i` は initial gauge choice を与えた後に `C_a=0` の4方程式で決定する。

---

# 20. gauge state

interaction v3.0 の common RHS evaluator が採用する gauge prescriptionと**同一**の target `F_a` を使う。

初期値：

\[
\boxed{H_a^{(0)}=F_a^{(0)}}
\]

さらに初期 `partial_t H_a=0` を project choice とする場合、同じ gauge-driver equation から `Theta_a^(0)` を逆算する。

この gauge choiceは `gauge_initialization_id` として恒等記録し、加法/乗法両経路で共通とする。

---

# 21. evolution representation への変換後監査

初期化 elliptic representationから evolution solver representationへ変換した後、**その最終 representation 上で**次を再評価する。

1. `ONE=1`
2. `T_tau,delta_tau,D_tau,kappa_tau` consistency
3. metric signature
4. lapse positivity
5. GH constraints
6. Hamiltonian constraint
7. momentum constraint
8. Maxwell Gauss constraints
9. horizon M/Q/J targets
10. ADM momentum
11. excision inclusion
12. all identity states finite

変換前にPASSしていても、変換後にFAILなら `Psi_0` を返さない。

---

# 22. 乗法/加法二経路への接続監査

initializer outputは一つだけ。

run開始直前に

\[
\Psi_0^{(\times)}=\mathrm{copy}(\Psi_0)
\]

\[
\Psi_0^{(+)}=\mathrm{copy}(\Psi_0)
\]

として、

\[
\boxed{
\|\Psi_0^{(\times)}-\Psi_0^{(+)}\|_\infty=0
}
\]

を bitwise または exact-copy level で要求する。

初期化器が二経路で再実行されることを禁止する。乱数seed、solver roundoff、iteration orderの差を初期状態差へ持ち込まないためである。

---

# 23. 必須 identity/config state

最低限、以下を `Psi_0` に恒等保持する。

## physical labels

- `N_BH`
- initial masses, charges, spins
- initial positions
- `P0,E0,phi0,chi0,tau0`
- charge ratios and audit labels

## time

- `delta_tau`
- `D_tau`
- `kappa_tau`

## evolution numerics

- `time_integrator_id`
- RK tableau
- spatial solver / formulation ID
- boundary treatment IDs
- gauge IDs
- dissipation IDs and coefficients
- refinement / grid configuration identifiers

## initializer numerics

- `initializer_method_id = SKN_XCTS_EM_V3`
- attenuation settings
- inner nonlinear/linear solver IDs
- all absolute/relative tolerances
- iteration maxima
- line-search parameters
- outer-control solver parameters
- excision surface specification
- interpolation/projection method ID

隠れ default を正式runで使用しない。

---

# 24. FAIL codes

最低限：

```text
INITIALIZATION_REJECTED_SUPEREXTREMAL
INITIALIZATION_INVALID_BINARY_OMEGA
INITIALIZATION_SEED_SUPERLUMINAL
XCTS_NEWTON_DIVERGED
XCTS_LINEAR_SOLVE_FAILED
XCTS_LINE_SEARCH_FAILED
XCTS_CONSTRAINT_FAIL
OUTER_CONTROL_DIVERGED
HORIZON_MASS_MATCH_FAIL
HORIZON_CHARGE_MATCH_FAIL
HORIZON_SPIN_MATCH_FAIL
ADM_MOMENTUM_MATCH_FAIL
EXCISION_NOT_INSIDE_HORIZON
GH_CONVERSION_FAIL
GH_CONSTRAINT_INIT_FAIL
MAXWELL_CONSTRAINT_INIT_FAIL
EVOLUTION_REPRESENTATION_CONVERSION_FAIL
TIME_STATE_INIT_FAIL
INITIAL_STATE_COPY_MISMATCH
```

FAIL時に target physical values、time state、solver physicsを自動変更して続行してはならない。

---

# 25. v2.0から撤回する項目

v3.0では次を撤回する。

1. initializerが interaction v2.0 の `D_i,V_*,P_*` state contractへ直接従うこと。
2. `N_BH=0,1,2` に binary `d0,C0,Omega_guess` を共通適用すること。
3. Newton iteration内部で `v0` を更新すること。
4. Newton iteration内部で `Omega_r,A` を更新すること。
5. seed mass updateごとに excision surfaceを seed horizonへ一致させること。
6. interaction v2.0 / CF4 / K-blockとのcompatibility auditを初期化完成条件にすること。
7. 乗法系用・加法系用に initializerを別実行する余地。

---

# 26. v3.0の完成条件

初期化仕様が完成しているとは、承認済み physical inputs と numerical configuration に対して、追加判断なしに

\[
\boxed{\Psi_0\quad\text{または明示FAIL}}
\]

を一意に返し、その `Psi_0` が interaction v3.0 の

\[
\Psi_0\to F_h(\Psi)\to
\begin{cases}
\text{multiplicative formal route}\\
\text{additive validation route}
\end{cases}
\]

の**共通初期状態**になることをいう。

初期化は万能相互作用の乗法性を作らない。

初期化の責任は、物理拘束を満たし、時間状態・恒等状態を含んだ一つの完全状態を作るところまでである。

---

# 27. 実装開始前監査

正式コード実装前に以下を全PASSする。

- [ ] interaction v3.0 state schemaと一致
- [ ] ONEあり
- [ ] T_tau, delta_tau, D_tau, kappa_tauあり
- [ ] time consistency式あり
- [ ] binaryだけにd0/Omega ruleを適用
- [ ] N_BH=0/1でbinary式を参照しない
- [ ] inner Newton中にcontrol parameterを更新しない
- [ ] outer control loopがinner solveの外側
- [ ] seed massとexcision geometryが矛盾しない
- [ ] one initializer output only
- [ ] multiplication/addition branches are exact copies of same Psi0
- [ ] evolution representation上でconstraintsを再監査
- [ ] hidden fallbackなし
- [ ] hidden dtなし

一つでも未定義なら正式実験コードを書いてはならない。

---

# 28. 最終定義

v3.0初期化は

\[
\boxed{
\text{physical targets + init config}
\rightarrow
\text{common constrained Einstein–Maxwell initial data}
\rightarrow
\text{evolution field representation}
\rightarrow
\Psi_0
}
\]

である。

そして

\[
\boxed{
\Psi_0^{(\times)}=\Psi_0^{(+)}=\Psi_0
}
\]

を唯一の実験開始条件とする。

**この初期化仕様以外の経路で作った state、または二経路で別々に初期化した state は正式実験に使用しない。**
