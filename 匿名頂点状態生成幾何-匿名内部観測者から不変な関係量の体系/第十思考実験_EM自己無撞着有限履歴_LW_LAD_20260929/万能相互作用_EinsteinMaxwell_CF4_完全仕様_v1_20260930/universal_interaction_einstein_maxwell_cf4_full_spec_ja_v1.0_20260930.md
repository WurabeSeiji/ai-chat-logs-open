# 万能相互作用 Einstein–Maxwell–GH–CF4 完全仕様
## ― 初期化・全状態ベクトル ψ・状態関数 S(ψ)・名前付き配列添字・読み出し・監査条件 ―

**著者:** 木原範昭 (Noriaki Kihara)  
**仕様版:** v1.0  
**日付:** 2026-09-30  
**目的:** 過去の第六・第七思考実験の初期条件を出発点として、ケース別の相互作用則を一切使わず、full nonlinear Einstein–Maxwell を唯一の物理方程式系とし、現在状態だけから毎 step 一つの万能相互作用配列 `S(ψ)` を生成して `ψ' = S(ψ) ψ` と更新する数値実験の完全仕様を定める。

---

# 0. 最上位原則

本仕様では、以下を絶対条件とする。

1. **唯一の永続情報は状態ベクトル `ψ` に明示された成分だけである。**
2. 初期化時に与える値は、原則すべて `ψ` の成分として保持する。
3. 初期化後に値が変わらないものも、外部定数として隠さず、`S(ψ)` の恒等行として保持する。
4. Einstein–Maxwell の現在場から導出できる量は、別の永続レジスタに複製せず、その都度 `ψ` から導出する。
5. `a,b` の現在質量・電荷・スピン・位置・軌道は **readout** とし、動力学へ戻さない。
6. ただし初期化時に指定した `M_a^(0),Q_a^(0),J_a^(0),...` は初期化記録として `ψ` に恒等保持する。
7. 相互作用則はケースごとに変えない。全ケースで同じ
   
   `K(ψ) -> CF4 -> S(ψ)`
   
   を使用する。
8. `S` は固定行列ではない。現在状態から生成される状態関数である。

\[
\boxed{S_n=S(\psi_n)}
\]

\[
\boxed{\psi_{n+1}=S(\psi_n)\,\psi_n}
\]

9. `S(ψ)` を生成する途中の `K^[1]...K^[4]`、CF4 中間状態、逆行列、Christoffel 記号、指数作用の Krylov 基底等は **一時量**であり、step 境界を越えて保持しない。
10. 数値解像度・基底次数・座標写像・ゲージ係数・拘束減衰係数・時間刻み・solver tolerance など、次状態生成または再現性に影響する初期設定はすべて `ψ` に保持する。
11. 連続場に無限次数が必要であっても、実験では有限 cutoff を初期化する。cutoff 自身を状態として保持し、cutoff を変えた独立走行の収束性を実験で判定する。
12. 読み出しは一方向である。

\[
\boxed{\psi\longrightarrow \mathcal R(\psi)}
\]

readout を `S` 生成へ帰還させない。

---

# 1. 採用する唯一の物理系

採用する物理系は **Einstein–Maxwell electrovac** である。

計量符号は

\[
(-,+,+,+)
\]

とする。

Einstein 方程式：

\[
G_{ab}=8\pi G\,T^{\rm EM}_{ab}.
\]

Maxwell 方程式：

\[
\nabla_a F^{ab}=0,
\qquad
\nabla_{[a}F_{bc]}=0.
\]

電磁 stress-energy：

\[
T^{\rm EM}_{ab}
=\frac{1}{4\pi}
\left(
F_{ac}F_b{}^c
-\frac14 g_{ab}F_{cd}F^{cd}
\right),
\]

かつ

\[
T^{\rm EM\,a}{}_a=0.
\]

重力発展は first-order generalized harmonic (GH) 系を使う。

\[
\Pi_{ab}=-n^c\partial_c g_{ab},
\qquad
\Phi_{iab}=\partial_i g_{ab}.
\]

ゲージ source `H_a` は外部関数にせず、first-order gauge driver の動的状態とする。

---

# 2. 添字規約

- spacetime index: `a,b,c,d,e,f = 0,1,2,3`
- spatial index: `i,j,k,l = 1,2,3`
- spectral domain index: `d = 1,...,N_domain`
- modal multi-index: `λ,μ ∈ Λ_d(N_cut)`
- 状態変数名をそのまま `S` と `K` の行・列添字に用いる。

例えば、

\[
S_{g_{01}[d,\lambda],\Pi_{23}[e,\mu]}(\psi)
\]

\[
K_{E^2[d,\lambda],B^3[e,\mu]}(\psi)
\]

のように書く。

数値実装時に整数 index へ変換してよいが、仕様・監査・ログでは必ず名前付き添字を保持する。

---

# 3. 完全状態ベクトル ψ

状態ベクトルを

\[
\boxed{
\psi=
\left(
\psi_{\rm unit},
\psi_{\rm time},
\psi_{\rm phys0},
\psi_{\rm numerics},
\psi_{\rm grid},
\psi_{\rm field}
\right)^T
}
\]

と定義する。

## 3.1 unit state

\[
\boxed{\mathrm{ONE}=1}
\]

を一成分として保持する。

非線形 local source を canonical affine lift で一意に行列化するために使用する。

---

## 3.2 時間状態

加法時刻を直接保持する代わりに、

\[
T_\tau=e^{\kappa_\tau\tau},
\qquad
D_\tau=e^{\kappa_\tau\Delta\tau}
\]

を状態に持つ。

`κ_τ` も恒等状態である。

raw 時刻と raw step は readout：

\[
\tau=\frac{\log T_\tau}{\kappa_\tau},
\qquad
\Delta\tau=\frac{\log D_\tau}{\kappa_\tau}.
\]

`D_τ` は v1.0 では恒等状態とする。

---

## 3.3 初期物理条件の恒等状態 ψ_phys0

初期化に使用した値を全て保存する。

\[
\psi_{\rm phys0}=\{
q_0,
G,c,
M_{\rm tot}^{(0)},
M_a^{(0)},M_b^{(0)},
Q_a^{(0)},Q_b^{(0)},
\mathbf J_a^{(0)},\mathbf J_b^{(0)},
\mathbf X_a^{(0)},\mathbf X_b^{(0)},
\mathbf P_a^{(0)},\mathbf P_b^{(0)},
P_0,E_0,\phi_0,\chi_0,\tau_0,
\nu_0,n_q,s_a,s_b,
\lambda_a^{(0)},\lambda_b^{(0)},C_0,D_0
\}.
\]

これらは初期化後、原則として恒等更新する。

重要：

`M_a^(0),Q_a^(0),J_a^(0)` 等は **現在値ではない**。

現在値は場から

\[
M_a^{\rm read}(\psi),
Q_a^{\rm read}(\psi),
\mathbf J_a^{\rm read}(\psi)
\]

として読む。

---

## 3.4 数値・方程式設定の恒等状態 ψ_numerics

次状態生成または再現性に影響する設定値を全て保持する。

### GH constraint damping

\[
\gamma_0,\gamma_1,\gamma_2,
\gamma_3,\gamma_4,\gamma_5,\epsilon_5.
\]

baseline は

\[
\gamma_3=\gamma_4=\gamma_5=0
\]

として Lindblom–Scheel–Kidder–Owen–Rinne 系を回復する。

### gauge driver

\[
\mu_H,\eta_H,\mu_L,\mu_S,p_H.
\]

必要なら gauge blend 用

\[
T_{\rm gauge-blend}
\]

も状態化する。

### finite expansion

\[
N_{\rm domain},
\{N_{r,d},L_d,M_d\}_{d=1}^{N_{\rm domain}},
\mathrm{basis\_id}_d,
\mathrm{quadrature\_id}_d.
\]

`N_cut` を単一値で用いる実装では

\[
N_{\rm cut}
\]

を恒等状態にする。

### 初期 constraint solver

\[
\epsilon_{\rm init},
I_{\rm init,max},
\epsilon_{\rm horizon-match},
I_{\rm horizon-match,max}.
\]

### matrix exponential action

\[
\epsilon_{\exp},
 m_{\exp,\max},
\mathrm{exp\_method\_id}.
\]

### run control（状態列の生成自体には帰還しない）

\[
N_{\rm step,max},
N_{\rm readout,stride},
r_{\rm readout,stop}.
\]

これらは監査記録として `ψ` に保持するが、`S` の物理 block に帰還しない。

---

## 3.5 グリッド・domain・座標写像状態 ψ_grid

固定グリッドを外部に置かない。

各 spectral domain の参照座標 `ξ^A` から inertial coordinate `x^i` への写像を

\[
x^i=\mathcal M^i_d(\xi)
\]

とし、その係数を全て状態にする。

\[
\mathcal M^i_d(\xi)
=\sum_{\lambda\in\Lambda_d}
\mathcal M^i[d,\lambda]\,\Phi_{d,\lambda}(\xi).
\]

`v1.0` では map 係数は **恒等状態**として保持する。

したがって grid は固定された外部物ではなく、`ψ` に明示された不変状態である。

外側 domain は rational/Chebyshev compactification で空間無限遠を表す。

内側は各 black hole の excision surface を domain map の状態として保持する。

inner boundary は horizon 内に置き、全 characteristic が domain 外へ流出することを validity gate とする。incoming characteristic が出現した場合は、その走行を `INVALID_INNER_BOUNDARY` として停止し、外部境界値を注入しない。

---

## 3.6 動的場状態 ψ_field

各 `d,λ` について、

### GH gravity

\[
g_{ab}[d,\lambda]\quad(10)
\]

\[
\Pi_{ab}[d,\lambda]\quad(10)
\]

\[
\Phi_{iab}[d,\lambda]\quad(30)
\]

### gauge driver

\[
H_a[d,\lambda]\quad(4)
\]

\[
\Theta_a[d,\lambda]\quad(4)
\]

### Maxwell

\[
E^i[d,\lambda]\quad(3)
\]

\[
B^i[d,\lambda]\quad(3)
\]

計 64 field components / mode とする。

map 3 components/mode を含めれば、状態として domain ごとに 67 components/mode を保持する。

---

# 4. 過去実験から継承する baseline 初期条件

第六思考実験補遺で検証済みの規格化を baseline にする。

\[
\boxed{G=q_0=c=1}
\]

\[
\boxed{M_{\rm tot}/q_0=20/3}
\]

等質量：

\[
M_a^{(0)}=M_b^{(0)}=10/3.
\]

\[
\nu_0=1/4.
\]

軌道条件：

\[
P_0=50,
\quad
E_0=0,
\quad
\phi_0=0,
\quad
\chi_0=0,
\quad
\tau_0=0.
\]

旧比較分解能：

\[
N_{\phi,0}=4000,
\qquad
h_{\phi,0}=2\pi/4000.
\]

この `N_phi0,h_phi0` も初期化記録 state として保持するが、full Einstein–Maxwell の時間発展 step とは同一視しない。

---

# 5. 過去5条件の電荷 baseline

\[
\lambda_A=Q_A/M_A,
\qquad
\lambda_B=Q_B/M_B.
\]

| case_id | n_q | s_a | s_b | λ_a | λ_b | C_0 | D_0 |
|---|---:|---:|---:|---:|---:|---:|---:|
| n1_attractive | 1 | +1 | -1 | +0.30 | -0.30 | 1.09 | 0.36 |
| n1_repulsive  | 1 | +1 | +1 | +0.30 | +0.30 | 0.91 | 0 |
| n3_attractive | 3 | +1 | -1 | +0.90 | -0.90 | 1.81 | 3.24 |
| n3_repulsive  | 3 | +1 | +1 | +0.90 | +0.90 | 0.19 | 0 |
| n4_attractive | 4 | +1 | -1 | +1.20 | -1.20 | 2.44 | 5.76 |

\[
C_0=1-\lambda_a\lambda_b,
\qquad
D_0=(\lambda_a-\lambda_b)^2.
\]

`C_0,D_0` は旧実験との audit 用恒等状態であり、full Einstein–Maxwell evolution の係数へ直接再注入しない。

---

# 6. Kerr–Newman 初期化可能性ゲート

full Einstein–Maxwell 初期化対象を regular Kerr–Newman black hole とする場合、各 hole について

\[
M_A^2\ge Q_A^2+a_A^2,
\qquad
a_A=J_A/M_A
\]

すなわち

\[
\boxed{
\left(\frac{Q_A}{M_A}\right)^2
+
\left(\frac{J_A}{M_A^2}\right)^2
\le1
}
\]

を満たさなければならない。

よって baseline spin=0 では、

- n1: `|Q|/M=0.30` → PASS
- n3: `|Q|/M=0.90` → PASS
- n4: `|Q|/M=1.20` → FAIL

である。

n4 条件は旧 point-particle 実験条件として `ψ_phys0` に保存可能だが、regular Kerr–Newman 初期化には使用しない。

勝手な `|Q|` 化、質量変更、複素化、extremality 超過の救済はしない。

spin extension では、

\[
\chi_A^{(0)}=J_A^{(0)}/M_A^{(0)2}
\]

として

\[
|\chi_A^{(0)}|
\le\sqrt{1-(Q_A^{(0)}/M_A^{(0)})^2}
\]

を必須 gate とする。

---

# 7. 旧 P0 から binary coordinate data への lift

旧 `P_0=50` は reduced one-body 軌道の dimensionless separation scale と解釈し、

\[
R_0=P_0 M_{\rm tot}^{(0)}
\]

を初期 coordinate separation とする。

等質量 center-of-mass frame：

\[
\mathbf X_a^{(0)}=(+R_0/2,0,0),
\qquad
\mathbf X_b^{(0)}=(-R_0/2,0,0).
\]

baseline は準円初期 guess として

\[
\Omega_{\rm guess}^2
=
\frac{G M_{\rm tot}^{(0)} C_0}{R_0^3}
\]

を用いる。

\[
\mu_0=\nu_0 M_{\rm tot}^{(0)}
\]

\[
p_{t,\rm guess}=\mu_0 R_0\Omega_{\rm guess}.
\]

\[
\mathbf P_a^{(0)}=(0,+p_{t,\rm guess},0),
\qquad
\mathbf P_b^{(0)}=(0,-p_{t,\rm guess},0).
\]

これは **XCTS solver の initial guess** であり、最終 Einstein–Maxwell 場を old LO 方程式で拘束しない。

---

# 8. Einstein–Maxwell 初期場生成

初期場は superposed Kerr–Newman free data + extended conformal thin sandwich (XCTS) + electromagnetic constraints を用いて生成する。

入力は `ψ_phys0, ψ_numerics, ψ_grid` の現在値だけである。

初期化 solver は次の拘束を同時に満たす。

Hamiltonian constraint：

\[
\boxed{
{}^{(3)}R+K^2-K_{ij}K^{ij}
=16\pi G\,\rho_{\rm EM}
}
\]

Momentum constraint：

\[
\boxed{
D_j(K^{ij}-\gamma^{ij}K)
=8\pi G\,S^i_{\rm EM}
}
\]

Maxwell Gauss constraints：

\[
\boxed{D_iE^i=0}
\]

\[
\boxed{D_iB^i=0}
\]

（electrovac computational domain 内）。

初期 individual horizon の target は

\[
(M_A^{\rm read},Q_A^{\rm read},\mathbf J_A^{\rm read})
\]

が

\[
(M_A^{(0)},Q_A^{(0)},\mathbf J_A^{(0)})
\]

へ tolerance 内で一致するよう initializer seed を反復補正する。

収束条件：

\[
\|\mathcal H\|_\infty<\epsilon_{\rm init},
\]

\[
\|\mathcal M_i\|_\infty<\epsilon_{\rm init},
\]

\[
\|D_iE^i\|_\infty<\epsilon_{\rm init},
\qquad
\|D_iB^i\|_\infty<\epsilon_{\rm init},
\]

かつ horizon parameter mismatch が `ε_horizon-match` 未満。

達しない場合は走行を開始しない。

---

# 9. ADM initial data から GH state への変換

XCTS が生成した

\[
\gamma_{ij},K_{ij},\alpha,\beta^i,E^i,B^i
\]

から spacetime metric `g_ab` を作る。

\[
g_{00}=-\alpha^2+\gamma_{ij}\beta^i\beta^j,
\]

\[
g_{0i}=\gamma_{ij}\beta^j,
\qquad
g_{ij}=\gamma_{ij}.
\]

\[
\Phi_{iab}=\partial_i g_{ab}.
\]

\[
\Pi_{ab}=-n^c\partial_cg_{ab}
\]

は ADM initial data と選択した初期 lapse/shift time derivative から構成する。

初期 gauge constraint を満たすよう

\[
\boxed{H_a^{(0)}=-\Gamma_a^{(0)}}
\]

とする。

first-order gauge driver の

\[
\partial_t\Theta_a+\eta_H\Theta_a
=-\eta_H\beta^k\partial_kH_a
\]

について初期 `∂t Θ=0` を採用する場合、

\[
\boxed{
\Theta_a^{(0)}=-\beta^k\partial_kH_a^{(0)}
}
\]

とする。

最後に全場を現在の有限 modal basis へ projection し、

\[
g_{ab}[d,\lambda],\Pi_{ab}[d,\lambda],\Phi_{iab}[d,\lambda],H_a[d,\lambda],\Theta_a[d,\lambda],E^i[d,\lambda],B^i[d,\lambda]
\]

を `ψ_0` に格納する。

---

# 10. 初期有限次数の定義

連続初期場 `u(x)` を

\[
u^{(N_0)}(x)
=\sum_{d,\lambda\in\Lambda_d(N_0)}
u[d,\lambda]\Phi_{d,\lambda}(x)
\]

で近似する。

初期 cutoff `N_0` は `ψ_numerics` の恒等状態である。

初期化成功条件には constraint residual だけでなく projection convergence を含める。

例えば独立走行

\[
N_0,2N_0,4N_0
\]

で horizon parameter と constraint norm が収束することを比較する。

ただし、必要次数が有限で十分か、時間発展中に高次が必要になるかは **事前仮定せず数値実験で判定する**。

---

# 11. `ψ` から導出する一時幾何量

以下は永続状態にしない。

\[
g^{ab},\gamma_{ij},\gamma^{ij},\alpha,\beta^i,n^a,n_a,
\sqrt\gamma,
{}^{(4)}\Gamma^a{}_{bc},
{}^{(3)}\Gamma^i{}_{jk},
K_{ij},K,
\Gamma_a,
\mathcal C_a,
\mathcal C_{iab}
\]

など。

GH gauge constraint：

\[
\boxed{\mathcal C_a=H_a+\Gamma_a}
\]

three-index constraint：

\[
\boxed{
\mathcal C_{iab}=D_ig_{ab}-\Phi_{iab}
}
\]

を毎 stage 現在状態から再計算する。

---

# 12. spectral differential operator の生成

参照 basis derivative matrix を

\[
D_A^{(d)}[N_d]
\]

とする。

現在の map state から Jacobian

\[
J^i{}_A=\partial x^i/\partial\xi^A
\]

と inverse Jacobian

\[
(J^{-1})^A{}_i
\]

を計算し、物理空間微分を

\[
\boxed{
\mathscr D_i(\psi)
=(J^{-1})^A{}_i\,D_A
}
\]

として生成する。

modal multiplication operator を

\[
\boxed{
\mathscr M[f(\psi)]_{[d,\lambda],[e,\mu]}
}
\]

と書く。

これは `f(x;ψ)u(x)` を modal projection した線形作用素である。

---

# 13. canonical affine generator K(ψ)

非線形方程式 RHS の行列化に恣意性を残さないため、規則を固定する。

各動的 field `u_A` について

\[
\partial_tu_A
=
\sum_B \mathcal L_A{}^B(\psi)u_B
+R_A(\psi)
\]

と整理する。

- differential operator が直接作用する項 → `L`
- current-state から局所評価できる非線形項 → `R`

そして

\[
\boxed{
K(\psi)=
\begin{pmatrix}
\mathcal L(\psi) & R(\psi)\\
0 & 0
\end{pmatrix}
}
\]

とする。

最後の列は `ONE` 列である。

したがって

\[
\boxed{K(\psi)\psi=F(\psi)}
\]

を stage ごとに機械監査できる。

---

# 14. 名前付き K 行列 ― 恒等状態

任意の恒等状態 `C` について continuous generator は

\[
\boxed{K_{C,B}=0\quad\forall B}
\]

とする。

したがって exponential により

\[
S_{C,C}=1,
\qquad
S_{C,B\ne C}=0
\]

となる。

対象：

- `ONE`
- `D_tau`
- `kappa_tau`
- 全 `ψ_phys0`
- 全数値設定
- 全 grid/map coefficients（v1.0）

である。

---

# 15. 時間状態の K 行

\[
\dot T_\tau=\kappa_\tau T_\tau.
\]

よって

\[
\boxed{
K_{T_\tau,T_\tau}(\psi)=\kappa_\tau
}
\]

のみ非零。

一 step `Δτ` では

\[
T_\tau'=e^{\kappa_\tau\Delta\tau}T_\tau=D_\tau T_\tau.
\]

---

# 16. GH evolution equations ― 完全採用形

mesh velocity は v1.0 では `0`（map state 恒等更新）とする。

## 16.1 g_ab

\[
\boxed{
\partial_tg_{ab}
=(1+\gamma_1)\beta^kD_kg_{ab}
-\alpha\Pi_{ab}
-\gamma_1\beta^i\Phi_{iab}
}
\]

よって modal named blocks：

\[
\boxed{
K_{g_{ab}[d,\lambda],g_{cd}[e,\mu]}
=
\delta_a^{(c}\delta_b^{d)}
\left[\mathscr M[(1+\gamma_1)\beta^k]\mathscr D_k\right]_{d\lambda,e\mu}
}
\]

local source：

\[
R_{g_{ab}}
=-\alpha\Pi_{ab}-\gamma_1\beta^i\Phi_{iab}
\]

\[
\boxed{
K_{g_{ab}[d,\lambda],\mathrm{ONE}}
=\mathcal P_{d\lambda}[R_{g_{ab}}]
}
\]

その他の `g` 行要素は 0。

---

## 16.2 Pi_ab

採用式：

\[
\begin{aligned}
\partial_t\Pi_{ab}
={}&\beta^kD_k\Pi_{ab}
-\alpha\gamma^{ki}D_k\Phi_{iab}
+\gamma_1\gamma_2\beta^kD_kg_{ab}\\
&+R^{\rm GH}_{\Pi,ab}
-16\pi G\alpha T^{\rm EM}_{ab}.
\end{aligned}
\]

Maxwell stress tensor は trace-free のため trace reverse は同一。

GH local source：

\[
\begin{aligned}
R^{\rm GH}_{\Pi,ab}
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
\frac{
\mathcal C_a\mathcal C_b-rac12g_{ab}\mathcal C_d\mathcal C^d
}{
\epsilon_5+2(n^d\mathcal C_d)^2+\mathcal C_d\mathcal C^d
}\\
&-\gamma_1\gamma_2\beta^i\Phi_{iab}.
\end{aligned}
\]

baseline は `γ3=γ4=γ5=0`。

named derivative blocks：

\[
\boxed{
K_{\Pi_{ab}[d,\lambda],\Pi_{cd}[e,\mu]}
=
\delta_a^{(c}\delta_b^{d)}
[\mathscr M[\beta^k]\mathscr D_k]_{d\lambda,e\mu}
}
\]

\[
\boxed{
K_{\Pi_{ab}[d,\lambda],\Phi_{icd}[e,\mu]}
=
-\delta_a^{(c}\delta_b^{d)}
[\mathscr M[\alpha\gamma^{ki}]\mathscr D_k]_{d\lambda,e\mu}
}
\]

\[
\boxed{
K_{\Pi_{ab}[d,\lambda],g_{cd}[e,\mu]}
=
\delta_a^{(c}\delta_b^{d)}
[\mathscr M[\gamma_1\gamma_2\beta^k]\mathscr D_k]_{d\lambda,e\mu}
}
\]

local column：

\[
\boxed{
K_{\Pi_{ab}[d,\lambda],\mathrm{ONE}}
=
\mathcal P_{d\lambda}
\left[
R^{\rm GH}_{\Pi,ab}-16\pi G\alpha T^{\rm EM}_{ab}
\right]
}
\]

`∇_(a H_b)` に必要な空間微分は `D_iH_a`、時間微分は gauge-driver RHS から現在 stage 内で一時的に評価する。

---

## 16.3 Phi_iab

\[
\boxed{
\begin{aligned}
\partial_t\Phi_{iab}
={}&\beta^kD_k\Phi_{iab}
-\alpha D_i\Pi_{ab}
+\alpha\gamma_2D_ig_{ab}\\
&+\frac12\alpha n^cn^d\Phi_{icd}\Pi_{ab}
+\alpha\gamma^{jk}n^c\Phi_{ijc}\Phi_{kab}
-\alpha\gamma_2\Phi_{iab}.
\end{aligned}
}
\]

named derivative blocks：

\[
\boxed{
K_{\Phi_{iab}[d,\lambda],\Phi_{jcd}[e,\mu]}
=
\delta_i^j\delta_a^{(c}\delta_b^{d)}
[\mathscr M[\beta^k]\mathscr D_k]_{d\lambda,e\mu}
}
\]

\[
\boxed{
K_{\Phi_{iab}[d,\lambda],\Pi_{cd}[e,\mu]}
=
-\delta_a^{(c}\delta_b^{d)}
[\mathscr M[\alpha]\mathscr D_i]_{d\lambda,e\mu}
}
\]

\[
\boxed{
K_{\Phi_{iab}[d,\lambda],g_{cd}[e,\mu]}
=
\delta_a^{(c}\delta_b^{d)}
[\mathscr M[\alpha\gamma_2]\mathscr D_i]_{d\lambda,e\mu}
}
\]

local：

\[
\boxed{
K_{\Phi_{iab}[d,\lambda],\mathrm{ONE}}
=
\mathcal P_{d\lambda}[R_{\Phi,iab}]
}
\]

with

\[
R_{\Phi,iab}
=
\frac12\alpha n^cn^d\Phi_{icd}\Pi_{ab}
+\alpha\gamma^{jk}n^c\Phi_{ijc}\Phi_{kab}
-\alpha\gamma_2\Phi_{iab}.
\]

---

# 17. gauge-driver equations

\[
\boxed{
\partial_tH_a-\beta^kD_kH_a
=-\mu_H(H_a-F_a)+\Theta_a
}
\]

\[
\boxed{
\partial_t\Theta_a+\eta_H\Theta_a
=-\eta_H\beta^kD_kH_a
}
\]

Damped-wave target：

\[
\boxed{
F_a
=\mu_L\log\left(\frac{\gamma^{p_H}}{\alpha}\right)n_a
-\mu_S\alpha^{-1}g_{ai}\beta^i
}
\]

where `γ=det γ_ij`。

named blocks：

\[
\boxed{
K_{H_a[d,\lambda],H_b[e,\mu]}
=\delta_a^b
[\mathscr M[\beta^k]\mathscr D_k]_{d\lambda,e\mu}
}
\]

\[
\boxed{
K_{H_a[d,\lambda],\mathrm{ONE}}
=
\mathcal P_{d\lambda}
[-\mu_H(H_a-F_a)+\Theta_a]
}
\]

\[
\boxed{
K_{\Theta_a[d,\lambda],H_b[e,\mu]}
=-\delta_a^b
[\mathscr M[\eta_H\beta^k]\mathscr D_k]_{d\lambda,e\mu}
}
\]

\[
\boxed{
K_{\Theta_a[d,\lambda],\mathrm{ONE}}
=\mathcal P_{d\lambda}[-\eta_H\Theta_a]
}
\]

---

# 18. Maxwell decomposition

normal observer に対して

\[
E_a=F_{ab}n^b,
\qquad
B_a={}^*F_{ab}n^b,
\qquad
E_an^a=B_an^a=0.
\]

spatial Levi-Civita tensor を `ε^{ijk}` とする。

Electrovac evolution：

\[
\boxed{
(\partial_t-\mathcal L_\beta)E^i
=\epsilon^{ijk}D_j(\alpha B_k)+\alpha K E^i
}
\]

\[
\boxed{
(\partial_t-\mathcal L_\beta)B^i
=-\epsilon^{ijk}D_j(\alpha E_k)+\alpha K B^i
}
\]

constraints：

\[
\boxed{D_iE^i=0}
\]

\[
\boxed{D_iB^i=0}
\]

を monitoring し、v1.0 では別の divergence-cleaning persistent state を追加しない。

simulation が有限 cutoff で constraints を保持できないことを示した場合、その結果自体を保存し、後続版で cleaning state 追加を検討する。

---

# 19. Maxwell named K blocks

Lie derivative：

\[
(\mathcal L_\beta E)^i
=\beta^jD_jE^i-E^jD_j\beta^i.
\]

したがって derivative-on-field block は

\[
\boxed{
K_{E^i[d,\lambda],E^j[e,\mu]}
=
\delta^i_j
[\mathscr M[\beta^k]\mathscr D_k]_{d\lambda,e\mu}
}
\]

curl B block：

\[
\boxed{
K_{E^i[d,\lambda],B^l[e,\mu]}
=
[\mathscr M[\alpha\epsilon^{ijk}\gamma_{kl}]\mathscr D_j]_{d\lambda,e\mu}
}
\]

local electric source：

\[
R_E^i
=-E^jD_j\beta^i
+\epsilon^{ijk}(D_j\alpha)B_k
+\alpha K E^i.
\]

\[
\boxed{
K_{E^i[d,\lambda],\mathrm{ONE}}
=\mathcal P_{d\lambda}[R_E^i]
}
\]

magnetic：

\[
\boxed{
K_{B^i[d,\lambda],B^j[e,\mu]}
=
\delta^i_j
[\mathscr M[\beta^k]\mathscr D_k]_{d\lambda,e\mu}
}
\]

\[
\boxed{
K_{B^i[d,\lambda],E^l[e,\mu]}
=
-[\mathscr M[\alpha\epsilon^{ijk}\gamma_{kl}]\mathscr D_j]_{d\lambda,e\mu}
}
\]

local magnetic source：

\[
R_B^i
=-B^jD_j\beta^i
-\epsilon^{ijk}(D_j\alpha)E_k
+\alpha K B^i.
\]

\[
\boxed{
K_{B^i[d,\lambda],\mathrm{ONE}}
=\mathcal P_{d\lambda}[R_B^i]
}
\]

---

# 20. Maxwell stress-energy の現在場からの構成

\[
\rho_{\rm EM}
=\frac{E_iE^i+B_iB^i}{8\pi}
\]

\[
S_i^{\rm EM}
=\frac{1}{4\pi}\epsilon_{ijk}E^jB^k
\]

\[
S_{ij}^{\rm EM}
=\frac{1}{4\pi}
\left[
-E_iE_j-B_iB_j
+\frac12\gamma_{ij}(E^2+B^2)
\right].
\]

これらから

\[
T^{\rm EM}_{ab}
=\rho n_an_b+2n_{(a}S_{b)}+S_{ab}
\]

を現在 stage の `g,E,B` だけから構成する。

初期化値 `Q_a^(0),Q_b^(0)` は使用しない。

---

# 21. 初期化専用恒等値から動力学への直接結合禁止

初期化後、以下を要求する。

任意の field row `u` と initial-parameter column `I0` について

\[
\boxed{K_{u,I0}=0}
\]

対象例：

\[
K_{g_{ab},M_a^{(0)}}=0,
\]

\[
K_{\Pi_{ab},J_{a,x}^{(0)}}=0,
\]

\[
K_{E^i,Q_a^{(0)}}=0,
\]

\[
K_{B^i,C_0}=K_{B^i,D_0}=0.
\]

物理情報はすべて初期場に写像済みであり、その後は場だけを進化させる。

---

# 22. K(ψ) 完全非零 block 一覧

field 部分で非零になり得る named block は次だけである。

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

恒等状態の generator row は全て zero。

この一覧外に非零 entry が出た場合、実装監査 FAIL とする。

---

# 23. 4次 commutator-free Lie/CF4 による S(ψ) 生成

低次暫定法は作らない。

最初から Celledoni–Marthinsen–Owren 4次 commutator-free scheme を唯一の time update とする。

step 幅：

\[
h=\Delta\tau=\log D_\tau/\kappa_\tau.
\]

stage 1：

\[
\psi^{[1]}=\psi_n,
\qquad
K^{[1]}=K(\psi^{[1]}).
\]

stage 2：

\[
\psi^{[2]}
=\exp\left(\frac h2K^{[1]}\right)\psi_n,
\]

\[
K^{[2]}=K(\psi^{[2]}).
\]

stage 3：

\[
\psi^{[3]}
=\exp\left(\frac h2K^{[2]}\right)\psi_n,
\]

\[
K^{[3]}=K(\psi^{[3]}).
\]

stage 4：

\[
\psi^{[4]}
=\exp\left[h\left(K^{[3]}-\frac12K^{[1]}\right)\right]\psi^{[2]},
\]

\[
K^{[4]}=K(\psi^{[4]}).
\]

first final exponential：

\[
A
=\frac h{12}
\left(
3K^{[1]}+2K^{[2]}+2K^{[3]}-K^{[4]}
\right).
\]

\[
\psi_{n+1/2}=e^A\psi_n.
\]

second final exponential：

\[
B
=\frac h{12}
\left(
-K^{[1]}+2K^{[2]}+2K^{[3]}+3K^{[4]}
\right).
\]

\[
\boxed{
\psi_{n+1}=e^B e^A\psi_n
}
\]

したがって万能相互作用配列そのものは

\[
\boxed{
S(\psi_n)=e^B e^A
}
\]

であり、名前付き配列係数は

\[
\boxed{
S_{A,B}(\psi_n)
=\left[e^Be^A\right]_{A,B}
}
\]

と定義する。

`e^B e^A` の積順序を逆転してはならない。

---

# 24. matrix exponential action

巨大 dense `S` を必ずしも全要素 materialize する必要はない。

必要なのは

\[
y=e^{A}x
\]

の高精度 action である。

Krylov / polynomial / rational approximation のいずれを使うかは `exp_method_id` 状態で固定する。

`ε_exp,m_exp,max` も状態。

指数作用の作業 vector、Krylov basis、Hessenberg matrix は step 内一時量であり、次 step へ持ち越さない。

ただし再現性ログとして、各指数作用の

- used dimension
- estimated error
- residual

を readout log に保存してよい。

---

# 25. 実装監査 identity

毎 CF4 stage で

\[
\boxed{
\|K(\psi)\psi-F_{\rm direct}(\psi)\|_\infty
<\epsilon_K
}
\]

を検査する。

`F_direct` は同じ Einstein–Maxwell–GH–gauge 方程式を独立な RHS evaluator として評価したもの。

この監査は数値一致確認であり、別の動力学経路として次状態生成には使用しない。

---

# 26. 読み出し ― apparent horizons

readout は状態へ帰還しない。

各時刻面の apparent horizon `H_a,H_b` を一時 iterative solver で探索する。

horizon finder の反復値は永続状態ではない。

初期化 target の `a,b` ラベルと現在 horizon component の対応は連続 world tube の readout 上の対応として管理し、`S` へ戻さない。

---

# 27. 電荷 readout

各 horizon で

\[
\boxed{
Q_H
=\frac{1}{4\pi}
\oint_H E^i s_i\,dA
}
\]

を読む。

従って

\[
Q_a^{\rm read}(\psi),
\qquad
Q_b^{\rm read}(\psi)
\]

を得る。

初期値 `Q_a^(0),Q_b^(0)` は比較用恒等状態。

---

# 28. spin readout

quasi-local rotational vector `\varphi^i` を horizon geometry から求め、

\[
\boxed{
J_H[\varphi]
=\frac{1}{8\pi G}
\oint_H
K_{ij}\varphi^i s^j\,dA
}
\]

を quasi-local angular momentum とする。

generic dynamical horizon では rotational vector の選択による数値的曖昧さを readout uncertainty としてログする。

\[
\mathbf J_a^{\rm read}(\psi),
\qquad
\mathbf J_b^{\rm read}(\psi)
\]

を初期 target と比較する。

**current spin の独立 persistent register は持たない。**

---

# 29. horizon mass readout

horizon area `A_H` から

\[
M_{\rm irr}=\sqrt{\frac{A_H}{16\pi G^2}}
\]

（自然単位では `G=1`）を求める。

Kerr–Newman Christodoulou relation：

\[
\boxed{
M_H^2
=\left(M_{\rm irr}+\frac{Q_H^2}{4M_{\rm irr}}\right)^2
+\frac{J_H^2}{4M_{\rm irr}^2}
}
\]

から

\[
M_a^{\rm read},M_b^{\rm read}
\]

を求める。

---

# 30. 軌道 readout

apparent horizon の coordinate centroid を

\[
\mathbf x_a(\tau),\mathbf x_b(\tau)
\]

とする。

coordinate separation：

\[
\boxed{
r_{\rm coord}=|\mathbf x_a-\mathbf x_b|
}
\]

relative phase：

\[
\boxed{
\phi_{\rm orb}
=\operatorname{atan2}(y_a-y_b,x_a-x_b)
}
\]

orbital frequency：

\[
\boxed{
\Omega_{\rm orb}=d\phi_{\rm orb}/d\tau
}
\]

を readout series から求める。

coordinate dependence を監査するため、同時に

- proper separation on slice
- GW phase frequency

も保存する。

readout の時間差分は observer の後処理であり、`S` に帰還しない。

---

# 31. gravitational-wave readout

extraction 2-surface で null tetrad を構成し、

\[
\boxed{
\Psi_4
=-C_{abcd}n^a\bar m^b n^c\bar m^d
}
\]

を読む。

spherical harmonic mode：

\[
\Psi_4^{\ell m}(\tau,r_{\rm ex})
\]

を保存。

instantaneous outgoing GW power は `Ψ4` の time integral を observer side で行って評価する。

積分履歴は readout log であり、生成状態へ戻さない。

---

# 32. electromagnetic-wave readout

Poynting vector：

\[
\boxed{
S^i_{\rm EM}
=\frac{1}{4\pi}\epsilon^{ijk}E_jB_k
}
\]

extraction sphere の outward normal `s_i` に対して

\[
\boxed{
P_{\rm EM}^{\rm out}
=\oint S^i_{\rm EM}s_i\,dA
}
\]

を読む。

必要なら EM Newman–Penrose scalar `Φ_2` も保存する。

---

# 33. horizon absorption readout

各 horizon について

\[
\Delta M_H(\tau),
\quad
\Delta Q_H(\tau),
\quad
\Delta J_H(\tau)
\]

を readout series として保持する。

外向き GW/EM flux と horizon 増分をエネルギー・角運動量収支として比較する。

---

# 34. constraint readout

毎 readout step で少なくとも

\[
\|\mathcal C_a\|,
\quad
\|\mathcal C_{iab}\|,
\quad
\|\mathcal H\|,
\quad
\|\mathcal M_i\|,
\quad
\|D_iE^i\|,
\quad
\|D_iB^i\|
\]

を保存する。

これらは correction へ使用しない。

constraint growth は実験結果である。

---

# 35. finite-cutoff 実験系列

同一 initial physical state を、異なる cutoff だけで初期化する。

例：

\[
N_0\in\{N,2N,4N,8N\}.
\]

各走行で同じ `K(ψ)` 生成規則・同じ CF4・同じ readout を使用する。

比較対象：

- initial constraint residual
- `r(τ),φ(τ),Ω(τ)`
- `M_a,M_b`
- `Q_a,Q_b`
- `J_a,J_b`
- `Ψ4`
- EM flux
- constraint norms

これにより有限次数で収束するかを数値的に判定する。

---

# 36. baseline 実験行列

最初の physical baseline は spinless：

\[
\mathbf J_a^{(0)}=\mathbf J_b^{(0)}=0.
\]

regular KN gate を通る4条件：

1. n1 attractive
2. n1 repulsive
3. n3 attractive
4. n3 repulsive

n4 attractive は `INITIALIZATION_REJECTED_SUPEREXTREMAL` を記録する。

次の系列として spin を変更：

\[
\mathbf J_a^{(0)},\mathbf J_b^{(0)}
\]

だけを変え、同じ initializer・同じ `S(ψ)` を用いる。

---

# 37. forbidden implementation

以下は禁止する。

1. case_id による Einstein/Maxwell/spin/RN/KT 等の dynamics 分岐。
2. `Q,M,J,C,D` を time evolution 中に外部引数として渡す。
3. `a,b` の current spin を別レジスタで進化させる。
4. 外部 adaptive time-step history。
5. 外部 mesh controller history。
6. 前 step の RK stage を保持する。
7. readout した orbit/spin/flux を `S` へ戻す。
8. 初期化失敗を係数修正で救済する。
9. n4 superextremal を勝手に別モデルへ写す。
10. CF4 の低次代用品を本実験へ混在させる。
11. final exponentials `e^B e^A` の順序変更。
12. `K` の local nonlinear source の因数分解方法を case ごとに変える。

---

# 38. pre-execution static audit

実行前にコードを静的監査し、以下を全 PASS した場合だけ数値走行を許可する。

### A1 state closure
次状態生成が `ψ` 以外の mutable/persistent 情報を読まない。

### A2 single K generator
`K(ψ)` が一つだけ。

### A3 single S generator
CF4 実装が一つだけ。

### A4 named-index coverage
第22節の非零 block 以外が存在しない。

### A5 initial-state coupling prohibition
初期化専用 state column から field row への `K` が全ゼロ。

### A6 readout non-feedback
readout 関数の戻り値を transition path が参照しない。

### A7 no persistent work state
CF4 stage、matrix exponential work、horizon finder work が step 境界を越えない。

### A8 all configuration stateized
solver tolerance, cutoff, basis, domain map, damping, gauge, `Δτ` が external config のまま残っていない。

### A9 RHS identity
stage sample に対し `K(ψ)ψ=F_direct(ψ)`。

### A10 initializer constraints
Einstein + Maxwell constraints と target horizon parameters が tolerance 内。

### A11 reproducibility
state schema、initial state dump、hash、source、compiler/library version を保存。

---

# 39. run-time invariants/audit

毎 step または指定 stride で以下を確認する。

- identity state drift = 0
- `ONE=1`
- `D_tau` drift = 0（v1.0）
- `N_cut` drift = 0（v1.0）
- initial parameter states drift = 0
- finite `ψ`
- metric signature `(-,+,+,+)` 維持
- `det γ > 0`
- horizon excision boundary outflow condition
- exponential action residual <= `ε_exp`
- constraint norms

違反した場合は補正せず、その時点の `ψ` と diagnostics を保存して停止する。

---

# 40. 保存すべき raw output

各走行で最低限以下を保存する。

1. `STATE_SCHEMA.json` ― 全名前付き状態変数と整数 index 対応
2. `INITIAL_STATE.npy/json` ― `ψ_0`
3. `INITIALIZATION_AUDIT.json`
4. `K_BLOCK_SCHEMA.json`
5. `STATIC_AUDIT.json`
6. `RUN_CONFIG_FROM_STATE.json` ― `ψ` から readout した設定値
7. `READOUT_TIMESERIES.csv/parquet`
8. `CONSTRAINT_TIMESERIES.csv/parquet`
9. `HORIZON_TIMESERIES.csv/parquet`
10. `GW_MODES.*`
11. `EM_FLUX.*`
12. `FINAL_STATE.*`
13. `RUN_HASHES.sha256`
14. source code + exact version

---

# 41. 仕様上の状態分類

## 41.1 identity state

- 初期物理条件
- 方程式係数
- 数値精度設定
- cutoff
- basis/domain/map
- `D_tau`

## 41.2 dynamic state

- `T_tau`
- `g_ab`
- `Pi_ab`
- `Phi_iab`
- `H_a`
- `Theta_a`
- `E^i`
- `B^i`

## 41.3 temporary derived quantities

- inverse metric
- lapse/shift
- connections
- constraints
- stress-energy
- differential/multiplication matrices
- `K^[1..4]`
- CF4 stage states
- exponential internal work

## 41.4 readout only

- current mass
- current charge
- current spin
- horizon geometry
- orbit
- GW
- EM flux
- absorption
- constraint diagnostics

---

# 42. 最終定義

本仕様における万能相互作用とは、単一の固定数値行列ではない。

\[
\boxed{
\text{現在の全状態 }\psi
\longmapsto
K(\psi)
\longmapsto
S(\psi)
}
\]

という**唯一の状態関数**である。

物理ケース差は初期 `ψ_0` にのみ存在する。

最終 transition：

\[
\boxed{
\psi_{n+1}=S(\psi_n)\psi_n
}
\]

with

\[
\boxed{
S(\psi_n)=e^{B(\psi_n)}e^{A(\psi_n)}
}
\]

and `A,B` は第23節の4次 CF4 stage から生成する。

これ以外の dynamics path は本実験には存在しない。

---

# 43. 参考文献・基準資料

1. L. Lindblom, M. A. Scheel, L. E. Kidder, R. Owen, O. Rinne, “A New Generalized Harmonic Evolution System,” Class. Quantum Grav. 23, S447 (2006), arXiv:gr-qc/0512093.
2. L. Lindblom, B. Szilágyi, “An Improved Gauge Driver for the Generalized Harmonic Einstein System,” Phys. Rev. D 80, 084019 (2009), arXiv:0904.4873.
3. S. Mukherjee, N. K. Johnson-McDaniel, W. Tichy, S. L. Liebling, “Conformally curved initial data for charged, spinning black hole binaries on arbitrary orbits,” Phys. Rev. D 105 (2022), arXiv:2202.12133.
4. O. Dreyer, B. Krishnan, E. Schnetter, D. Shoemaker, “Introduction to Isolated Horizons in Numerical Relativity,” Phys. Rev. D 67, 024018 (2003), arXiv:gr-qc/0206008.
5. E. Celledoni, A. Marthinsen, B. Owren, “Commutator-free Lie group methods,” Future Generation Computer Systems 19 (2003) 341–352.
6. SpECTRE documentation, `gh::TimeDerivative`, current reference implementation of first-order GH RHS and matter source structure.
7. 木原範昭, 第六思考実験 v1.1 および補遺：荷電準円二体系の状態閉包と `G=q0=c=1, M/q0=20/3` 初期化同値検証。
8. 木原範昭, 第七思考実験 v1.1：`aa,bb` 自己相互作用初期化と旧条件の継承。

---

# 44. v1.0 で未だ「結果」であって仕様で決めないもの

以下は仕様段階で答えを仮定しない。

1. 有限 `N_cut` で長時間収束するか。
2. constraint drift がどの程度出るか。
3. Maxwell divergence cleaning state が後続版で必要になるか。
4. `Δτ` を動的 state として変更する必要があるか。
5. grid/map state を動的に進化させる必要があるか。
6. spin precession、charge transfer、horizon absorption がどのように現れるか。
7. GW/EM 放射が old LO model とどの領域で一致・不一致になるか。
8. 高次 mode が発散するか、飽和するか。

これらはすべて**万能相互作用を実際に走らせて得る実験結果**とする。

---

# 45. 完了条件

この仕様の実装版が「完成」と呼べるのは、最低限以下を満たした場合だけである。

\[
\boxed{\text{initializer PASS}}
\]

\[
\boxed{\text{static audit PASS}}
\]

\[
\boxed{K(\psi)\psi=F_{\rm direct}(\psi)}
\]

\[
\boxed{\psi_{n+1}=S(\psi_n)\psi_n}
\]

\[
\boxed{\text{same }S\text{ for every charge/spin case}}
\]

\[
\boxed{\text{orbit, }M,Q,J,GW,EM\text{ all read from }\psi}
\]

が確認された時点である。

