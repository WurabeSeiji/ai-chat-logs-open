# 万能相互作用 Einstein–Maxwell 初期化仕様 v3.1
## — 乗法正式系・加法検算系へ同一の完全初期状態 Ψ₀ を供給する共通初期化仕様（水素型二体 a,b） —

**著者:** 木原範昭 (Noriaki Kihara)  
**仕様整理:** ChatGPT（v3.0 まで）、Claude（v3.1 改定）  
**版:** v3.1  
**日付:** 2026-10-01  
**上位相互作用仕様:** `universal_interaction_einstein_maxwell_multiplicative_spec_ja_v3.2_20261001.md`  
**置換対象:** `universal_initialization_einstein_maxwell_spec_ja_v3.0_20260930.md`

---

# 0A. v3.1 改定（2026-10-01、木原範昭の指示による）

本改定の目的・原則・単位・`δτ`・注意書きは、上位相互作用仕様 v3.2 の §0A に従う。初期化に関わる変更は次のとおりである。

1. **単位**：\(G=q_0=c=1\)（\(q_0\) は電荷素量の 1/3）、電荷は Gauss 幾何単位（\(4\pi\epsilon_0=1\)）。質量を単位とする \(G=M=c=1\) と、それに基づく初期分離 \(d_0=P_0M_{tot}\) は誤りであり、撤回する。
2. **電荷の量子化**：電荷は \(q_0\) の整数倍 \(Q=n_q q_0\) だけを取る。第 1 実験の対象は水素型二体（a：電子型 \(Q_a=-3\)、b：陽子型 \(Q_b=+3\)）である。
3. **地平面前提の除去**：Einstein–Maxwell 方程式と Kerr–Newman の厳密式は変更しない。ただし、点の発散を避けるために入れた正則性ゲート、excision、地平面境界条件、地平面での目標値照合は外す。背景の空間も時間もない系では「点」に意味がない。
4. **回転（スピン）**：電子も陽子もフェルミオン（スピン 1/2）なので、最初から回転が必要である。回転しない a,b には内部時計がない。初期スピンは \(J_a=J_b=\hbar/2\) とする（§7A）。
5. **a,b の量の定義**：a,b の質量・電荷・角運動量は、a（または b）だけを囲む閉曲面上の閉じた位相（Gauss 電束、慣性系の引きずりの位相、Komar 型積分）で定義する（§17A）。
6. **δτ**：恒等状態ではない。初期値は系の時計の 1 周期が \(N_{\rm tick}=31\) 刻みとなるように定める（§2.2）。

**注意（将来の誤読防止）**：初期スピン \(J=\hbar/2=9/(2\alpha)\)、\(N_{\rm tick}=31\)、1 刻み内の 4 段は、\(\alpha\) の値と有限位数根 \(m/n=23/124\) から予測値として与えたものである。本研究は最小の公理系が自己無撞着になることを条件と考えている。そのため循環のように見えるが、これは隠れたパラメータの導入ではない。任意に与えたパラメータが自己無撞着になるには特定の値でなければならないことを示そうとするものである（相互作用仕様 v3.2 §0A.5）。

---

# 0. 改定理由（v3.0 の記録）

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

# 2. 上位 interaction v3.2 との state contract

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

`kappa_tau` は初期化後の恒等 state である。`delta_tau,D_tau` は恒等 state ではなく、発展中は相互作用仕様 v3.2 §3.2.1 の規則で毎 step 定まる。

初期値は、系の時計の 1 周期がちょうど \(N_{\rm tick}=31\) 刻みになるように

\[
\boxed{
\delta\tau_0=\frac{2\pi}{N_{\rm tick}\,\Omega_0},
\qquad
D_{\tau,0}=e^{\kappa_\tau\delta\tau_0}
}
\]

とする。\(\Omega_0\) は §5.3 の初期の相対回転の角振動数である。\(N_{\rm tick}=31\) は恒等 state として保持する。

initializer は

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

interaction v3.2 は、採用した数値 Einstein–Maxwell solver の空間離散化・境界処理を `F_h` の一部としてそのまま使用する。

したがって initializer v3.0 は、interaction v2.0 のように独自の `D_i,V_*,P_*` を上位物理原則として要求しない。

代わりに次を要求する。

1. 初期場は、正式発展 solver が受理する **同一の field representation と grid representation** へ変換して格納する。
2. interpolation / projection / prolongation が必要なら、採用 solver の公式変換 routine を使う。
3. initializer が独自に作った hidden FD stencil / hidden RBF stencil を正式発展 field の定義へ混入しない。
4. evolution grid と initialization elliptic grid が異なる場合、その変換 algorithm・tolerance・version/hash を `Psi_initcfg` に恒等記録する。
5. 変換後に physical constraints を evolution representation 上で再評価する。

---

# 4. 全ケース共通の物理初期化法

二体 a,b の実験の物理初期化法は

\[
\boxed{
\text{superposed boosted Kerr–Newman free data}
\rightarrow
\text{coupled XCTS + electromagnetic constraints}
}
\]

とする。Kerr–Newman の厳密式は、任意の \(Q,J\) に対してそのまま用いる。\(|Q|/M>1\) や \(J/M^2>1\)（地平面を持たない領域）でも式は Einstein–Maxwell 方程式の厳密解であり、式を変更しない。

重要：これは `Q=0`, `J=0`, 同符号電荷, 逆符号電荷, 非対称電荷, 質量比の異なる二体の全ケースで同じである。

禁止：

- MP exact solution への case switch
- double-RN exact solution への case switch
- neutral case だけ別 puncture initializer へ切替
- opposite-charge case だけ別式へ切替
- solver failure を別 initializer で救済

**v3.1 撤回**：v3.0 の Kerr–Newman regularity gate

\[
\left(\frac{Q_A^{(0)}}{M_A^{(0)}}\right)^2
+
\left(\frac{|J_A^{(0)}|}{M_A^{(0)2}}\right)^2
\le 1
\]

と `INITIALIZATION_REJECTED_SUPEREXTREMAL` は撤回する。このゲートは a,b を地平面を持つブラックホールと仮定し、点の発散を避けるために置かれたものであり、方程式の要請ではない。電子型・陽子型の a,b（\(|Q|/M\gg1\)、\(J/M^2\gg1\)）はこのゲートを満たさないが、初期化の対象である。

---

# 5. N_body=0,1,2 の扱いを分離して矛盾を除く

`N_body`（v3.0 の `N_BH`。a,b をブラックホールと仮定しないため改名）は topology input であり、case-specific physics switch ではない。

## 5.1 N_body=0

Minkowski validation 専用。

- 種の和は空和。
- `d0,C0,Omega0,Omega_r,A` は **未使用**。
- free data は flat / zero EM。
- `u_0` は flat-space fixed point を作る。

## 5.2 N_body=1

single Kerr–Newman validation（電子型または陽子型の一体、回転あり）。

- `d0,C0,Omega_guess` は **未使用**。
- center は `X_1^(0)`。
- translational boost `v_1` と spin `J_1` は明示 input からのみ生成。
- excision・地平面境界は持たない。

## 5.3 N_body=2

正式な二体実験。

### 5.3.1 電荷と質量

\[
Q_a^{(0)}=-3,\qquad Q_b^{(0)}=+3\qquad(\text{単位 }q_0)
\]

質量は \(M_a^{(0)}\)（電子型）、\(M_b^{(0)}\)（陽子型）とし、次の二つの比で指定する。

- 質量比 \(M_b/M_a\)：実値 1836.15 または縮小値
- クーロン力と重力の比

\[
C_0=1-\lambda_a^{(0)}\lambda_b^{(0)},
\qquad
\lambda_A=\frac{Q_A}{M_A}
\]

を \(10^2\le C_0\le10^6\) の範囲で選ぶ。実際の水素は \(C_0\approx2.3\times10^{39}\) だが、重力による歪みが倍精度の丸め誤差に埋もれず読み出せるよう、この範囲で選ぶ。質量の選び方は重力の強さと計算時間に効くが、以下の比（軌道半径/回転半径 \(=2/\alpha\)、\(v/c=\alpha\)）には影響しない。

参考例（質量比 1836.15、\(C_0=10^4\)）：\(M_a\approx7.00\times10^{-4}\)、\(M_b\approx1.285\)。採用値は run manifest に exact 値で記録する。

### 5.3.2 初期分離と相対回転

\(P_0=50\) から \(d_0=P_0M_{tot}\) を作る v3.0 の規則は撤回する（質量を長さの単位にしており、\(G=q_0=c=1\) に反する）。

初期の軌道角運動量を、スピンと同じ単位で

\[
\boxed{L_0=n\hbar=2nJ,\qquad n=1}
\]

とし（\(J=\hbar/2\)、§7A）、円軌道の力の釣り合い（Gauss 幾何単位、\(G=1\)）

\[
\mu\,\Omega_0^2\,d_0=\frac{K_0}{d_0^2},
\qquad
K_0=M_aM_bC_0=M_aM_b-Q_aQ_b,
\qquad
\mu=\frac{M_aM_b}{M_{tot}}
\]

から

\[
\boxed{
d_0=\frac{L_0^2}{\mu K_0}
}
\]

\[
\boxed{
\Omega_0^2=\frac{M_{tot}C_0}{d_0^3}
}
\]

とする。クーロン力が支配的なとき \(K_0\approx9\) で、\(d_0\) は Bohr 半径 \(\hbar^2/(\mu e^2)\) に一致する。これは軌道の初期推定であり、発展方程式へ注入しない。

`C0<0` で `Omega0` が実数にならない場合は、準円軌道 baseline の domain 外として明示 FAIL し、絶対値化等をしない。

---

# 6. binary baseline coordinate state

重心を原点とする（v3.0 の等質量 baseline は撤回）：

\[
\mathbf X_a^{(0)}=\left(+\frac{M_b}{M_{tot}}d_0,0,0\right),
\qquad
\mathbf X_b^{(0)}=\left(-\frac{M_a}{M_{tot}}d_0,0,0\right)
\]

center of mass：

\[
\mathbf r_{CM}
=\frac{\sum_A M_A^{(0)}\mathbf X_A^{(0)}}{M_{tot}^{(0)}}
\]

baseline：

\[
\boxed{\Omega_0\text{ は §5.3.2 の値}},\qquad
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

各体 A の seed parameters：

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
\mathbf J_A^{seed}=\mathbf J_A^{(0)}
}
\]

（v3.0 の \(\mathbf J=\boldsymbol\chi M^2\) と \(|\chi|\le\sqrt{1-(Q/M)^2}\) の制約は撤回する。）

## 7A. 初期スピン（v3.1）

電子も陽子もフェルミオン（スピン 1/2）なので、最初から回転が必要である。回転しない a,b の場は静的で、内部時計として読める量がない。

大きさは両者で同じ

\[
\boxed{
|\mathbf J_a^{(0)}|=|\mathbf J_b^{(0)}|=\frac{\hbar}{2}
}
\]

とする。電荷 \(Q=3\)（単位 \(q_0\)）と \(\alpha=e^2/(\hbar c)\)（Gauss 規約）から

\[
\hbar=\frac{9}{\alpha}\approx1233.3,
\qquad
\frac{\hbar}{2}=\frac{9}{2\alpha}\approx616.7
\]

である。\(\alpha\) は CODATA 2022 の \(\alpha^{-1}=137.035999177\) を用い、exact な入力値として manifest に記録する。

向きは軌道の回転軸（z 軸）にそろえる。

- 基準 run：電子型と陽子型のスピンを反平行（全スピン 0、水素の超微細構造の基底状態に対応）
- 比較 run：平行（全スピン 1）
- 対照 run：\(\mathbf J_a=\mathbf J_b=0\)（内部時計のない場合）

a,b の内部時計は状態として追加せず、場から読み出す導出量とする（\(\omega_{\rm int}=M/J\)。Kerr–Newman の回転半径 \(J/M\)、磁気回転比 \(g=2\)）。読み出しは読出し仕様 v3.1 に定める。

**注意**：この初期値は \(\alpha\) から予測値として与えたものであり、隠れたパラメータではない（§0A の注意書き、相互作用仕様 v3.2 §0A.5）。

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

**v3.1**：a,b の位置に内側境界（地平面境界、excision 面）を置かない。v3.0 の quasi-equilibrium horizon boundary condition、apparent-horizon condition、inner lapse gauge、horizon charge の Neumann condition は撤回する。

a,b の種（Kerr–Newman の厳密式）が持つ特異的な部分は、厳密式のまま解析的に保持し、楕円型方程式で解くのは種からの補正場だけとする（標準の puncture 法と同じ扱い）。背景の空間も時間もない系では「点」に意味がないので、点の近くの値を理由に方程式や領域を変更しない。

電荷は種の厳密式が持つ値のまま保持し、Gauss 補正 \(\varphi\) は外側境界条件と正則性だけで決める。

---

# 14. excision surface（v3.1 で撤回）

v3.0 の excision surface \(\mathcal S_A^{exc}\)、apparent horizon の存在と包含の監査、`EXCISION_NOT_INSIDE_HORIZON` はすべて撤回する（§0A の 3）。

# 14A. a,b を囲む閉曲面（v3.1）

a,b の量を読むために、a（または b）だけを囲み、他方を含まない閉曲面

\[
\boxed{\Sigma_A\quad(A=a,b)}
\]

を `Psi_initcfg/Psi_grid` の固定設定として与える。\(\Sigma_A\) は境界条件ではなく、閉じた位相を読むための面である。面の選び方はゲージであり、選んだ面の定義（形・大きさ・中心の決め方）を恒等記録する。

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
M_{\Sigma,a}-M_a^{(0)},
M_{\Sigma,b}-M_b^{(0)},
\mathbf P_{ADM},
\mathbf J_{\Sigma,a}-\mathbf J_a^{(0)},
\mathbf J_{\Sigma,b}-\mathbf J_b^{(0)}
)
}
\]

ここで \(M_{\Sigma,A},\mathbf J_{\Sigma,A}\) は §14A の閉曲面 \(\Sigma_A\) 上の閉じた位相（Komar 型積分、慣性系の引きずり）として読む量である（v3.0 の地平面上の \(M_H,J_H\) を置き換える）。

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

charge matchingは constraint solveの結果として、閉曲面 \(\Sigma_A\) を貫く電束（Gauss）で監査する。

\[
|Q_{\Sigma,A}-Q_A^{(0)}|<\epsilon_{charge_match}
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

closed-surface targets（v3.0 の horizon targets を置き換える）：

\[
|M_{\Sigma,A}-M_A^{(0)}|<\epsilon_{mass_match}
\]

\[
|Q_{\Sigma,A}-Q_A^{(0)}|<\epsilon_{charge_match}
\]

\[
\|\mathbf J_{\Sigma,A}-\mathbf J_A^{(0)}\|<\epsilon_{spin_match}
\]

ADM momentum：

\[
\|\mathbf P_{ADM}\|<\epsilon_{ADM_P}
\]

これらがすべてPASSするまで `Psi_0` を返さない。これは初期データが指定した値と拘束式を満たしているかの確認であり、物理結果についての予想ではない。

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

interaction v3.2 の common RHS evaluator が採用する gauge prescriptionと**同一**の target `F_a` を使う。

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
9. closed-surface M/Q/J targets（§14A の \(\Sigma_A\) 上）
10. ADM momentum
11. all identity states finite

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

- `N_body`
- 単位規約 `G=q0=c=1`、電荷規約 `GAUSSIAN_GEOMETRIZED`（4πε₀=1）
- initial masses, charges（\(q_0\) の整数倍 \(n_q\) として exact 記録）, spins（\(J=\hbar/2\)、使用した \(\alpha\) の exact 値）
- initial positions
- `L0, n, d0, Omega0, E0, phi0, chi0, tau0`
- 質量比、\(C_0\)、charge ratios and audit labels

## time

- `delta_tau`（初期値。発展中は状態から定まる）
- `D_tau`（同上）
- `kappa_tau`
- `N_tick = 31`

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
- closed-surface \(\Sigma_A\) specification（§14A）
- interpolation/projection method ID

隠れ default を正式runで使用しない。

---

# 24. FAIL codes

最低限：

```text
INITIALIZATION_INVALID_BINARY_OMEGA
INITIALIZATION_SEED_SUPERLUMINAL
XCTS_NEWTON_DIVERGED
XCTS_LINEAR_SOLVE_FAILED
XCTS_LINE_SEARCH_FAILED
XCTS_CONSTRAINT_FAIL
OUTER_CONTROL_DIVERGED
CLOSED_SURFACE_MASS_MATCH_FAIL
CLOSED_SURFACE_CHARGE_MATCH_FAIL
CLOSED_SURFACE_SPIN_MATCH_FAIL
ADM_MOMENTUM_MATCH_FAIL
GH_CONVERSION_FAIL
GH_CONSTRAINT_INIT_FAIL
MAXWELL_CONSTRAINT_INIT_FAIL
EVOLUTION_REPRESENTATION_CONVERSION_FAIL
TIME_STATE_INIT_FAIL
INITIAL_STATE_COPY_MISMATCH
```

FAIL時に target physical values、time state、solver physicsを自動変更して続行してはならない。

v3.1 では `INITIALIZATION_REJECTED_SUPEREXTREMAL` と `EXCISION_NOT_INSIDE_HORIZON` を撤回した（§4、§14）。

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

v3.1では次を撤回する。

1. 単位 \(G=M=c=1\) 系の初期分離 \(d_0=P_0M_{tot}\) と等質量 baseline。
2. Kerr–Newman regularity gate と spin 制約 \(|\chi|\le\sqrt{1-(Q/M)^2}\)。
3. excision surface、地平面境界条件、地平面上の目標値照合。
4. `delta_tau,D_tau` を恒等 state とすること。

---

# 26. v3.1の完成条件

初期化仕様が完成しているとは、承認済み physical inputs と numerical configuration に対して、追加判断なしに

\[
\boxed{\Psi_0\quad\text{または明示FAIL}}
\]

を一意に返し、その `Psi_0` が interaction v3.2 の

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

- [ ] interaction v3.2 state schemaと一致
- [ ] ONEあり
- [ ] T_tau, delta_tau, D_tau, kappa_tauあり
- [ ] time consistency式あり
- [ ] 二体だけにd0/Omega ruleを適用
- [ ] N_body=0/1で二体の式を参照しない
- [ ] 単位 G=q0=c=1、Gauss 規約、電荷が q0 の整数倍
- [ ] 初期スピン J=ħ/2 と向き（反平行/平行/対照 0）が確定
- [ ] delta_tau の初期値が N_tick=31 から定まる
- [ ] inner Newton中にcontrol parameterを更新しない
- [ ] outer control loopがinner solveの外側
- [ ] excision・地平面境界条件・正則性ゲートを使っていない
- [ ] 閉曲面 Σ_A の定義が確定
- [ ] one initializer output only
- [ ] multiplication/addition branches are exact copies of same Psi0
- [ ] evolution representation上でconstraintsを再監査
- [ ] hidden fallbackなし
- [ ] hidden dtなし

一つでも未定義なら正式実験コードを書いてはならない。

---

# 28. 最終定義

v3.1初期化は

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
