# Einstein–Maxwell 万能相互作用：水素型二体 第1実験・読出・図化仕様 v3.1
## ― 背景時空も粒子も先に置かない二体 a,b の発展、系の時計 δτ と a,b の内部時計、放射・吸収収支、無損失記録 ―

**著者:** 木原範昭 (Noriaki Kihara)  
**仕様整理:** ChatGPT（v3.0 まで）、Claude（v3.1 改定）  
**仕様版:** v3.1  
**日付:** 2026-10-01  
**上位相互作用仕様:** `universal_interaction_einstein_maxwell_multiplicative_spec_ja_v3.2_20261001.md`
**上位初期化仕様:** `universal_initialization_einstein_maxwell_spec_ja_v3.1_20261001.md`
**置換対象:** `experiment_readout_visualization_spec_einstein_maxwell_ja_v3.0_20260930.md`

---

# 0A. v3.1 改定（2026-10-01、木原範昭の指示による）

本改定の目的・原則・単位・`δτ`・注意書きは、上位相互作用仕様 v3.2 の §0A に従う。読出しに関わる変更は次のとおりである。

1. **実験の立場**：v3.0 は「天体スケールの古典・連続系として扱う」「第 1 段階では素粒子状態が得られたと解釈しない」「continuous classical first」「量子化仮説を持ち込まない」としていた。これは誤りである。本研究の原則は、外から素粒子状態を**与えてはならない**ことであり、結果として素粒子状態に類似したものが読み出された場合にそう解釈することを禁じるものではない。また本研究はスケール不変を原則とし、天体スケールで成り立つ厳密な方程式が原子スケール（Planck 長）でも成り立つという仮説に立つ。
2. **対象**：第 1 実験の対象は水素型二体（a：電子型 \(Q_a=-3\)、b：陽子型 \(Q_b=+3\)、単位 \(G=q_0=c=1\)、Gauss 規約）である。電子・陽電子型の対は対消滅するので対象にしない。安定して存在できる最小の構成は水素型である。
3. **m/n の読み出しの解禁**：v3.0 §3・§30・§43 が禁止していた \(m/n\) の読み出し、刻みと周期の比の読み出し、有限位数根 \(\cos^2(\pi m/n)\) との比較は、第 1 実験から読み出しに含める。
4. **δτ**：`δτ` は数値刻みではなく系の唯一の時計の刻みであり、系の時計の 1 周期が \(N_{\rm tick}=31\) 刻みとなるよう状態から定まる。v3.0 の \(h,h/2,h/4\) による \(h\to0\) 収束の要求は撤回する。
5. **地平面の除去**：地平面での吸収・地平面上の \(M,Q,J\) の読み出しは、a,b だけを囲む閉曲面 \(\Sigma_A\) 上の閉じた位相の読み出しに置き換える。
6. **内部時計とロック**：a,b の内部時計 \(\omega_{\rm int}=M/J\)、系の時計 \(\Omega_{\rm sys}\) とのロック比、\(\hbar_{\rm eff}=2J\)、\(\alpha_{\rm eff}\) を読み出しに加える。
7. **予想は成功条件にしない**：分類ラベル・準定常判定・ロックの有無などは、生の読み出しとして保存し、結果が出てから解釈する。実験の成功条件は §42 に定める。

**注意（将来の誤読防止）**：1 刻み内の 4 段、\(N_{\rm tick}=31\)、初期スピン \(J=\hbar/2=9/(2\alpha)\) は、\(\alpha\) の値と有限位数根 \(m/n=23/124\)（\(124=31\times4\)）から予測値として与えたものである。本研究は最小の公理系が自己無撞着になることを条件と考えている。そのため循環のように見えるが、これは隠れたパラメータの導入ではない。任意に与えたパラメータが自己無撞着になるには特定の値でなければならないことを示そうとするものである（相互作用仕様 v3.2 §0A.5）。

---

# 0. 本仕様の位置づけ

本仕様は、相互作用仕様 v3.2 で定義した正式な乗法更新

\[
\boxed{\Psi_{n+1}^{(\times)}=S(\Psi_n^{(\times)})\Psi_n^{(\times)}}
\]

を実際に走らせたときの**第1実験の読出・図化・検算仕様**を定める。

同時に、同一の初期状態、同一の Einstein–Maxwell 離散 RHS、同一の空間離散化・境界処理・時間積分係数を用いる加法参照系

\[
\Psi_{n+1}^{(+)}=\Psi_n^{(+)}+\Delta\Psi_n^{(+)}
\]

を**実装検算専用**として並走させる。これは実装が意図どおりの計算をしているかの確認であり、実験の問いではない。

正本となる状態列は乗法系だけである。

\[
\boxed{\Psi^{(\times)}\text{ が正式状態列}}
\]

加法系は物理判定・主図・論文用の正式値へ採用しない。

初期化仕様 v3.1 が生成する唯一の初期状態 \(\Psi_0\) を bitwise copy して

\[
\boxed{\Psi_0^{(\times)}=\Psi_0^{(+)}}
\]

から両経路を開始する。

raw data、signed residual、exact parameter、全時間系列を保存し、狭い谷を見落とさない記録とする。

---

# 1. 第1実験の目的と読み出し

目的は、背景時空も粒子も先に置かず、厳密な Einstein–Maxwell 方程式を乗法更新で発展させたとき、水素型二体 a,b の状態の変化から、時空と粒子様の動力学がどのように読み出されるかを見ることである（相互作用仕様 v3.2 §0A.1）。

物理上の中心readoutは、**乗法正式系からのみ**

\[
\boxed{
T\mapsto
\left(
\delta\tau(T),\ 
\Omega_{\rm sys}(T),\ 
r_{\rm orb}(T),\ 
\Delta F_{\rm GW}(T),\ 
\Delta F_{\rm EM}(T),\ 
J_a(T),\ J_b(T)
\right)
}
\]

を読む。

加法系についても同一readout関数 \(\mathcal R\) を適用するが、目的は

\[
\Delta\mathcal R=\mathcal R(\Psi^{(\times)})-\mathcal R(\Psi^{(+)})
\]

の実装検算だけである。

結果を後から整理するための記述用ラベルとして、次を用いる。これは予想でも成功条件でもない。

### A. 継続縮退
放射・吸収を含むエネルギー散逸が続き、\(r_{\rm orb}(T)\) の減少が止まらず、`δτ` が縮み続ける。

### B. 有限半径準定常
初期過渡後に \(dr_{\rm orb}/dT\to0\)、`δτ` がほぼ一定になる。

### C. 準定常へ向かうが弱い放射が残る
有限半径へ非常にゆっくり近づきながら、GW/EM の弱い flux が長時間残る。

### D. 膨張・散乱
\(r_{\rm orb}(T)\) が増加する。

### E. 非定常・振動状態
\(r_{\rm orb}\)、`δτ`、GW flux、EM flux が有限振幅で長時間振動する。

どのラベルに当たるかは結果として記録する。事前に特定の結論を仮定しない。

---

# 2. スケールに対する実験上の立場

本研究はスケール不変を原則とする。比だけが物理であり、絶対的なサイズには意味がない。主要読出しは

\[
\frac{Q}{M},\quad
\frac{r}{J/M},\quad
\frac{\Omega_{\rm sys}}{\omega_{\rm int}},\quad
\frac{\delta\tau}{P_{\rm sys}},\quad
\frac{v}{c}
\]

などの無次元比とする。

天体スケールで成り立つ厳密な方程式が原子スケールでも成り立つという仮説に立つので、結果を天体系・原子系のどちらとして解釈することも禁止しない。raw dimensionless data を残し、後に別スケールへ写像できる形式とする。

---

# 3. 第1実験で行わないこと

以下は禁止する。

1. smooth / interpolation により raw extrema を消すこと
2. 特定の結果を期待した post-selection（都合のよい run だけを残すこと）
3. 読み出した値を使って状態を補正すること
4. 既知の谷位置に合わせて、実行後に run 条件を差し替えること

\(m/n\) の読み出し、刻みと周期の比の読み出し、有限位数根 \(\cos^2(\pi m/n)\) との比較は禁止しない（v3.0 の禁止は撤回）。ただし比較は読み出し後の解釈であり、発展へ帰還させない。

---

# 4. ただし最初から保存する「将来用」情報

各 run の manifest には最低限次を保存する。

\[
M_a^{(0)},M_b^{(0)},
Q_a^{(0)},Q_b^{(0)},
\mathbf J_a^{(0)},\mathbf J_b^{(0)}
\]

\[
r_0,\quad
\phi_0,\quad
\Omega_0,\quad
e_0
\]

\[
\Delta\tau,\quad
\kappa_\tau,\quad
N_{\rm cut}
\]

\[
\gamma_0,\gamma_1,\gamma_2,
\mu_H,\eta_H,\mu_L,\mu_S,p_H
\]

\[
\delta\tau,\quad D_\tau,\quad\kappa_\tau,\quad T_{\tau,0}
\]

\[
\epsilon_{\rm equiv},\quad
\epsilon_{\rm readout-equiv},\quad
\epsilon_{\rm time-state}
\]

```text
time_integrator_id = RK4_CLASSICAL
multiplicative_update_spec = v3.2
unit_convention = G=q0=c=1
charge_convention = GAUSSIAN_GEOMETRIZED (4*pi*eps0=1)
N_tick = 31
alpha_input = 1/137.035999177 (exact 入力値)
J_input = hbar/2 = 9/(2*alpha)
```

RK4 tableau の \(a_{rs},b_r,c_r\) も恒等状態として保存する。

\[
\{R_{{\rm ext},k}\}
\]

および全 initializer / gauge / basis / numerical tolerance。

値は可能な限り、

- exact integer
- exact rational
- decimal representation

の両方を残す。

例：

```text
Q_a_exact = -3   (単位 q0)
Q_b_exact = +3
mass_ratio_exact = 1836.15267343
```

有理比・整数比としての解析は禁止しない（v3.0 の禁止は撤回）。

---

# 5. 状態変数として追加・確認すべき実験設定

親仕様の原則

\[
\boxed{\text{初期化するものは原則すべて状態}}
\]

に従い、本実験で初期設定する値も唯一の \(\Psi_0\) に恒等状態として保持する。乗法系と加法系で別設定を作らない。

少なくとも以下を確認・追加する。

## 5.1 終了条件設定

\[
T_{\max}
\]

\[
N_{\rm step,max}
\]

\[
r_{\rm stop}
\]

\[
N_{\rm period,min}\ge2
\]

これらは実験runの停止条件であり、初期化時に決める。v3.0 の \(r_{\rm common-horizon}\)（地平面の合体）は撤回し、a,b の閉曲面 \(\Sigma_a,\Sigma_b\) が互いを分離して囲めなくなった時点を \(r_{\rm stop}\) とする。

ロック位数 \(n\) を読むには、少なくとも 2 周期分の記録が必要である（有限位数共鳴論文：位数 \(n\) の根を読むのに \(2n\) の標本）。系の時計の 1 周期は 31 刻み（段で 124）なので、最低 2 周期（62 刻み、段で 248）以上を記録する。

## 5.2 読み出し設定

\[
N_{\rm readout\_stride}
\]

\[
N_{\rm checkpoint\_stride}
\]

\[
R_{{\rm ext},1},R_{{\rm ext},2},\ldots,R_{{\rm ext},N_R}
\]

\[
N_R
\]

## 5.3 plateau 判定設定

plateau 判定値は物理発展へ帰還させないが、再現性のため初期化値として保持する。

\[
W_{\rm avg}
\]

\[
N_{\rm orbit,judge}
\]

\[
C_{\rm tol,r},\quad
C_{\rm tol,flux}
\]

これらは結果を記述するためのラベル付けの設定であり、実験の成功条件ではない。

---

# 6. 第1実験の初期物理条件

電子も陽子もフェルミオン（スピン 1/2）なので、最初から回転が必要である（初期化仕様 v3.1 §7A）。

\[
\boxed{|\mathbf J_a^{(0)}|=|\mathbf J_b^{(0)}|=\frac{\hbar}{2}=\frac{9}{2\alpha}}
\]

向きは軌道の回転軸にそろえる。v3.0 の \(\mathbf J_a=\mathbf J_b=0\) の基準は撤回し、内部時計のない対照 run としてだけ残す。

spin は独立レジスタとして発展させず、閉曲面 \(\Sigma_A\) 上の閉じた位相（慣性系の引きずり、磁気双極子）から読み出す。

質量比は実値 1836.15 を基準とし、縮小値の run を追加してよい。v3.0 の等質量基準は撤回する。

---

# 7. 電荷条件と run の組

単位は \(G=q_0=c=1\)（\(q_0\) は電荷素量の 1/3）、電荷は Gauss 規約である。電荷は \(q_0\) の整数倍だけを取る。v3.0 の連続的な \(q=Q/M=0.3,0.6,0.9,0.95,0.97,0.99\)、極限 \(q=1\) の Majumdar–Papapetrou anchor は、この量子化と地平面前提に反するので撤回する。

## Run H0：水素型・スピン反平行（基準）

\[
Q_a=-3,\quad Q_b=+3,\quad
\mathbf J_a=-\mathbf J_b,\quad |\mathbf J|=\frac{\hbar}{2}
\]

質量比 1836.15、\(C_0\) は初期化仕様 v3.1 §5.3.1 の範囲。

## Run H1：水素型・スピン平行（比較）

\[
\mathbf J_a=\mathbf J_b,\quad |\mathbf J|=\frac{\hbar}{2}
\]

H0 との差として、スピン–スピン結合（水素の超微細構造に対応）の効果を読む。

## Run H2：回転なし（対照）

\[
\mathbf J_a=\mathbf J_b=0
\]

内部時計のない場合として記録する。

## Run H3：質量比・\(C_0\) の変更

H0 と同じ電荷・スピンで、質量比（縮小値）と \(C_0\)（\(10^2\)〜\(10^6\)）を変えた run。比 \(r/(J/M)\)、\(v/c\)、ロック比が質量の選び方に依存するかを読む。

---

# 8. 初期軌道条件

初期軌道は初期化仕様 v3.1 §5.3.2 に従う。軌道角運動量を \(L_0=n\hbar=2nJ\)（\(n=1\)）とし、

\[
d_0=\frac{L_0^2}{\mu K_0},\qquad
\Omega_0^2=\frac{M_{tot}C_0}{d_0^3}
\]

である。初期化器は \(d_0\)、\(\Omega_0\)、\(L_0\)、\(n\) を必ず report する。

v3.0 の旧 \(P_0=50\) から \(d_0=P_0M_{tot}\) を作る規則は撤回する（質量を長さの単位にしており、\(G=q_0=c=1\) に反する）。表示用の軌道半径は重心からの距離

\[
r_{{\rm orb},A}=|\mathbf x_A-\mathbf r_{CM}|
\]

とし、分離 \(d=|\mathbf x_a-\mathbf x_b|\) と併せて保存する。

---

# 9. \(\delta\tau\) の扱い

`δτ` は数値計算の刻み幅ではなく、系の唯一の時計の 1 刻みである（相互作用仕様 v3.2 §3.2.1）。

\[
\delta\tau_n=\frac{P_{\rm sys}(\Psi_n)}{N_{\rm tick}},
\qquad
N_{\rm tick}=31
\]

系の時計の 1 周期がちょうど 31 刻み（1 刻み内の 4 段で数えて 124 段）になるよう状態から定まる。放射によって系の時計が変われば `δτ` も変わる。

v3.0 の「同一の物理初期条件について \(h,h/2,h/4\) の 3 系列を走らせ、\(h\to0\) で収束することを必須とする」は撤回する。系の時間は状態から生まれるので、`δτ` より細かい時間は存在しない。`δτ` を細かくして連続極限を見ることは、背景の連続時間を前提にすることになる。与えた \(N_{\rm tick}\) で得た結果がそのまま結果である。

`δτ` の系列は中心の読み出しの一つとして全 step 保存する。

**注意**：\(N_{\rm tick}=31\) と 4 段は、\(\alpha\) に対応する有限位数根 \(m/n=23/124\) から予測値として与えたものであり、隠れたパラメータではない（§0A の注意書き）。

---

# 10. 空間 cutoff / mode order の収束確認

同じ物理runを、

\[
N_{\rm cut}^{(L)}
<
N_{\rm cut}^{(M)}
<
N_{\rm cut}^{(H)}
\]

の3水準で比較する。

特に、

\[
r_{\rm orb}(T)
\]

\[
\Delta F_{\rm GW}(T)
\]

\[
\Delta F_{\rm EM}(T)
\]

の差を保存する。

深い谷状のstructureが後に現れたとき、物理と spectral cutoff artifact を区別する材料として、最初からこのデータを残す。これは記録であり、実験の成功条件ではない。空間の格子はゲージのラベルであって背景ではない。

---

# 11. 時間 readout

時間の正本はstep番号や外部clockではなく、各状態が持つ

\[
T_\tau=e^{\kappa_\tau\tau}
\]

である。

正式な物理時間は乗法状態から

\[
\boxed{T^{(\times)}=\mathcal R_T(\Psi^{(\times)})=\frac{\log T_\tau^{(\times)}}{\kappa_\tau}}
\]

と読む。加法検算系にも同じ関数を適用し、

\[
T^{(+)}=\frac{\log T_\tau^{(+)}}{\kappa_\tau}
\]

を求める。

毎readoutで

\[
\boxed{\Delta T=T^{(\times)}-T^{(+)}}
\]

および

\[
\boxed{C_{\delta\tau}=\delta\tau-\frac{\log D_\tau}{\kappa_\tau}}
\]

を保存する。

必須raw項目：

```text
T_tau_mul
T_tau_add
T_read_mul
T_read_add
delta_tau
D_tau
kappa_tau
T_read_difference
time_state_consistency_error
```

以後、主図の `T` は常に `T_read_mul`、すなわち乗法正式系の時間を意味する。

## 11.1 系の時計で数えた時間（v3.1）

`δτ` は系の時計の刻みなので、`T` は刻みの積算

\[
T_n=T_0+\sum_{k<n}\delta\tau_k
\]

であり、step 番号の定数倍ではない。あわせて、系の時計の周期数

\[
N_{\rm period}(T)=\frac{\phi_{\rm sys}(T)}{2\pi}
\]

（\(\phi_{\rm sys}\) は a,b の相対回転の位相）を保存する。周期数を横軸にした図も作る。

## 11.2 a,b それぞれの固有時（v3.1）

一般相対論では、経過する固有時は場所ごと・世界線ごとに違う。v3.0 は系全体に一つの `T` しか持たず、時間の歪みを読み出していなかった。v3.1 では a,b それぞれの位置の lapse \(\alpha_{\rm lapse}\) と shift から、1 刻みあたりの固有時

\[
\delta\tau_A^{\rm proper}=\sqrt{-g_{ab}\,\dot x_A^a\dot x_A^b}\;\delta\tau
\]

を積算した \(\tau_A^{\rm proper}(T)\) と、赤方偏移 \(z_A=\delta\tau_A^{\rm proper}/\delta\tau\) を保存する。ここで \(\dot x_A\) は a（または b）の閉曲面の中心の座標速度であり、読み出し量である。発展へ帰還させない。

追加の必須raw項目：

```text
N_tick
Omega_sys
P_sys
phi_sys
N_period
tau_proper_a
tau_proper_b
redshift_a
redshift_b
```

# 11A. readout の正本と検算経路

任意のreadout関数を \(\mathcal R_A(\Psi)\) とする。

正式値は

\[
\boxed{R_A^{\rm official}=\mathcal R_A(\Psi^{(\times)})}
\]

である。加法系についても同じreadout実装を独立に適用し、

\[
R_A^{(+)}=\mathcal R_A(\Psi^{(+)})
\]

\[
\boxed{\Delta R_A=R_A^{(\times)}-R_A^{(+)}}
\]

を検算用に保存する。

加法readoutを正式物理値として採用したり、差を使って乗法状態を補正したりしてはならない。

---

# 12. 軌道 readout

最低限次を保存する。

\[
\mathbf x_a(T),\quad
\mathbf x_b(T)
\]

\[
d(T)
=
|\mathbf x_a-\mathbf x_b|
\]

equal-mass baseline では

\[
\boxed{
r_{\rm orb}(T)=\frac{d(T)}2
}
\]

を主図に使う。

さらに、

\[
\phi(T)
\]

\[
\Omega(T)=\frac{d\phi}{dT}
\]

\[
\dot r_{\rm orb}(T)
\]

を保存する。

座標依存性を監査するため、可能なら追加で

- proper separation
- 閉曲面 \(\Sigma_a,\Sigma_b\) の中心間の分離

も保存する。

主図は一つの定義に固定し、定義を途中で変えない。

---

# 13. GW readout

## 13.1 外向きGW

複数の抽出半径

\[
R_{{\rm ext},k}
\]

で

\[
\Psi_4(T,R_{{\rm ext},k},\theta,\phi)
\]

を読み出す。

球面調和モード

\[
\Psi_4^{\ell m}(T,R_{{\rm ext},k})
\]

を保存する。

主成分だけでなく、使用した全モードをraw保存する。

GW energy flux at large radius を

\[
F_{\rm GW}^{\rm out}(T;R_{{\rm ext},k})
\]

として計算する。

最終的には複数抽出半径から finite-radius error を監査する。

## 13.2 a,b への吸収（閉曲面上の閉じた位相、v3.1）

v3.0 の地平面での吸収は撤回する。a,b を地平面を持つブラックホールと仮定しないためである。

ゲージ理論は、電荷の内側を見ずに、電荷を囲む閉曲面・閉ループ上で位相が閉じるかだけで電荷を定義する（Gauss の法則、ホロノミー）。一般相対論の地平面も内側を隠すという同じ役割を果たしていた。本仕様では地平面の代わりに、a（または b）だけを囲む閉曲面 \(\Sigma_A\)（初期化仕様 v3.1 §14A）上の閉じた位相で a,b の量を定義し、**吸収をその時間変化**として定義する。

\[
\boxed{
F^{\rm abs}_{A}(T)=\frac{dE_{\Sigma,A}}{dT}
}
\]

\(E_{\Sigma,A}\) は \(\Sigma_A\) 上の Komar 型積分で読む a（または b）の質量・エネルギーである。GW と EM への分解は、\(\Sigma_A\) を貫く重力波のエネルギー流（\(\Psi_4\) 型、または Isaacson 型の有効エネルギー流）と Poynting flux で行う。

\[
F_{\rm GW}^{\rm abs}=F_{{\rm GW},a}^{\Sigma}+F_{{\rm GW},b}^{\Sigma}
\]

ここで `abs` は \(\Sigma_A\) を内向きに通る正のenergy flux magnitude とする。

閉曲面の選び方はゲージである。質量と角運動量は、系が時間変化していると囲む面によって値が変わりうる（電荷は真空中なら面によらない）。面の選び方への依存は、複数の面で読んで記録する。依存の大小についての予想は成功条件にしない。

---

# 14. EM readout

## 14.1 外向きEM

各 extraction sphere でPoynting fluxを積分し、

\[
F_{\rm EM}^{\rm out}(T;R_{{\rm ext},k})
\]

を得る。

可能なら electromagnetic Newman–Penrose radiative scalars も保存する。

## 14.2 a,b への吸収（閉曲面、v3.1）

各閉曲面 \(\Sigma_A\) を内向きに貫く Poynting flux を積分し、

\[
F_{{\rm EM},a}^{\Sigma}(T),
\qquad
F_{{\rm EM},b}^{\Sigma}(T)
\]

\[
\boxed{
F_{\rm EM}^{\rm abs}
=
F_{{\rm EM},a}^{\Sigma}
+
F_{{\rm EM},b}^{\Sigma}
}
\]

とする。

---

# 15. ユーザー指定の「放射–吸収差分」

本実験の主図では、指定通り次を図化する。

\[
\boxed{
\Delta F_{\rm GW}
=
F_{\rm GW}^{\rm out}
-
F_{\rm GW}^{\rm abs}
}
\]

\[
\boxed{
\Delta F_{\rm EM}
=
F_{\rm EM}^{\rm out}
-
F_{\rm EM}^{\rm abs}
}
\]

これらは**放射と吸収のbalance residual**として扱う。

正負を保存し、

\[
\Delta F=0
\]

を中心線として図示する。

---

# 16. 重要：軌道散逸量は差分とは別に保存する

エネルギー収支では、

- 無限遠へ出ていくenergy
- a,b の閉曲面 \(\Sigma_A\) を内向きに通るenergy（§13.2、§14.2）

は、どちらもorbital sectorから見ればenergy sink になり得る。

そこで正のmagnitude conventionを採る場合、

\[
\boxed{
D_{\rm GW}
=
F_{\rm GW}^{\rm out}
+
F_{\rm GW}^{\rm abs}
}
\]

\[
\boxed{
D_{\rm EM}
=
F_{\rm EM}^{\rm out}
+
F_{\rm EM}^{\rm abs}
}
\]

も必ず保存する。

総散逸readout：

\[
\boxed{
D_{\rm rad}
=
D_{\rm GW}+D_{\rm EM}
}
\]

したがって、

\[
\Delta F_{\rm GW}\simeq0
\]

だけを見て「GW散逸が消えた」と判断してはならない。

主図はユーザー指定の差分を使い、**エネルギー監査図では和も使う。**

---

# 17. cumulative energy readout

時間積分して、

\[
E_{\rm GW}^{\rm out}(T)
=
\int_0^T
F_{\rm GW}^{\rm out}(t)\,dt
\]

\[
E_{\rm GW}^{\rm abs}(T)
=
\int_0^T
F_{\rm GW}^{\rm abs}(t)\,dt
\]

\[
E_{\rm EM}^{\rm out}(T)
=
\int_0^T
F_{\rm EM}^{\rm out}(t)\,dt
\]

\[
E_{\rm EM}^{\rm abs}(T)
=
\int_0^T
F_{\rm EM}^{\rm abs}(t)\,dt
\]

を保存する。

さらに、

\[
E_{\rm rad}^{\rm total}
=
E_{\rm GW}^{\rm out}
+
E_{\rm GW}^{\rm abs}
+
E_{\rm EM}^{\rm out}
+
E_{\rm EM}^{\rm abs}
\]

を保存する。

---

# 18. a,b の量・内部時計・ロックの readout（v3.1）

v3.0 の horizon readout は撤回し、閉曲面 \(\Sigma_A\) 上の閉じた位相からの読み出しに置き換える。

## 18.1 a,b の量

各stepまたは指定readout cadenceで、

\[
M_{\Sigma,a}(T),\quad M_{\Sigma,b}(T)
\qquad(\text{Komar 型積分})
\]

\[
Q_{\Sigma,a}(T),\quad Q_{\Sigma,b}(T)
\qquad(\text{Gauss 電束、電磁ホロノミー})
\]

\[
\mathbf J_{\Sigma,a}(T),\quad \mathbf J_{\Sigma,b}(T)
\qquad(\text{慣性系の引きずりの位相、磁気双極子})
\]

を保存する。

## 18.2 内部時計（導出量）

a,b の内部時計は状態として持たず、場から読む導出量とする。

\[
\boxed{
\omega_{{\rm int},A}(T)=\frac{M_{\Sigma,A}}{|\mathbf J_{\Sigma,A}|}
}
\]

（Kerr–Newman の回転半径 \(J/M\) を光速で回る振動数。電子のパラメータでは Zitterbewegung 振動数 \(2Mc^2/\hbar\) に一致し、磁気回転比は \(g=2\) である。）

磁気双極子モーメント \(\mu_A\) と \(g_A=2M_A\mu_A/(Q_AJ_A)\) も保存する。

## 18.3 系の時計とのロック

\[
\boxed{
\rho_A(T)=\frac{\Omega_{\rm sys}(T)}{\omega_{{\rm int},A}(T)}
}
\]

を保存し、その有理近似 \(m/n\)（連分数展開の各段）を記録する。有限位数根 \(R_{n,m}=\cos^2(\pi m/n)\)（\(23/124\) など）との比較は読み出し後の解釈として行う。

## 18.4 作用の単位と結合定数の読み出し

\[
\hbar_{{\rm eff},A}(T)=2|\mathbf J_{\Sigma,A}(T)|
\]

\[
\alpha_{{\rm eff},A}(T)=\frac{Q_{\Sigma,A}^2}{\hbar_{{\rm eff},A}(T)}
\qquad(\text{Gauss 規約、}c=1)
\]

を保存する。Heaviside–Lorentz 規約の論文と比較するときは \(e_{\rm HL}=\sqrt{4\pi}\,e_{\rm Gauss}\) で換算する。

**注意**：初期スピン \(J=\hbar/2\) は \(\alpha\) から与えている（§0A の注意書き）。本系では導出は必ず循環する。ここで見るのは、与えた値と、動力学の結果として読み出した値が矛盾なく一周して閉じるか（自己無撞着性）である。

## 18.5 放射・吸収と閉曲面上の量の対応

\[
\Delta M_{\Sigma},\quad
\Delta J_{\Sigma},\quad
\Delta Q_{\Sigma}
\]

と吸収 flux の積分との対応を監査できるようにする。

---

# 19. constraint / numerical / equivalence readout

物理constraintは正式な乗法状態 \(\Psi^{(\times)}\) から最低限、

\[
\|{\cal C}_{\rm GH}\|,\qquad\|D_iE^i\|,\qquad\|D_iB^i\|
\]

を保存する。可能なら \(L_2\) と \(L_\infty\) の両方を保存する。

加法状態についても同じconstraint readoutを行い、readout差を検算として保存する。

旧 `Kpsi-F_direct` residual は使用しない。相互作用仕様 v3.2 の二経路同値監査を、実装が意図どおりの計算をしているかの**実装検算**として行う。これは実験の問いではない（相互作用仕様 v3.2 §0A.1、§23）。

\[
E_{stage,r}=\|u_r^{(\times)}-u_r^{(+)}\|_\infty
\]

\[
E_{rhs,r}=\|k_r^{(\times)}-k_r^{(+)}\|_\infty
\]

\[
\boxed{E_{state}=\|\Psi_{n+1}^{(\times)}-\Psi_{n+1}^{(+)}\|_\infty}
\]

\[
\boxed{E_{state,rel}=\frac{E_{state}}{\max(1,\|\Psi_{n+1}^{(+)}\|_\infty)}}
\]

主要readoutには

\[
\boxed{E_{read,A}=|R_A^{(\times)}-R_A^{(+)}|}
\]

を保存する。

さらに identity-state drift、time-state consistency、initial-state equality、numerical-config equality を監査する。

\[
E_{state,rel}<\epsilon_{equiv}
\]

かつ主要readoutについて

\[
E_{read,A}<\epsilon_{readout-equiv,A}
\]

を必須とする。不一致時に片方を補正してはならない。

---

# 20. 第1主図：ユーザー指定3段図

**Figure 01: Orbit radius and radiation–absorption residuals**

Figure 01 を含む物理主図は、すべて乗法正式系 \(\Psi^{(\times)}\) から得たreadoutだけを描く。加法系は物理主曲線へ混在させない。

共通横軸：

\[
\boxed{T}
\]

### Panel (a)

\[
\boxed{r_{\rm orb}(T)}
\]

linear scaleを基本とする。

必要に応じて補助版としてlog yも出す。

### Panel (b)

\[
\boxed{\Delta F_{\rm GW}(T)}
\]

ゼロ線を明示。

signed symlog を標準追加版として作る。

### Panel (c)

\[
\boxed{\Delta F_{\rm EM}(T)}
\]

ゼロ線を明示。

signed symlog を標準追加版として作る。

ファイル名：

```text
figure01_orbit_radius_gw_em_balance_linear.png
figure01_orbit_radius_gw_em_balance_linear.svg
figure01_orbit_radius_gw_em_balance_symlog.png
figure01_orbit_radius_gw_em_balance_symlog.svg
```

---

# 21. なぜ symlog 版を最初から作るか

後に非常に深い谷や極小残差が現れた場合、

通常のlinear plotでは

\[
10^{-4},10^{-8},10^{-12}
\]

などが全てゼロに見える。

一方、通常のlogでは符号とゼロ交差が失われる。

したがって、

\[
\boxed{\text{signed symlog}}
\]

を標準の補助図とする。

linear版も必ず残し、symlogだけで判断しない。

---

# 22. Figure 02：物理的散逸と軌道縮退

共通横軸 \(T\)。

### Panel (a)

\[
r_{\rm orb}(T)
\]

### Panel (b)

\[
D_{\rm GW}(T)
\]

### Panel (c)

\[
D_{\rm EM}(T)
\]

### Panel (d)

\[
D_{\rm rad}(T)
\]

ファイル名：

```text
figure02_orbit_and_positive_energy_drain.png
figure02_orbit_and_positive_energy_drain.svg
```

目的：

\[
r_{\rm orb}\text{ の縮小}
\]

と

\[
\text{正味energy drain}
\]

の相関を見る。

---

# 23. Figure 03：放射と吸収の分解

### Panel (a)

\[
F_{\rm GW}^{\rm out}(T)
\]

と

\[
F_{\rm GW}^{\rm abs}(T)
\]

### Panel (b)

\[
F_{\rm EM}^{\rm out}(T)
\]

と

\[
F_{\rm EM}^{\rm abs}(T)
\]

positive fluxについてはlog y版も作る。

ファイル名：

```text
figure03_gw_em_out_abs_decomposition_linear.png
figure03_gw_em_out_abs_decomposition_log.png
```

---

# 24. Figure 04：軌道縮退速度

### Panel (a)

\[
r_{\rm orb}(T)
\]

### Panel (b)

\[
\dot r_{\rm orb}(T)
\]

### Panel (c)

\[
\frac{d\log r_{\rm orb}}{dT}
\]

raw finite differenceを保存する。

表示用には orbit-average 版を重ねてもよいが、rawを消してはならない。

目的：

- 本当に plateau へ向かうか
- 単に非常に遅い inspiral なのか

を区別する。

---

# 25. Figure 05：位相・角速度

\[
\phi(T)
\]

\[
\Omega(T)
\]

\[
\dot\Omega(T)
\]

を図化する。

有限半径に見えても、

\[
\Omega(T)
\]

が長期ドリフトしている場合は真の準定常ではない可能性がある。

---

# 26. Figure 06：a,b の閉曲面上の量（v3.1）

最低限、

\[
M_{\Sigma,a},M_{\Sigma,b}
\]

\[
Q_{\Sigma,a},Q_{\Sigma,b}
\]

\[
|\mathbf J_{\Sigma,a}|,|\mathbf J_{\Sigma,b}|
\]

を共通時間軸で図化する。複数の閉曲面で読んだ値を重ねる。

放射吸収との対応を確認する。

---

# 27. Figure 07：系の時計・内部時計・ロック（v3.1）

v3.0 の Figure 07（\(h,h/2,h/4\) の \(\Delta\tau\) convergence）は撤回する（§9）。

共通横軸 \(T\)（および系の時計の周期数 \(N_{\rm period}\) を横軸にした版）。

### Panel (a)

\[
\delta\tau(T)
\]

### Panel (b)

\[
\Omega_{\rm sys}(T),\quad \omega_{{\rm int},a}(T),\quad \omega_{{\rm int},b}(T)
\]

### Panel (c)

\[
\rho_a(T)=\frac{\Omega_{\rm sys}}{\omega_{{\rm int},a}},\quad
\rho_b(T)=\frac{\Omega_{\rm sys}}{\omega_{{\rm int},b}}
\]

### Panel (d)

\[
\hbar_{{\rm eff},A}(T),\quad \alpha_{{\rm eff},A}(T)
\]

raw 点を全表示し、平滑化しない。

ファイル名：

```text
figure07_system_clock_internal_clock_lock.png
figure07_system_clock_internal_clock_lock.svg
```

---

# 28. Figure 08：N_cut convergence

\[
N_L,N_M,N_H
\]

について、

\[
r_{\rm orb}(T)
\]

\[
\Delta F_{\rm GW}(T)
\]

\[
\Delta F_{\rm EM}(T)
\]

を比較する。

高周波structureが \(N_{\rm cut}\) とともに移動・消失する場合はspectral artifactを疑う。

---

# 29. Figure 09：抽出半径 convergence

複数の

\[
R_{{\rm ext},k}
\]

について、

\[
F_{\rm GW}^{\rm out}(T;R_{{\rm ext},k})
\]

\[
F_{\rm EM}^{\rm out}(T;R_{{\rm ext},k})
\]

を比較する。

retarded-time alignmentを行う場合は、raw time と aligned time の双方を保存する。

---

# 30. Figure 10：run 間比較（v3.1）

Run H0（スピン反平行）、H1（平行）、H2（回転なし）、H3（質量比・\(C_0\) の変更）を比較する。

横軸候補：run 条件（スピンの向き、質量比、\(C_0\)）

縦軸候補：

\[
r_{\rm late}/(J/M_a),\quad
\langle\delta\tau\rangle_{\rm late}/P_{\rm sys},\quad
\langle\rho_A\rangle_{\rm late},\quad
\langle\alpha_{{\rm eff},A}\rangle_{\rm late},\quad
\langle|\Delta F_{\rm GW}|\rangle_{\rm late},\quad
\langle|\Delta F_{\rm EM}|\rangle_{\rm late}
\]

をscatterとして出す。

**線で滑らかに補間しない。** 全run点を表示する。\(m/n\) 軸への変換は禁止しない（v3.0 の禁止は撤回）。

---

# 30A. Figure 11：乗法–加法同値性監査

物理主図とは分離して、検算専用図を作る。

- Panel (a): \(E_{state,rel}(T)\)
- Panel (b): \(E_{stage,2},E_{stage,3},E_{stage,4}\)
- Panel (c): \(\Delta r_{\rm orb},\Delta(\Delta F_{\rm GW}),\Delta(\Delta F_{\rm EM})\)
- Panel (d): \(|C_{\delta\tau}|\) と identity-state drift

誤差が0を含むため必要に応じてlinearまたはsymlogを使う。これは物理結果図ではなく、万能相互作用の表現同値性を監査する図である。

---

# 31. 「後で深い谷を見落とさない」図化ルール

第1実験から次を徹底する。

1. raw scatter pointを全表示
2. smoothing curveは主図に使わない
3. interpolationは可視化補助に限り、raw点を必ず重ねる
4. signed residualを保存
5. linear版とsymlog版を両方作る
6. zero crossingの前後を保持
7. exact initial parameter representationを保存
8. local minimumを自動削除しない
9. parameter binningしない
10. runごとの最終値だけでなくtime seriesを残す
11. cadence down-sampling前のraw seriesを保管
12. later analysis用に run ID と全state initializationを一意に結びつける

---

# 32. 準定常ラベルの付け方

準定常かどうかは、結果を記述するためのラベルであり、実験の成功条件ではない。予想（ロックすれば準定常になる、など）が外れても系の動作には影響しない。

固定した任意値 \(|\dot r|<10^{-k}\) だけでラベルを付けない。`δτ` の系列、\(\Omega_{\rm sys}\) の変化、空間 cutoff の比較、抽出半径の比較、constraint residual を並べて記録し、late-time window で

\[
\left|\frac{d\log r}{dT}\right|,\qquad
\left|\frac{d\log\delta\tau}{dT}\right|
\]

が小さい状態が複数周期続く場合を plateau candidate として記録する。

これは**readout 上のラベル**であり、発展へ帰還させない。

---

# 33. outcome分類

各runについて、最終的に次のラベルを出す。

## `INSPIRAL`

\[
\dot r<0
\]

が長時間継続し、plateau条件を満たさない。

## `MERGER`

閉曲面 \(\Sigma_a,\Sigma_b\) が互いを分離して囲めなくなる（v3.0 の common horizon 判定は撤回）。

## `QUASI_STATIONARY`

\[
r\to r_*>0
\]

かつlate-time slopeがnumerical floor内。

## `QUASI_STATIONARY_RADIATIVE`

半径はplateauに入るが、

\[
D_{\rm rad}
\]

が数値floorより有意に残る。

## `EXPANDING`

長時間で

\[
\dot r>0
\]

## `OSCILLATORY`

有限振幅の長期振動が残り、単純plateauでない。

## `NUMERICALLY_UNRESOLVED`

constraint residual が大きい、または空間 cutoff の比較で大きく結果が変わるなど、記録から数値上の問題が疑われる。

以上のラベルは、結果を記述するために後から付けるものであり、予想でも成功条件でもない。

---

# 34. 「放射が止まる」と「放射–吸収差がゼロ」を区別する

次の3状態を明確に区別する。

### Case 1

\[
F^{\rm out}\to0,\quad F^{\rm abs}\to0
\]

真にradiation channel自体が静まる。

### Case 2

\[
F^{\rm out}\simeq F^{\rm abs}\neq0
\]

差分はゼロだが、flux自体は残る。

### Case 3

\[
F^{\rm out}\neq F^{\rm abs}
\]

明確なbalance residual。

したがって、

\[
\Delta F\to0
\]

だけで「放射停止」と書いてはならない。

---

# 35. 第1実験の最小run順序

計算コストを抑えつつ物理と数値を切り分けるため、以下の順で走らせる。

## Stage V0: Minkowski / zero-field

万能演算が自明状態を壊さないことを確認する（実装検算）。

## Stage V1: 回転する一体

電子型（\(Q=-3\)、\(J=\hbar/2\)）と陽子型（\(Q=+3\)、\(J=\hbar/2\)）の一体を、それぞれ単独で発展させる。閉曲面上の \(M,Q,J\)、内部時計 \(\omega_{\rm int}=M/J\)、磁気回転比 \(g\) を読み出して記録する。

## Stage H0: 水素型・スピン反平行（基準）

## Stage H1: 水素型・スピン平行

## Stage H2: 水素型・回転なし（対照）

## Stage H3: 質量比・\(C_0\) の変更

各 run は最低 2 周期（62 刻み、段で 248）以上を記録する（§5.1）。

v3.0 の Stage V2（Majumdar–Papapetrou の静的均衡 anchor）、B1〜B4（連続的な \(q\) の family）、C（\(\Delta\tau\) 3 水準の収束）は撤回する。MP の均衡は同符号・\(Q=M\) の極限ブラックホールを前提とし、\(q_0\) の量子化とも合わない。

---

# 36. 第1実験での「均衡」の意味

本実験で重力–電磁の均衡を見るとき、単一のNewtonian force equalityだけを使わない。

少なくとも、

1. `δτ`（系の時計）の変化
2. 軌道半径と相対回転の角振動数の変化
3. GW outward flux
4. GW の a,b への吸収（閉曲面）
5. EM outward flux
6. EM の a,b への吸収（閉曲面）
7. 閉曲面上の \(M,Q,J\) の変化と内部時計
8. 系の時計と内部時計の比
9. constraints

を同時に記録する。どの状態を「均衡」と呼ぶかは、結果が出てから解釈する。

---

# 37. raw time-series file schema

正式な物理time-seriesは乗法系を正本とし、1 row / readout step とする。

```text
run_id
step
T_tau_mul
T_tau_add
T_read_mul
T_read_add
delta_tau
D_tau
kappa_tau
T_read_difference
time_state_consistency_error
N_cut

xa_x
xa_y
xa_z
xb_x
xb_y
xb_z
separation_d
orbit_radius_r
orbit_phase_phi
orbit_omega
dr_dT

gw_flux_out
gw_flux_abs_sigma_a
gw_flux_abs_sigma_b
gw_flux_abs_total
gw_out_minus_abs
gw_out_plus_abs

em_flux_out
em_flux_abs_sigma_a
em_flux_abs_sigma_b
em_flux_abs_total
em_out_minus_abs
em_out_plus_abs

total_positive_drain

N_tick
Omega_sys
P_sys
phi_sys
N_period
tau_proper_a
tau_proper_b
redshift_a
redshift_b

Ma_sigma
Mb_sigma
Qa_sigma
Qb_sigma
Jax_sigma
Jay_sigma
Jaz_sigma
Jbx_sigma
Jby_sigma
Jbz_sigma
sigma_surface_id
omega_int_a
omega_int_b
lock_ratio_a
lock_ratio_b
lock_ratio_a_convergent_m
lock_ratio_a_convergent_n
lock_ratio_b_convergent_m
lock_ratio_b_convergent_n
hbar_eff_a
hbar_eff_b
alpha_eff_a
alpha_eff_b
magnetic_moment_a
magnetic_moment_b
g_factor_a
g_factor_b

gh_constraint_L2
gh_constraint_Linf
maxwell_divE_L2
maxwell_divE_Linf
maxwell_divB_L2
maxwell_divB_Linf

stage2_state_error_inf
stage3_state_error_inf
stage4_state_error_inf
k1_error_inf
k2_error_inf
k3_error_inf
k4_error_inf
final_state_error_inf
final_state_error_rel
identity_state_error_inf
initial_state_error_inf
numerical_config_error_inf
```

上記の物理readout列は乗法正式系の値である。

主要readoutの加法値と差分は別ファイル `readout_equivalence.*` に保存する。最低限、orbit radius、GW/EM balance residual、閉曲面上の M/Q/J、`δτ` の mul/add/diff を残す。

---

# 38. waveform/raw-mode file schema

GW は正式な乗法系について `gw_modes_mul.*` に

```text
T
R_ext
ell
m
Re_Psi4
Im_Psi4
```

を保存する。加法検算系も同じschemaで `gw_modes_add.*` に保存し、差を `gw_modes_diff.*` に保存する。

EM も同様に

```text
T
R_ext
ell
m
Re_EM_mode
Im_EM_mode
```

を `em_modes_mul.*`, `em_modes_add.*`, `em_modes_diff.*` に保存する。

可能ならbinary format/HDF5を用い、CSVはsummary用とする。

---

# 39. run manifest

各runにつき `run_manifest.json` を作り、最低限、

```text
interaction_spec_version = v3.2
initialization_spec_version = v3.1
readout_spec_version = v3.1
unit_convention = G=q0=c=1
charge_convention = GAUSSIAN_GEOMETRIZED
alpha_input
J_input
N_tick
closed_surface_definitions
code_version
state_schema_version
rhs_evaluator_version
rhs_evaluator_hash
time_integrator_id = RK4_CLASSICAL
rk4_tableau
multiplicative_operator_version
all_initial_state_values
initial_state_hash
exact_parameter_forms
basis_definition
domain_definition
boundary_definition
extraction_radii
delta_tau
D_tau
kappa_tau
N_cut
solver_tolerances
epsilon_equiv
epsilon_readout_equiv
epsilon_time_state
stopping_conditions
readout_cadence
checkpoint_cadence
```

を保存する。`dt` や外部current timeを正本値としてmanifestへ置かない。時間刻みの正本は状態の `delta_tau,D_tau,kappa_tau` であり、`delta_tau,D_tau` は系の時計から毎 step 定まる（初期値のみ manifest に記録する）。

---

# 40. 推奨出力ファイル構成

```text
run_<RUNID>/
    run_manifest.json
    initial_state_master.*
    initial_state_mul.*
    initial_state_add.*
    final_state_mul.*
    final_state_add.*
    timeseries_raw_mul.*
    readout_equivalence.*
    state_equivalence_history.*
    time_state_audit.*
    closed_surface_readout_mul.*
    closed_surface_readout_add.*
    closed_surface_readout_diff.*
    clock_lock_readout_mul.*
    gw_modes_mul.*
    gw_modes_add.*
    gw_modes_diff.*
    em_modes_mul.*
    em_modes_add.*
    em_modes_diff.*
    energy_budget_mul.*
    constraint_history_mul.*
    constraint_history_add.*
    figures/
        figure01_orbit_radius_gw_em_balance_linear.png
        figure01_orbit_radius_gw_em_balance_symlog.png
        figure02_orbit_and_positive_energy_drain.png
        figure03_gw_em_out_abs_decomposition_linear.png
        figure04_orbit_decay_rate.png
        figure05_phase_omega.png
        figure06_closed_surface_quantities.png
        figure07_system_clock_internal_clock_lock.png
        figure08_spatial_convergence.png
        figure09_extraction_radius_convergence.png
        figure10_run_comparison.png
        figure11_multiplicative_additive_equivalence.png
```

物理主図は乗法系だけを用いる。`figure11` だけが二経路比較を主目的とする。

---

# 41. 図化の共通規則

1. 物理主図は乗法正式系のreadoutだけを使用する
2. 横軸は原則 `T_read_mul`
3. raw data を基準にする
4. smoothingした場合は必ず別series
5. raw seriesを消さない
6. flux differenceはzero lineを描く
7. signed quantityに通常logを使わない
8. signed residualにはsymlogを使う
9. positive fluxにはlog版を追加可
10. 軸単位とnormalizationを図中に明記
11. run IDをcaptionかfilenameに含める
12. transient除外を行った場合、除外区間を明示
13. late-time averageだけを見せずfull historyも残す
14. interpolation curveだけの図は禁止
15. parameter sweepは全raw pointを表示する
16. y-rangeを手動で狭めて変動を誇張しない
17. 極小構造を潰さないためsymlog補助図を作る
18. 加法系は物理主図に重ねず、equivalence監査図またはdifference図に限定する

---

# 42. 第1実験の成功条件（v3.1）

v3.0 の完了条件（\(\delta\tau\to0\) での収束、乗法–加法同値性、各 family の分類）は撤回する。予想や分類は成功条件にしない（§0A の 7）。

第1実験の成功条件は、系の動作そのものに関する次の項目に限る。

1. **方程式を変えずに適用していること**：Einstein–Maxwell 方程式、Kerr–Newman の厳密式を変更していない。正則性ゲート・excision・地平面境界条件を使っていない。
2. **隠れた仮定がないこと**：単位 \(G=q_0=c=1\)、Gauss 規約、`δτ` の規則（\(N_{\rm tick}=31\)）、初期スピン \(J=\hbar/2\) など、与えた値はすべて状態または manifest に exact 値で記録されている。
3. **規則どおりに状態を更新していること**：正式系は乗法更新 \(\Psi_{n+1}=S(\Psi_n)\Psi_n\) だけで生成されている。`δτ` は状態から定まっている。readout を発展へ帰還させていない。
4. **実装検算**：乗法系と加法系が丸め誤差の範囲で一致し、time-state consistency と identity-state の監査を満たす。これは実装が意図どおりかの確認である。
5. **再現性**：同じ条件で同じ結果が再現される（code、初期状態、raw data、ハッシュを保存）。
6. **全読み出しの保存**：§11〜§19 の読み出しと raw time series がすべて保存されている。

ロックするかどうか、ロック比、\(\alpha_{\rm eff}\)・\(\hbar_{\rm eff}\) の値、閉曲面の選び方への依存、準定常かどうかは、結果として記録し、後から解釈する。

---

# 43. 第2段階へ進む条件（v3.1 で撤回）

v3.0 は「continuum run の収束確認、\(\delta\tau\) を数値刻みとして小さくした極限の確認、Majumdar–Papapetrou anchor の再現などを満たすまで \(m/n\) 探索へ進まない」としていた。これは撤回する。\(m/n\) の読み出しは第1実験から行う（§3、§18.3）。

---

# 44. 本実験の中心図

本仕様で最重要の図は次の3段図である。

\[
\boxed{
\begin{array}{c}
r_{\rm orb}(T)\\
\Delta F_{\rm GW}(T)\\
\Delta F_{\rm EM}(T)
\end{array}
}
\]

と、系の時計・内部時計・ロックの図（Figure 07）である。

これにより、

- 軌道が縮む
- 放射–吸収balance residualが減る
- 有限半径へ近づく
- residualが残る
- mergerへ向かう

を同一時間軸上で直接比較する。

ただし物理的energy drainの判定には必ず

\[
D_{\rm GW}=F_{\rm GW}^{\rm out}+F_{\rm GW}^{\rm abs}
\]

\[
D_{\rm EM}=F_{\rm EM}^{\rm out}+F_{\rm EM}^{\rm abs}
\]

も併記または別図で監査する。

---

# 45. 第1実験で問うこと

\[
\boxed{
\text{背景時空も粒子も先に置かず、厳密な Einstein–Maxwell 方程式を乗法更新で発展させたとき、水素型二体 }a,b\text{ の状態の変化から、何が読み出されるか}
}
\]

特に、系の時計 `δτ`、a,b の内部時計、その比、\(\hbar_{\rm eff}\)、\(\alpha_{\rm eff}\)、放射と吸収の収支を、完全な raw time series として残す。与えた値（\(\alpha\)、\(J=\hbar/2\)、\(N_{\rm tick}=31\)）と読み出した値が自己無撞着に閉じるかは、結果として記録し、後から解釈する。

---

# 46. 参考になる標準物理上の対応

- Generalized harmonic numerical relativity では \(\Psi_4\) を用いたGW readoutが標準的に使用できる。
- 電荷は、囲む閉曲面を貫く電束（Gauss の法則）およびループの位相（ホロノミー、Aharonov–Bohm 位相）で、内側を見ずに定義できる。質量・角運動量には Komar 型積分がある。
- Kerr–Newman 解に電子のパラメータを入れると、磁気回転比は \(g=2\) で Dirac 電子と一致する（Carter 1968）。回転半径 \(J/M\) は換算 Compton 波長の半分になる（Burinskii の Kerr–Newman 電子模型）。
- 電子の内部時計（Compton 振動数）と軌道の位相が揃う条件から Bohr の量子化条件が出る（de Broglie 1924）。
- 束縛した二体系の時計はらせん Killing ベクトルに対応して一つであり、各天体の固有時計は赤方偏移 \(z_A\) で結ばれる（Friedman–Uryū–Shibata 2002、Detweiler 2008、Le Tiec–Blanchet–Whiting 2012）。
- 反復交換散乱の有限位数共鳴（木原、Concept DOI 10.5281/zenodo.21421366）：\(R_{n,m}=\cos^2(\pi m/n)\)、\(\alpha^{-1}\approx137\) に対応する根 \(m/n=23/124\)。
- 本研究ではこれらを外部近似力として dynamics に注入せず、読み出しの解釈にのみ使用する。

参考：
1. SpECTRE General Relativity documentation, `gr::psi_4`.
2. B. Carter, Phys. Rev. 174, 1559 (1968).
3. L. de Broglie, Recherches sur la théorie des quanta (1924).
4. S. Detweiler, Phys. Rev. D 77, 124026 (2008).
5. A. Le Tiec, L. Blanchet, B. F. Whiting, Phys. Rev. D 85, 064039 (2012).
6. J. L. Friedman, K. Uryū, M. Shibata, Phys. Rev. D 65, 064035 (2002).
7. 木原範昭, 反復交換散乱における有限位数共鳴の発見, Zenodo, Concept DOI 10.5281/zenodo.21421366 (2026).

---

# 47. 仕様上の最終固定事項

第1実験では、

\[
\boxed{
\text{方程式は変えない。背景時空も粒子も先に置かない。与えた値はすべて記録する。}
}
\]

を守る。

保存・図化は、

\[
\boxed{
\text{raw data 完全保存、狭い谷を見落とさない記録}
}
\]

とする。

v3.0 の「continuous classical first」「量子化仮説を初期の parameter selection へ持ち込まない」は撤回する。本系の初期値（\(J=\hbar/2\)、\(N_{\rm tick}=31\)、1 刻み内の 4 段）は \(\alpha\) から予測値として与えたものであり、隠れたパラメータではない。任意に与えたパラメータが自己無撞着になるには特定の値でなければならないことを示そうとするものである（§0A の注意書き）。
