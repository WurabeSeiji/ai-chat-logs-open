# Einstein–Maxwell に厳密に従う二体 a,b の動力学の導出 v3
## ― 半角の担い手を状態とし、4 段を weight 1/2→1→3/2→2 の冪として、E–M の各項を weight の段に置く。δτ は節点時刻の変分で決まる ―

**著者:** 木原範昭 (Noriaki Kihara)  
**整理:** Claude  
**版:** v3  
**日付:** 2026-10-01  
**置換対象:** `einstein_maxwell_strict_two_body_expansion_ja_v2_20261001.md`（v1・v2 は残置）  
**位置づけ:** 木原氏の指示「Einstein–Maxwell に厳密に従う」に基づく導出。相互作用仕様 v3.2・初期化仕様 v3.1・読出し仕様 v3.1、固定点 v1〜v4 のうち本書と食い違う箇所は本書を優先する。

---

# 0A. v2 からの訂正と決定

1. **段の構造（木原氏合意）**：状態に半角の担い手 \(V=e^{i\phi/2}\)（a, b の内部位相、weight 1/2）を持つ。1 刻みの 4 段は \(V,V^2,V^3,V^4\)（weight 1/2, 1, 3/2, 2）であり、E–M の各項は媒介する場のスピンに等しい weight の段に置く（Maxwell は weight 1、Einstein は weight 2）。Thomas 歳差は weight 1/2 の段から出る。段の位相 \(\sigma_s=i^s\) は、担い手がその段で獲得する四分の一回転である。
2. **読み A（木原氏決定）**：有限位数共鳴の一回の反復 = 一段。124 段 = 31 刻みが状態の再帰（半角の \(4\pi\)）。軌道一周は 62 段 = **15.5 刻み**。一周の終わりは刻みの第 2 段に落ち、そこで \(\sigma=-1\)、\(V\to-V\)。
3. **δτ の変分式**：節点時刻の変分から \(E_d(k-1)=E_d(k)\)（エネルギーの連続）。v2 の「\(\partial L_d/\partial\delta\tau_k=0\)」は誤り。
4. **ALD**：相対座標の加速度と力の因子 \(\mu\) を区別。
5. **二体のスピン–軌道**：Barker–O'Connell の二体式。Thomas 項は段の構造から、磁気の項は weight 1 の段に。
6. **Kerr–Newman の四重極**：\(-m_Aa_A^2\)（weight 2）、\(-Q_Aa_A^2\)（weight 1）を追加。\(\alpha^2\) の次数。
7. **B1・B2**：反射する境界の内側の放射場を定在モードの振幅（波の状態）として持つ。有限の履歴の窓は使わない。
8. **速度**：増分で定義。履歴は 3 階以上の微分だけ。
9. **陰性**：速度依存の項を含む離散 E–L 方程式は \(z_{k+1}\) について陰的。初期化で小さい項について解を展開し、陽的な積和にする。
10. **初期の円軌道**：離散 E–L 方程式を \(\Delta\phi\) 一定・弦一定で解いて決める。

---

# 0. E–M に厳密に従うと、a, b の運動はどこから出るか

Einstein–Maxwell 電磁真空

\[
G_{ab}=8\pi T^{\rm EM}_{ab},\qquad \nabla_aF^{ab}=0,\qquad dF=0
\]

には物質場がない。電荷・質量・スピンを持つ a, b は、場の方程式の解の中で**位相が閉じない場所**（Kerr–Newman 型の場の塊）としてだけ存在する。粒子を別に置かないという原則 2 は、E–M そのものの性質である。

a, b の運動も E–M が決める。Einstein–Infeld–Hoffmann（1938）は重力場の方程式だけから特異点の運動方程式を導き、Infeld–Wallace（1940）は同じ方法を Maxwell 場に適用して Lorentz 力と放射反作用を場の方程式から導いた。導出の形は、**場を初期配置の周りで展開し、各次数の係数を求める**ことである。これが「初期化で展開する」ことの E–M 側の実体であり、得られる項が離散作用 \(L_d\) の中身になる。

連続の軌道との一致は、事前に確かめる基準にしない。離散系がそれ自身で閉じた力学であり、結果を読み出して解釈する。

---

# 1. 展開の小さい量と、L_d に入れる項の次数

展開の小さい量は、初期条件（\(Q=\mp3\)、\(J=\hbar/2\)、\(L_0=2J\)）が決める軌道速度 \(v_0=K/L_0\approx\alpha\) と、Kerr–Newman の回転半径と分離の比 \(a_A/R\)（\(a_A=J_A/m_A\)。電子型で \(\approx\alpha/2\)）である。主力（クーロン）に対する相対的な大きさ（重力／クーロン \(=1/C_0\approx10^{-4}\)）：

| 項 | weight（段） | 相対的な大きさ | 値 |
|---|---|---|---|
| クーロン | 1 | 1 | — |
| Newton 重力 | 2 | \(1/C_0\) | \(10^{-4}\) |
| Darwin、磁気スピン–軌道 | 1 | \(\alpha^2\) | \(5.3\times10^{-5}\) |
| Thomas 歳差 | 1/2 | \(\alpha^2\) | \(5.3\times10^{-5}\) |
| KN 電荷四重極 \(-Q_aa_a^2\) | 1 | \((a_a/R)^2\) | \(1.3\times10^{-5}\) |
| スピン–スピン | 1 | \(\alpha^2m_a/(4m_b)\) | \(7\times10^{-9}\) |
| 電磁の放射反作用（ALD） | 1 | \(\tfrac23\alpha^3\) | \(2.6\times10^{-7}\) |
| 電磁の \(1/c^4\) | 1 | \(\alpha^4\) | \(2.8\times10^{-9}\) |
| EIH 1PN、重力–電磁 1PN 交差項、Lense–Thirring、KN 質量四重極 | 2 | \(\alpha^2/C_0\) 程度 | \(5\times10^{-9}\) |
| 重力波の放射反作用（Burke–Thorne） | 2 | \(\alpha^5/C_0\) 程度 | \(10^{-11}\) |
| 電磁の \(1/c^6\)、KN の \((a/R)^4\) | 1 | \(\alpha^6\)、\((a/R)^4\) | \(10^{-13}\)、\(2\times10^{-10}\) |

相対的な大きさ \(\varepsilon\) の項を L_d から落としたときの位相のずれは N 周（軌道）で約 \(2\pi\varepsilon N\)。読み出しの解像度 1 段 \(=2\pi/62\)（一周 62 段）に対して

\[
\boxed{\varepsilon\,N<\frac1{62}}
\]

を満たす次数まで、**その次数の全項**を入れる。放射で軌道が縮む run（B3、\(N\approx3\times10^5\) 周）では \(\alpha^4\) まで。ロックを読む run（状態の再帰 2 回 = 248 段 = 4 周）では \(\alpha^2\) まで。重力波の反作用は解像度の下だが、E–M の項として入れる。

---

# 2. 状態

| 区分 | 状態 | 表現・備考 |
|---|---|---|
| 半角の担い手 | \(V=e^{i\phi/2}\) | weight 1/2。\(V^2=e^{i\phi}=n\)（相対ベクトルの向き）、\(V^4=e^{2i\phi}\)（aa・ab・bb の位相）は導出 |
| 分離 | \(e^{\rho}\)、\(e^{-\rho}\) | \(r=e^{\rho}\)、\(1/r^n=e^{-n\rho}\) |
| 直前の節点 | 前段・前々段・三段前の \((V,e^{\rho})\) | 増分（速度）と 3 階差分（ALD、Burke–Thorne の \(d^5I/dt^5\) は五つ）に使う |
| スピンの向き | \(e^{i\psi_a}\)、\(e^{i\psi_b}\)、軌道軸成分の符号 | 大きさ \(\hbar/2\) は定数 |
| 空洞のモード（B1・B2） | \(A_{\ell=1,n}\)（電磁）、\(A_{\ell=2,n}\)（重力）、\(n\) は波長 \(\ge c\,\delta\tau\) の動径モード | 複素振幅（波の状態）。B3 では持たない |
| 時間 | \(\delta\tau_k\)（節点時刻の変分で決まる）、\(1/\delta\tau_k\)（補助）、\(\sigma_s=i^s\)、\(T_\tau\)、\(D_\tau\) | \(T_\tau'=D_\tau T_\tau\)、\(D_\tau\) は \(\delta\tau\) の相対変化の級数で更新 |
| 定数 | ONE、\(Q_a\)、\(Q_b\)、\(m_a\)、\(m_b\)、\(\alpha\)、\(J\)、\(a_A=J_A/m_A\)、\(f\)、境界距離 \(D\)、\(\kappa\)、展開の次数 | 恒等 |

a, b の個別の位置は、相対位置と 1PN の重心関係から読み出す。

---

# 3. L_d に入る項と、置く段（\(G=c=1\)、Gauss。\(n=V^2\)、\(r=e^{\rho}\)、\(\mathbf L\) は軌道角運動量）

## 3.1 weight 1/2 の段（\(V\)）

- 担い手の位相の進み \(V\to Ve^{i\Delta\phi/2}\)。
- スピンの Thomas 歳差：各体のスピンに対し、二体スピン–軌道（Barker–O'Connell）

\[
H_{SO}=-\sum_A\frac{Q_AQ_B}{2r^3}\left[\frac{g_A-1}{m_A^2}+\frac{g_A}{m_Am_B}\right]\mathbf L\cdot\mathbf S_A
\]

のうち、運動学的な Thomas 項 \(+\dfrac{Q_AQ_B}{2m_A^2r^3}\mathbf L\cdot\mathbf S_A\)（括弧の \(-1/m_A^2\) に当たる）。これは半角の位相が一周で \(\pi\) しか進まないことの表れであり、別の項として足さず、この段の担い手の回転から生成する。

## 3.2 weight 1 の段（\(V^2\)、Maxwell）

\[
L_C=-\frac{Q_aQ_b}{r}
\]

\[
L_D=+\frac{Q_aQ_b}{2r}\left[v_a\cdot v_b+(v_a\cdot n)(v_b\cdot n)\right]
\quad(\text{Landau–Lifshitz §65})
\]

磁気スピン–軌道（上の二体式の \(g_A/m_A^2+g_A/(m_Am_B)\) の部分）：

\[
H_{SO}^{\rm mag}=-\sum_A\frac{Q_AQ_B\,g_A}{2r^3}\left[\frac1{m_A^2}+\frac1{m_Am_B}\right]\mathbf L\cdot\mathbf S_A,\qquad g_A=2
\]

スピン–スピン：

\[
H_{SS}=\frac1{r^3}\left[\boldsymbol\mu_a\cdot\boldsymbol\mu_b-3(\boldsymbol\mu_a\cdot n)(\boldsymbol\mu_b\cdot n)\right],\qquad
\boldsymbol\mu_A=\frac{Q_A}{m_A}\mathbf S_A
\]

Kerr–Newman の電荷四重極 \(Q^{(2)}_A=-Q_Aa_A^2\)（軸はスピンの向き）と相手の電荷の結合：

\[
H_{Q2}=\sum_A\frac{Q_B\,Q^{(2)}_A}{2r^3}\left[3(\hat s_A\cdot n)^2-1\right]
\]

電磁の \(1/c^4\) 項（Golubenkov–Smorodinskii）：E–M から導かれた閉じた式として文献から取り、初期化で係数にする。

放射反作用（双極子次数）：相対座標の加速度

\[
\ddot{\mathbf z}_{RR}=\tfrac23\,\mu\left(\lambda_a-\lambda_b\right)^2\dot{\mathbf a},\qquad
\lambda_A=\frac{Q_A}{m_A}
\]

（力は \(\mu\) 倍。各体の自己力と相手の放射場の反作用の和。\(\dot{\mathbf a}\) は直前 4 節点の 3 階差分。）B3 ではこの局所項。B1・B2 では §3.4 のモードの反作用がこれを含む。

## 3.3 weight 3/2 の段（\(V^3\)）

E–M にスピン 3/2 の場はない。この段は自由運動（運動項だけ）。

## 3.4 weight 2 の段（\(V^4\)、Einstein）

\[
L_N=\frac{m_am_b}{r}
\]

\[
L_{1PN}=\tfrac18\left(m_av_a^4+m_bv_b^4\right)
+\frac{m_am_b}{2r}\left[3(v_a^2+v_b^2)-7v_a\cdot v_b-(v_a\cdot n)(v_b\cdot n)\right]
-\frac{m_am_b(m_a+m_b)}{2r^2}
\]

重力–電磁の 1PN 交差項（\(Q_aQ_b(m_a+m_b)/r^2\) 型）：Khalil–Buonanno–Steinhoff–Vines 2018 の Einstein–Maxwell 極限から、初期化で係数にする。

Lense–Thirring（重力のスピン–軌道）：Barker–O'Connell 1975 の重力二体式

\[
H_{LT}=\sum_A\frac1{r^3}\left[2+\frac32\frac{m_B}{m_A}\right]\mathbf L\cdot\mathbf S_A
\]

（Barker–O'Connell 1975、Damour–Schäfer 1988 の 1PN スピン–軌道。\(G=c=1\)。）この係数のうち運動学的な Thomas 部分は §3.1 の段で生成されるので、重複して足さない。

Kerr–Newman の質量四重極 \(M^{(2)}_A=-m_Aa_A^2\) と相手の質量の結合：

\[
H_{M2}=-\sum_A\frac{m_B\,M^{(2)}_A}{2r^3}\left[3(\hat s_A\cdot n)^2-1\right]
\]

放射反作用（Burke–Thorne）：

\[
\Phi_{RR}=\frac15x^ix^j\frac{d^5I_{ij}}{dt^5},\qquad I_{ij}\ \text{は跡なし質量四重極}
\]

（\(d^5/dt^5\) は直前 6 節点の 5 階差分。）B3 ではこの局所項。

## 3.5 境界と空洞のモード（B1・B2）

境界（距離 \(D\)、反射率 \(1-f\)）の内側の放射場は定在モードになる。モード振幅 \(A_{\ell n}\) は波の状態で、各段で位相因子 \(e^{i\omega_{\ell n}\delta\tau/4}\) を掛け、往復時間 \(2D\) あたり \((1-f)\) に減衰する。駆動は \(\ell=1\) が \(d^{(3)}\)（weight 1 の段）、\(\ell=2\) が \(I^{(5)}_{ij}\)（weight 2 の段）。モードから a, b への反作用は、モード振幅の固定の線形結合で、同じ段に入る。

- B1（\(f=0\)）：減衰なし。放射エネルギーは全部モードに蓄積し、a, b に戻る。
- B2（\(0<f<1\)）：往復ごとに \((1-f)\)。
- B3（\(f=1\)）：モードを持たず、§3.2・§3.4 の局所の放射反作用項を使う。モード展開の \(f\to1\) の極限がこの局所項に一致することは初期化で確認する。

モードの数は、波長 \(\ge c\,\delta\tau\) の動径モード数（\(D/(c\,\delta\tau)\) 程度）。

---

# 4. 離散作用と変分

## 4.1 節点・刻み・段

節点 \(j\) は段の境界。刻み \(k\) は 4 段（\(j=4k,\dots,4k+3\)）からなり、各段の長さは \(\delta\tau_k/4\)。段 \(s=j\bmod4\) の位相は \(\sigma_s=i^s\)。段 \(s\) では weight \((s+1)/2\) の項（§3.1〜3.4）だけが働く。

離散作用

\[
S_d=\sum_j L_d^{(s)}(q_j,q_{j+1},\delta\tau_k/4)
\]

\(q=(V,e^{\rho},e^{i\psi_a},e^{i\psi_b},A_{\ell n})\)。運動項は 1 段の変位（弦）

\[
c_j^2=|z_{j+1}-z_j|^2=e^{2\rho_j}+e^{2\rho_{j+1}}-2e^{\rho_j+\rho_{j+1}}\cos\Delta\phi_j,\qquad
\cos\Delta\phi_j=\mathrm{Re}\!\left(V_{j+1}^2\,\overline{V_j^2}\right)
\]

で \(\mu c_j^2/(2\cdot\delta\tau_k/4)\)。相互作用項は段の両端の状態の対称な平均で評価する（時間反転対称を保つ選択。E–M ではなく離散作用の定義の選択であり、記録する）。速度は増分 \((z_{j+1}-z_j)/(\delta\tau_k/4)\)。

## 4.2 増分についての変分（位相刻みに対する傾斜）

\[
D_2L_d^{(s-1)}(q_{j-1},q_j)+D_1L_d^{(s)}(q_j,q_{j+1})=F_d^-(j-1)+F_d^+(j)
\]

右辺は放射反作用（B3）またはモードの反作用（B1・B2）の離散 Lagrange–d'Alembert の力。\(D_1L_d,D_2L_d\) が位相刻みに対する傾斜（離散運動量）。速度依存の項を含む段では \(q_{j+1}\) について陰的なので、初期化で小さい項（\(\alpha^2\) 以下）について解を §1 の次数まで展開し、\(q_{j+1}\) を OLD 状態の固定の積和で与える。\(1/\delta\tau_k\) は補助状態。

## 4.3 節点時刻の変分（δτ）

刻み \(k\) の 4 段は \(\delta\tau_k/4\) を共有する。節点時刻 \(t_k\)（刻みの境界）を変分すると、\(t_k\) は \(\delta\tau_{k-1}\) と \(\delta\tau_k\) の両方に入るので

\[
\boxed{E_d(k-1)=E_d(k)},\qquad
E_d(k)\equiv-\sum_{s=0}^{3}\frac{\partial L_d^{(s)}}{\partial(\delta\tau_k/4)}
\]

（Kane–Marsden–Ortiz 1999）。これが刻みの長さ \(\delta\tau_k\) を決める離散エネルギーの式で、E の値は初期データで決まる。放射があるときは離散 Lagrange–d'Alembert により \(E_d(k)-E_d(k-1)\) が反作用の離散の仕事に等しく、\(\delta\tau_k\) はそれに従って変わる。実装では比 \(\delta\tau_k/\delta\tau_{k-1}\) を、ほぼ 1 に近い量の固定項数の級数で更新する。\(\delta\tau_k\) の実解が存在しない刻みは FAIL として記録する（補正しない）。

**\(\delta\tau\) は設定しない。** 初期値 \(\delta\tau_0\) だけを、初期の離散円軌道で一周 = 62 段（15.5 刻み）、状態の再帰 = 124 段（31 刻み）となる値に置く。以後は上の式で決まる。一周が何段になるか、再帰が何段になるかは結果であり、読み出す。

## 4.4 時間と乗法更新

\[
T_\tau'=D_\tau T_\tau,\qquad
D_\tau'=D_\tau\left(1+\kappa\,\Delta\delta\tau/4+\tfrac12(\kappa\,\Delta\delta\tau/4)^2+\cdots\right)
\]

\[
V'=V\left(1+i\tfrac{\Delta\phi}{2}-\tfrac12\!\left(\tfrac{\Delta\phi}{2}\right)^2+\cdots\right),\qquad
e^{\rho'}=e^{\rho}\left(1+\Delta\rho+\tfrac12\Delta\rho^2+\cdots\right)
\]

級数は固定の項数で倍精度の丸め以下。遷移の中に \(\sqrt{\ }\)、\(\exp\)、\(\log\)、状態による割り算は現れない。

---

# 5. 読み出し（発展へ帰還させない）

- 担い手の位相 \(\phi_j/2\)、一周の段数、再帰（\(V\to V\)）の段数、一周の終わりが落ちる段と \(\sigma\) の値
- 1 段あたりの位相の進み \(\Delta\phi_j/(4\pi)\)（状態の再帰を 1 とする比）とその連分数 \(m/n\)
- \(\delta\tau_k\) の系列、\(E_d(k)\)
- 各段の三つ組 \((c_{A,j},\ \delta\tau_k/4,\ \delta\tau_{A,j})\)：各体の弦、共通の段、各体の固有時の増分 \(\delta\tau_{A,j}=(\delta\tau_k/4)(1-\tfrac12v_A^2-U_A)\)。これらの関係は結果として読む
- 内部時計 \(\omega_{{\rm int},A}=m_A/J_A\)（定数、Kerr–Newman のリングの回転）と系の角振動数の比、その連分数
- 放射量（反作用項またはモードへの仕事）、吸収量（モードから戻る仕事）、モードのエネルギー
- 1PN の重心関係から読む a, b の個別の位置、\(\mathbf L\)、\(\mathbf S_A\) の向き
- \(T_\tau\)

---

# 6. 初期化で行うこと

1. N（走らせる周数）から L_d に入れる項の次数を決める（§1）。
2. その次数までの全項を §3 の式で、weight の段ごとに \(L_d^{(s)}\) に入れ、§4.2 の陰的な方程式を小さい項について展開して固定の係数にする。B1・B2 では境界条件からモードの周波数・減衰・結合係数を決め、B3 では \(f\to1\) の極限が局所項に一致することを確認する。
3. 初期の離散円軌道：\(\Delta\phi\) 一定・弦一定で §4.2 を解き、\(L_0=2J\) から \(R_0\)、\(E_0\) を決める。\(\delta\tau_0\) は一周 = 62 段となる値。
4. 初期状態：\(Q=\mp3\)、\(m_a\)、\(m_b\)、\(J=\hbar/2\)（反平行／平行／0）、\(a_A=J_A/m_A\)、\(f\)、\(D\)、\(V_0=1\)、\(\sigma_0=1\)。
5. 相互作用は §4 の積和だけ。

---

# 参考文献

- A. Einstein, L. Infeld, B. Hoffmann, Ann. Math. 39, 65 (1938)
- L. Infeld, P. R. Wallace, Phys. Rev. 57, 797 (1940)
- L. D. Landau, E. M. Lifshitz, The Classical Theory of Fields, §65（Darwin Lagrangian）
- V. N. Golubenkov, Ya. A. Smorodinskii, JETP 4, 442 (1957)
- B. M. Barker, R. F. O'Connell, Phys. Rev. D 12, 329 (1975)（二体のスピン–軌道・スピン–スピン、電磁と重力）
- W. L. Burke, J. Math. Phys. 12, 401 (1971); K. S. Thorne, Astrophys. J. 158, 997 (1969)
- M. Khalil, A. Buonanno, J. Steinhoff, J. Vines, Phys. Rev. D 98, 104010 (2018)
- C. Kane, J. E. Marsden, M. Ortiz, J. Math. Phys. 40, 3353 (1999)
- J. E. Marsden, M. West, Acta Numerica 10, 357 (2001)
- 木原範昭, 反復交換散乱における有限位数共鳴の発見, Zenodo Concept DOI 10.5281/zenodo.21421366 (2026)
- 木原範昭, 厳密 Einstein–Maxwell アンカー aa/ab/bb/spin 結果 v1、厳密ルール再構築 spin hierarchy（2026-09-29、同リポジトリ）
