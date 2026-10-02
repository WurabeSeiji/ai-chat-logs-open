# Einstein–Maxwell に厳密に従う二体 a,b の動力学の導出 v2
## ― 初期化で E–M を展開して離散作用の係数を作り、相互作用中は離散変分方程式（位相刻みに対する傾斜）を固定の積和で評価する ―

**著者:** 木原範昭 (Noriaki Kihara)  
**整理:** Claude  
**版:** v2  
**日付:** 2026-10-01  
**置換対象:** `einstein_maxwell_strict_two_body_expansion_ja_v1_20261001.md`（v1 は残置）  
**位置づけ:** 木原氏の指示「Einstein–Maxwell に厳密に従う」に基づく導出。相互作用仕様 v3.2・初期化仕様 v3.1・読出し仕様 v3.1、固定点 v1〜v4 のうち本書と食い違う箇所は本書を優先する。

---

# 0A. v1 からの訂正

1. **式の誤り**：Darwin 項の符号、スピン–軌道項（\(Q_A\) と符号の欠落）、スピン–スピン項の符号、Burke–Thorne の四重極（跡なし）、ALD の形（各体の自己力だけでは相互の放射場の反作用が抜ける）。§3 で訂正。
2. **次数の表**：EIH 1PN は重力項の補正なので、主力（クーロン）に対して \(\alpha^2/C_0\)。§1 で訂正。
3. **B1〜B3 の扱い**：v1 は「放射反作用に吸収率 \(f\) を掛ける（Wheeler–Feynman）」と書いた。これは E–M の導出ではない。E–M に従うと、反射する境界からは往復時間だけ遅れて戻る場が a, b に働く。§3 末尾で訂正。
4. **時間の刻み**：v1 は E–M の展開を連続時間の運動方程式として扱い、\(\delta\tau_k=2\pi/(31\Omega_k)\) を毎刻み課していた。これは微分形への逆戻りであり、ロックを結果ではなく前提にしていた。v2 では離散作用を有限の位相刻みで定義し、増分についての変分（位相刻みに対する傾斜）と \(\delta\tau_k\) についての変分（離散エネルギー式）で次の状態と刻みを決める。31 刻み/周期は \(\delta\tau\) の初期値だけである。§4 を全面的に置き換え。
5. **\(\Theta\) の削除**：1 周期 = 124 段を課す \(\Theta'=e^{2\pi i/124}\Theta\) は削除。系の位相は \(e^{i\phi}\) から読む。
6. **内部時計**：§3 の項には内部時計 \(m_A/J_A\) の位相に結合する項がないことを明記（§5）。

---

# 0. E–M に厳密に従うと、a, b の運動はどこから出るか

Einstein–Maxwell 電磁真空

\[
G_{ab}=8\pi T^{\rm EM}_{ab},\qquad \nabla_aF^{ab}=0,\qquad dF=0
\]

には物質場がない。電荷と質量を持つ a, b は、場の方程式の解の中で**位相が閉じない場所**（Kerr–Newman 型の場の塊）としてだけ存在する。粒子を別に置かないという原則 2 は、E–M そのものの性質である。

a, b の運動も E–M が決める。Einstein–Infeld–Hoffmann（1938）は重力場の方程式だけから特異点の運動方程式を導き、Infeld–Wallace（1940）は同じ方法を Maxwell 場に適用して Lorentz 力と放射反作用を場の方程式から導いた。導出の形は、**場を初期配置の周りで展開し、各次数の係数を求める**ことである。これが「初期化で展開する」ことの E–M 側の実体であり、得られる項が離散作用 \(L_d\) の中身になる。

連続の軌道との一致は、事前に確かめる基準にしない。離散系がそれ自身で閉じた力学であり、結果を読み出して解釈する。

---

# 1. 展開の小さい量と、L_d に入れる項の次数

展開の小さい量は、初期条件（\(Q=\mp3\)、\(J=\hbar/2\)、\(L_0=2J\)）が決める軌道速度 \(v_0=K/L_0\approx\alpha\) である。各項の、主力（クーロン）に対する相対的な大きさ（重力／クーロン \(=1/C_0\approx10^{-4}\)）：

| 項 | 相対的な大きさ | 値 |
|---|---|---|
| Newton 重力＋クーロン | 1 | — |
| Darwin（電磁の \(1/c^2\)）、スピン–軌道（電磁） | \(\alpha^2\) | \(5.3\times10^{-5}\) |
| スピン–スピン（電磁） | \(\alpha^2\,m_a/(4m_b)\) | \(7\times10^{-9}\) |
| 電磁の放射反作用（ALD） | \(\tfrac23\alpha^3\) | \(2.6\times10^{-7}\) |
| 電磁の \(1/c^4\) | \(\alpha^4\) | \(2.8\times10^{-9}\) |
| EIH 1PN、重力–電磁の 1PN 交差項、Lense–Thirring | \(\alpha^2/C_0\) | \(5\times10^{-9}\) |
| 重力波の放射反作用（Burke–Thorne） | \(\alpha^5/C_0\) 程度 | \(10^{-11}\) |
| 電磁の \(1/c^6\) | \(\alpha^6\) | \(1.5\times10^{-13}\) |

相対的な大きさ \(\varepsilon\) の項を L_d から落としたときの位相のずれは N 周期で約 \(2\pi\varepsilon N\)。読み出しの解像度 1 段 \(=2\pi/124\) に対して

\[
\boxed{\varepsilon\,N<\frac1{124}}
\]

を満たす次数まで、**その次数の全項**を入れる。放射で軌道が縮む run（B3、\(N\approx3\times10^5\) 周期）では \(\alpha^4\) まで、ロックを読む run（2 周期）では \(\alpha^2\) まで。重力波の反作用は解像度の下だが、E–M の項として入れる（効かないことは読み出しで分かる）。

---

# 2. 状態

| 区分 | 状態 | 表現 |
|---|---|---|
| 相対位置 | \(z=e^{\rho+i\phi}\) | 乗法状態 \(e^{\rho}\)、\(e^{-\rho}\)、\(e^{i\phi}\)（\(1/r^n=e^{-n\rho}\)、単位ベクトル \(n=e^{i\phi}\)） |
| 相対位置の履歴 | 直近の刻みの \(z\) の倍音係数 | 速度・加速度・加速度の変化率は倍音係数の固定の線形結合（第九思考実験） |
| 放射の履歴（B1・B2） | 双極子 \(d\)・四重極 \(I_{ij}\) の、境界までの往復時間ぶんの倍音係数 | 戻る場の評価に使う（§3） |
| スピン | 大きさ \(\hbar/2\)（定数）、向き | \(e^{i\psi_a}\)、\(e^{i\psi_b}\)、軌道軸成分の符号 |
| 時間 | \(\delta\tau_k\)（変分で決まる）、\(1/\delta\tau_k\)（補助）、\(T_\tau\)、\(D_\tau\)、段の位相 \(\sigma=i^s\) | \(T_\tau'=D_\tau T_\tau\)、\(D_\tau\) は \(\delta\tau\) の相対変化の級数で更新 |
| 定数 | ONE、\(Q_a\)、\(Q_b\)、\(m_a\)、\(m_b\)、\(\alpha\)、\(J\)、\(f\)、境界距離 \(D\)、\(\kappa\)、展開の次数 | 恒等 |

a, b の個別の位置は、相対位置と 1PN の重心関係から読み出す。

---

# 3. L_d に入る項（\(G=c=1\)、Gauss。\(n=e^{i\phi}\)、\(r=e^{\rho}\)）

**0 次**

\[
L_0=\tfrac12m_av_a^2+\tfrac12m_bv_b^2+\frac Kr,\qquad K=m_am_b-Q_aQ_b
\]

**Darwin（Landau–Lifshitz §65）**

\[
L_D=+\frac{Q_aQ_b}{2r}\left[v_a\cdot v_b+(v_a\cdot n)(v_b\cdot n)\right]
\]

（クーロン項 \(-Q_aQ_b/r\) と逆符号。）

**EIH 1PN**

\[
L_{1PN}=\tfrac18\left(m_av_a^4+m_bv_b^4\right)
+\frac{m_am_b}{2r}\left[3(v_a^2+v_b^2)-7v_a\cdot v_b-(v_a\cdot n)(v_b\cdot n)\right]
-\frac{m_am_b(m_a+m_b)}{2r^2}
\]

**重力–電磁の 1PN 交差項**（電磁場のエネルギーが重力源になる項、\(Q_aQ_b(m_a+m_b)/r^2\) 型）と**電磁の \(1/c^4\) 項**（Golubenkov–Smorodinskii）は、E–M から導かれた閉じた式として文献から取り、初期化で係数にする（Khalil–Buonanno–Steinhoff–Vines 2018 の Einstein–Maxwell 極限）。

**スピン–軌道（電磁、\(g=2\)、Thomas 歳差込み）**

\[
H_{SO}=-\sum_{A}\frac{Q_AQ_B}{2m_A^2r^3}\,(g_A-1)\,\mathbf L\cdot\mathbf S_A
\qquad(g=2\ \text{で }g-1=1)
\]

（水素型 \(Q_aQ_b=-9\) では \(+9/(2m_A^2r^3)\,\mathbf L\cdot\mathbf S_A\)。）

**スピン–スピン**

\[
H_{SS}=\frac1{r^3}\left[\boldsymbol\mu_a\cdot\boldsymbol\mu_b-3(\boldsymbol\mu_a\cdot n)(\boldsymbol\mu_b\cdot n)\right],
\qquad
\boldsymbol\mu_A=\frac{Q_A}{m_A}\mathbf S_A
\]

重力側の Lense–Thirring 項は \(1/C_0\) 倍で同じ形。

**放射反作用（電磁、双極子次数）**：相対座標に対する力

\[
\mathbf F_{RR}^{\rm EM}=\tfrac23\,\mu\left(\frac{Q_a}{m_a}-\frac{Q_b}{m_b}\right)^2\dot{\mathbf a}
\]

（各体の自己力 \(\tfrac23Q_A^2\dot{\mathbf a}_A\) と、相手の放射場の反作用の和。\(\dot{\mathbf a}\) は履歴倍音の固定の線形結合。）

**放射反作用（重力、Burke–Thorne）**

\[
\Phi_{RR}=\frac15x^ix^j\frac{d^5I_{ij}}{dt^5},\qquad I_{ij}\ \text{は跡なし質量四重極}
\]

**吸収 B1〜B3（E–M の境界条件として）**：外へ出る放射に対する反作用は上の二つの項である。境界（距離 \(D\)）が反射すると、放射した場は往復時間 \(2D\) だけ遅れて a, b に戻って働く。戻る場の項は、放射反作用の項を往復時間だけ過去の履歴で評価し、反射率 \((1-f)\) と反射の位相因子を掛けた形になる。係数（位相因子）は境界条件から初期化で決める。

- B1（\(f=0\)）：全反射。戻る項が全部入る。往復時間ぶんの履歴が必要。
- B2（\(0<f<1\)）：部分反射。戻る項に \((1-f)\)。
- B3（\(f=1\)）：全吸収。戻る項なし。放射反作用だけ。

---

# 4. 離散作用と変分（相互作用の全体）

## 4.1 離散作用

刻み k の状態を \(z_k\)（相対位置）、\(\psi_{A,k}\)（スピンの向き）、刻みの長さを \(\delta\tau_k\) とし、離散作用を

\[
S_d=\sum_k L_d(z_k,z_{k+1},\psi_k,\psi_{k+1},\delta\tau_k)
\]

で定義する。\(L_d\) は有限の増分で書く。相対運動の運動項は 1 刻みの変位（弦）

\[
c_k^2=|z_{k+1}-z_k|^2=e^{2\rho_k}+e^{2\rho_{k+1}}-2e^{\rho_k+\rho_{k+1}}\cos\Delta\phi_k
\]

を使って \(\mu c_k^2/(2\delta\tau_k)\)。§3 の相互作用項は、刻みの両端の状態の対称な平均で評価する（時間反転対称を保つ選択。これは E–M ではなく離散作用の定義の選択であり、記録する）。

1 刻みは 4 段に分け、各段 \(s=0..3\) は小増分 \(\delta\tau_k/4\) で、段の位相 \(\sigma_s=i^s\) を状態として持つ。以下の変分方程式は段ごとに成り立つ。

## 4.2 増分についての変分（位相刻みに対する傾斜）

\(z_k\) についての変分から離散 Euler–Lagrange 方程式

\[
D_2L_d(z_{k-1},z_k,\delta\tau_{k-1})+D_1L_d(z_k,z_{k+1},\delta\tau_k)=F^{-}_d(k-1)+F^{+}_d(k)
\]

を得る。右辺は放射反作用と戻る場の離散 Lagrange–d'Alembert の力（Marsden–West 2001 §3）。\(D_1L_d\)、\(D_2L_d\) が「位相刻みに対する傾斜」（離散運動量）である。運動項については

\[
\mu\frac{z_k-z_{k-1}}{\delta\tau_{k-1}}-\mu\frac{z_{k+1}-z_k}{\delta\tau_k}
\]

であり、\(1/\delta\tau\) は補助状態として持つ。この式を \(z_{k+1}\) について解いた形が次の状態を作る積和で、係数は固定である。スピンの向き \(\psi_{A}\) についても同じ形。

## 4.3 \(\delta\tau_k\) についての変分（離散エネルギー式）

\(\delta\tau_k\) を変分の変数とすると（Kane–Marsden–Ortiz 1999）

\[
\frac{\partial L_d}{\partial\delta\tau_k}=0
\quad\Longrightarrow\quad
\frac{\mu c_k^2}{2\delta\tau_k^2}+V_d(k)=E_k
\]

が離散のエネルギー式になり、\(\delta\tau_k\) を決める。放射があるときは離散 Lagrange–d'Alembert により \(E_{k+1}-E_k\) が放射反作用（と戻る場）の離散の仕事に等しく、\(\delta\tau_k\) はそれに従って変わる。実装では \(\delta\tau_k/\delta\tau_{k-1}\) を、ほぼ 1 に近い比の級数（固定の項数）で更新する。

**\(\delta\tau\) は設定しない。** 初期値 \(\delta\tau_0\) だけを、初期軌道で 1 周期 = 31 刻みとなる値に置く。以後は上の式で決まる。1 周期が何刻みになるか（31 のまま保つか、\(23/124\) に対応する比に落ちるか）は結果であり、読み出す。

## 4.4 時間

\[
T_\tau'=D_\tau T_\tau,\qquad D_\tau'=D_\tau\left(1+\kappa\,\Delta\delta\tau/4+\tfrac12(\kappa\,\Delta\delta\tau/4)^2+\cdots\right)
\]

乗法状態の更新：\(e^{\rho'}=e^{\rho}(1+\Delta\rho+\tfrac12\Delta\rho^2+\cdots)\)、\(e^{i\phi'}=e^{i\phi}(1+i\Delta\phi-\tfrac12\Delta\phi^2+\cdots)\)。級数は固定の項数で倍精度の丸め以下。

遷移の中に \(\sqrt{\ }\)、\(\exp\)、\(\log\)、状態による割り算は現れない。

---

# 5. 読み出し（発展へ帰還させない）

- 系の位相 \(\phi_k\)、1 周期あたりの刻み数、\(\delta\tau_k\) の系列
- 各刻みの三つ組 \((c_{A,k},\ \delta\tau_k,\ \delta\tau_{A,k})\)：各体の弦、共通の刻み、各体の固有時の増分 \(\delta\tau_{A,k}=\delta\tau_k(1-\tfrac12v_A^2-U_A)\)（\(U_A\) は他方の重力＋クーロンのポテンシャル。運動と重力の両方の時間の歪みが同じ展開から出る）。これらの間の関係は結果として読む
- 内部時計 \(\omega_{{\rm int},A}=m_A/J_A\)（定数）と系の角振動数の比、その連分数。§3 の項には \(\omega_{\rm int}\) の位相に結合する項がないので、この比は「定数と \(\Omega\) の比」として読む
- 放射量（放射反作用の項の仕事）、吸収量（戻る場の項の仕事）、エネルギー \(E_k\)
- \(T_\tau\)

---

# 6. 初期化で行うこと

1. N（走らせる周期数）から L_d に入れる項の次数を決める（§1）。
2. その次数までの全項を §3 の式で \(L_d\) に入れ、離散 E–L 方程式と離散エネルギー式を §4 の形に展開して、固定の係数（履歴倍音の線形結合係数を含む）を数値化する。境界条件から戻る場の項の係数を決める。
3. 初期状態：\(Q=\mp3\)、\(m_a\)、\(m_b\)、\(J=\hbar/2\)（反平行／平行／0）、\(L_0=2J\)、\(f\)、\(D\)、\(\delta\tau_0\)（初期軌道で 31 刻み/周期）。
4. 相互作用は §4 の積和だけ。

---

# 参考文献

- A. Einstein, L. Infeld, B. Hoffmann, Ann. Math. 39, 65 (1938)
- L. Infeld, P. R. Wallace, Phys. Rev. 57, 797 (1940)
- L. D. Landau, E. M. Lifshitz, The Classical Theory of Fields, §65（Darwin Lagrangian）
- V. N. Golubenkov, Ya. A. Smorodinskii, JETP 4, 442 (1957)
- W. L. Burke, J. Math. Phys. 12, 401 (1971); K. S. Thorne, Astrophys. J. 158, 997 (1969)
- M. Khalil, A. Buonanno, J. Steinhoff, J. Vines, Phys. Rev. D 98, 104010 (2018)
- C. Kane, J. E. Marsden, M. Ortiz, J. Math. Phys. 40, 3353 (1999)（時間刻みを変分変数とする離散力学）
- J. E. Marsden, M. West, Acta Numerica 10, 357 (2001)（離散変分力学、離散 Lagrange–d'Alembert）
