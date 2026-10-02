# Einstein–Maxwell に厳密に従う二体 a,b の動力学の導出
## ― 初期化で E–M を展開して固定係数を作り、相互作用中は積和一回だけで次の状態を作る ―

**著者:** 木原範昭 (Noriaki Kihara)  
**整理:** Claude  
**版:** v1  
**日付:** 2026-10-01  
**位置づけ:** 木原氏の指示「Einstein–Maxwell に厳密に従う」に基づく導出。相互作用仕様 v3.2・初期化仕様 v3.1・読出し仕様 v3.1、および固定点 v1〜v4 のうち本書と食い違う箇所は、本書を優先する。

---

# 0. E–M に厳密に従うと、a, b の運動はどこから出るか

Einstein–Maxwell 電磁真空

\[
G_{ab}=8\pi T^{\rm EM}_{ab},\qquad \nabla_aF^{ab}=0,\qquad dF=0
\]

には物質場がない。電荷と質量を持つ a, b は、場の方程式の解の中で**位相が閉じない場所**（Kerr–Newman 型の場の塊）としてだけ存在する。粒子を別に置かないという原則 2 は、E–M そのものの性質である。

その a, b がどう動くかも E–M が決める。Einstein–Infeld–Hoffmann（1938）は重力場の方程式だけから特異点の運動方程式を導いた。Infeld–Wallace（1940）は同じ方法を Maxwell 場に適用し、Lorentz 力と放射反作用を場の方程式から導いた。粒子の運動法則を別に仮定せず場の方程式から導くという Einstein 自身の方針であり、原則 2 と同じである。

この導出の形は、**場を初期配置の周りで展開し、各次数の係数を求める**ことである。展開の次数は精度要求で決まり、求めた係数は固定の多項式になる。これが「初期化で展開して、相互作用中は数式を固定する」ことの E–M 側の実体である。

空間格子の上で ℓ≤2 の波だけを進める案は、この展開を飛ばしていたので、a と b の結合（セルより小さい尺度のクーロン力・重力）がどこにも入っていなかった。

撤回された軌道の規則（吸収された放射の分だけ R が減る）が E–M と整合しない理由も同じである。E–M では放射反作用は**力の項**（電磁：Abraham–Lorentz–Dirac、重力：Burke–Thorne）として運動方程式に入り、エネルギー収支の等式は長時間平均でしか成り立たない。

---

# 1. 展開の小さい量と、必要な次数

展開の小さい量は、初期条件（\(Q=\mp3\)、\(J=\hbar/2\)、\(L_0=2J\)）が決める軌道速度 \(v_0=K/L_0\approx\alpha\) である。各項の相対的な大きさ（重力／クーロン \(=1/C_0\approx10^{-4}\) の場合）：

| 項 | 相対的な大きさ | 1 周期あたり |
|---|---|---|
| Newton 重力＋クーロン | 1 | — |
| Darwin（電磁の \(1/c^2\)）、EIH 1PN、スピン–軌道 | \(\alpha^2\) | \(5.3\times10^{-5}\) |
| 電磁の放射反作用（ALD） | \(\alpha^3\) | \(3.9\times10^{-7}\) |
| 電磁の \(1/c^4\)、重力 1PN × \(1/C_0\) | \(\alpha^4\)、\(\alpha^2/C_0\) | \(2.8\times10^{-9}\)、\(5\times10^{-9}\) |
| 重力波の放射反作用 | \(\alpha^5/C_0\) 程度 | \(10^{-11}\) |
| 電磁の \(1/c^6\) | \(\alpha^6\) | \(1.5\times10^{-13}\) |

相対的な大きさ \(\varepsilon\) の項を落としたときの位相の誤差は、N 周期で約 \(2\pi\varepsilon N\)。解像度は 1 段 \(=2\pi/124\) なので、条件は

\[
\boxed{\varepsilon\,N<\frac1{124}}
\]

- 放射で軌道が縮む run（B3、\(N\approx3\times10^5\) 周期）：\(\varepsilon<2.7\times10^{-8}\) → \(\alpha^4\) の項まで要る（\(\alpha^2\) では足りない）。
- ロックを読む run（\(2n=248\) 段 \(=2\) 周期）：\(\varepsilon<4\times10^{-3}\) → \(\alpha^2\) まででよい。

必要な次数は初期化で N から決め、その次数までの**すべての項**を入れる（項を選んで落とさない）。重力波の反作用は \(10^{-11}\) で解像度の下だが、E–M の項として存在するので係数に入れる（効かないことは読み出しで分かる）。

---

# 2. 各変分で保持する状態

| 変分 | E–M から出る項 | 状態 |
|---|---|---|
| \(\delta A\)（Maxwell） | クーロン、Darwin、\(1/c^4\)、放射反作用（ALD）、磁気双極子（\(g=2\)）のスピン–軌道・スピン–スピン | 相対位置 \(z=e^{\rho+i\phi}\) を \(e^{\rho}\)、\(e^{-\rho}\)、\(e^{i\phi}\) の三つの乗法状態で持ち、その有限履歴の倍音係数（速度・加速度・加速度の変化率は倍音係数の固定の線形結合） |
| \(\delta g\)（Einstein） | EIH 1PN、重力と電磁の交差項（Reissner–Nordström 型 \(Q^2/r^2\)）、Lense–Thirring 型スピン–軌道、Burke–Thorne 放射反作用 | 同じ状態（\(1/r^n\) は \(e^{-n\rho}\)、単位ベクトル \(n\) は \(e^{i\phi}\)） |
| スピン | 大きさ \(\hbar/2\)（定数）、向きの位相 | \(e^{i\psi_a}\)、\(e^{i\psi_b}\)（軌道面内の歳差）、軌道軸成分の符号 |
| 時間 | — | \(\sigma=i^k\)、\(\Theta\)、\(\delta\tau\)、\(D_\tau\)、\(T_\tau\) |
| 定数 | — | ONE、\(Q_a\)、\(Q_b\)、\(m_a\)、\(m_b\)、\(\alpha\)、\(J\)、\(f\)、\(N_{\rm tick}=31\)、\(\kappa\)、展開の次数 |

位置そのものではなく、履歴の倍音係数が状態である（第九思考実験の読み出し）。a, b の個別の位置は、相対位置と 1PN の重心関係から読み出す。

---

# 3. 入る項（\(G=c=1\)、Gauss。\(n=e^{i\phi}\)、\(r=e^{\rho}\)）

**0 次**

\[
L_0=\tfrac12m_av_a^2+\tfrac12m_bv_b^2+\frac Kr,\qquad K=m_am_b-Q_aQ_b
\]

**Darwin（Landau–Lifshitz §65）**

\[
L_D=-\frac{Q_aQ_b}{2r}\left[v_a\cdot v_b+(v_a\cdot n)(v_b\cdot n)\right]
\]

（符号は引力の K の定義に合わせる。）

**EIH 1PN**

\[
L_{1PN}=\tfrac18\left(m_av_a^4+m_bv_b^4\right)
+\frac{m_am_b}{2r}\left[3(v_a^2+v_b^2)-7v_a\cdot v_b-(v_a\cdot n)(v_b\cdot n)\right]
-\frac{m_am_b(m_a+m_b)}{2r^2}
\]

**重力と電磁の 1PN 交差項**（電磁場のエネルギーが重力源になる項、\(Q_aQ_b(m_a+m_b)/r^2\) 型）と**電磁の \(1/c^4\) 項**（Golubenkov–Smorodinskii）は、E–M から導かれた閉じた式として文献から取り、初期化で係数にする（Khalil–Buonanno–Steinhoff–Vines 2018 の Einstein–Maxwell 極限に 1PN の荷電二体 Lagrangian がまとまっている）。

**スピン–軌道（電磁、\(g=2\)、Thomas 歳差込み）**

\[
H_{SO}=\sum_A\frac{Q_B}{2m_A^2r^3}\,(g_A-1)\,\mathbf L\cdot\mathbf S_A
\qquad(g=2\ \text{で係数 }1)
\]

**スピン–スピン**

\[
H_{SS}=\frac1{r^3}\left[3(\boldsymbol\mu_a\cdot n)(\boldsymbol\mu_b\cdot n)-\boldsymbol\mu_a\cdot\boldsymbol\mu_b\right],
\qquad
\boldsymbol\mu_A=\frac{Q_A}{m_A}\mathbf S_A
\]

重力側の Lense–Thirring 項は \(1/C_0\) 倍で同じ形。

**放射反作用（電磁、ALD の最低次）**：各体に

\[
\mathbf F_A=\tfrac23Q_A^2\,\dot{\mathbf a}_A
\]

相対運動では \(\tfrac23(Q_a/m_a-Q_b/m_b)^2\mu\,\dot{\mathbf a}\) の形。\(\dot{\mathbf a}\) は履歴倍音の固定の線形結合。

**放射反作用（重力、Burke–Thorne）**

\[
\Phi_{RR}=\frac15x^ix^j\frac{d^5I_{ij}}{dt^5},\qquad I_{ij}\ \text{は質量四重極}
\]

5 階微分も履歴倍音から。

**吸収 B1〜B3**：Wheeler–Feynman（1945）により、放射反作用は遠方の吸収体の応答そのものである。吸収体がなければ（\(f=0\)）遅延と先進が半分ずつで正味の反作用は消え、完全吸収（\(f=1\)）で ALD が出る。したがって B1〜B3 は、**放射反作用の項に吸収率 \(f\) を掛ける**ことに対応する。これは E–M（Maxwell 方程式＋境界条件）の結果であり、軌道の規則ではない。

---

# 4. 時間の刻みと積和

- 1 刻み = 4 段。各段で一回の積和。\(\sigma'=i\sigma\)。
- 共通時計 \(\Theta\) は相対位置の位相 \(e^{i\phi}\) の回転（a, b が同じ角速度で重心を回る、系全体の回転）。\(\Theta'=e^{2\pi i/124}\Theta\)。
- \(\delta\tau\) は「1 周期 = 31 刻み」になるよう、\(\phi\) の履歴倍音から読む角速度の逆数で決める（逆数は補助状態、小さい相対変化の級数で更新）。
- 乗法状態の更新：

\[
e^{\rho'}=e^{\rho}\left(1+\Delta\rho+\tfrac12\Delta\rho^2+\cdots\right),\qquad
e^{i\phi'}=e^{i\phi}\left(1+i\Delta\phi-\tfrac12\Delta\phi^2+\cdots\right)
\]

\(\Delta\rho\)、\(\Delta\phi\) は 1 段あたり \(10^{-2}\) 程度なので、級数は固定の項数で倍精度の丸め以下。

- 履歴倍音の更新は 9/29 の厳密ルール再構築と同じ（位相シフト＋新しい標本の書き込み）。
- \(T_\tau'=D_\tau T_\tau\)、\(D_\tau'=D_\tau(1+\kappa\Delta\delta\tau/4+\cdots)\)。

すべて固定の係数と状態の積和である。\(1/r^n\)、\(\sqrt{\ }\)、\(\exp\)、\(\log\) は遷移の中に現れない。

---

# 5. 時間の歪みの読み出し

同じ展開から、各体の固有時の進み

\[
\frac{d\tau_A}{dt}=1-\tfrac12v_A^2-U_A
\]

（\(U_A\) は他方が作るポテンシャル、重力項とクーロン項の両方）が出る。運動による遅れと重力による遅れの両方が、E–M の展開の係数として入っている。

記録するもの：共通時計 \(\delta\tau\) に対する各体の時計の比 \(z_A=d\tau_A/dT\)、\(x^2+y^2+z^2-t^2\) の各刻みの値、内部時計 \(\omega_{{\rm int},A}=m_A/J_A\) との比とその連分数、放射量（ALD と Burke–Thorne の仕事）、吸収量（\(f\times\) 放射量）、\(\delta\tau\) の系列、\(T_\tau\)。

---

# 6. 初期化で行うこと

1. N（走らせる周期数）から必要な展開次数を決める（§1）。
2. その次数までの全項の係数を、§3 の E–M から導かれた式で数値化し、履歴倍音の線形結合係数と合わせて配列 M にする。
3. 初期状態：\(Q=\mp3\)、\(m_a\)、\(m_b\)、\(J=\hbar/2\)、\(L_0=2J\)、\(f\)、\(\delta\tau_0=2\pi/(31\Omega_0)\)。
4. 相互作用は \(z'=M\cdot(z\otimes z)\) だけ。

---

# 参考文献

- A. Einstein, L. Infeld, B. Hoffmann, Ann. Math. 39, 65 (1938)
- L. Infeld, P. R. Wallace, Phys. Rev. 57, 797 (1940)
- J. A. Wheeler, R. P. Feynman, Rev. Mod. Phys. 17, 157 (1945)
- L. D. Landau, E. M. Lifshitz, The Classical Theory of Fields, §65（Darwin Lagrangian）
- V. N. Golubenkov, Ya. A. Smorodinskii, JETP 4, 442 (1957)
- W. L. Burke, J. Math. Phys. 12, 401 (1971); K. S. Thorne, Astrophys. J. 158, 997 (1969)
- M. Khalil, A. Buonanno, J. Steinhoff, J. Vines, Phys. Rev. D 98, 104010 (2018)
