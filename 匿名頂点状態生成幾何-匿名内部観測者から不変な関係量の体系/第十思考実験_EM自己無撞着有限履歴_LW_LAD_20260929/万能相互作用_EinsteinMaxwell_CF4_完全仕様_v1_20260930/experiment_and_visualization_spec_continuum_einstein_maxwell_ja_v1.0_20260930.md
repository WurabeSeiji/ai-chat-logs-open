# Einstein–Maxwell 万能相互作用：連続天体系 第1実験・図化仕様
## ― 古典連続系ベースライン、放射・吸収収支、準定常／継続縮退判定、将来の離散谷検出に備えた無損失記録 ―

**著者:** 木原範昭 (Noriaki Kihara)  
**仕様版:** v1.0  
**日付:** 2026-09-30  
**親仕様:** `universal_interaction_einstein_maxwell_cf4_full_spec_ja_v1.0_20260930.md`

---

# 0. 本仕様の位置づけ

本仕様は、親仕様で定義した

\[
\boxed{\psi_{n+1}=S(\psi_n)\psi_n}
\]

を実際に走らせるための**第1段階の実験計画と図化計画**を定める。

第1段階では、系を**天体スケールの古典・連続 Einstein–Maxwell 系**として扱う。

ここではまだ、

- \(m/n\) の有理比探索
- 離散電荷番号による谷探索
- \(\Delta\tau\) の物理的量子化探索
- 「粒子準位」の同定
- \(m/n\) に対する \(\cos^2\) 型反復谷の直接探索

は行わない。

ただし、本研究では後に

\[
\frac{m}{n}
\]

型の演算が反復されたとき、非常に狭く深い谷が現れる可能性を重視している。

そのため第1段階から、

1. raw data を一切平滑化せず保存する
2. 正負を失わない差分を保存する
3. 非常に小さい値を見られる symlog / log 表示を用意する
4. 初期値を丸めず保持する
5. \(\Delta\tau\)、cutoff、基底、抽出半径などを run metadata と状態に保持する
6. 途中でデータを間引いて不可逆に捨てない
7. 後から \(m/n\) を再構成できるだけの初期・演算情報を残す

ことを必須とする。

つまり、

\[
\boxed{
\text{第1実験では量子化を探さないが、量子化を後から見落とさない記録を行う}
}
\]

ことが本仕様の基本方針である。

---

# 1. 第1実験の主目的

主目的は、full nonlinear Einstein–Maxwell を単一の万能演算

\[
K(\psi)\rightarrow S(\psi)
\]

で発展させたとき、荷電二体系が古典連続系としてどの長時間挙動を示すかを調べることである。

特に、

\[
\boxed{
T\mapsto
\left(
r_{\rm orb}(T),
\Delta F_{\rm GW}(T),
\Delta F_{\rm EM}(T)
\right)
}
\]

を中心読出しとする。

ここで期待される大分類は次の通りである。

### A. 継続縮退
放射・吸収を含むエネルギー散逸が続き、

\[
r_{\rm orb}(T)\downarrow
\]

が止まらず、最終的に merger / common horizon / \(r\to0\) 相当へ向かう。

### B. 有限半径準定常
初期過渡後に

\[
\frac{dr_{\rm orb}}{dT}\to0
\]

となり、

\[
r_{\rm orb}(T)\to r_*>0
\]

へ近づく。

### C. 準定常へ向かうが弱い放射が残る
有限半径へ非常にゆっくり近づきながら、GW/EM の弱い flux が長時間残る。

### D. 膨張・散乱
電磁反発等が優勢となり、

\[
r_{\rm orb}(T)\uparrow
\]

となる。

### E. 非定常・振動状態
\(r_{\rm orb}\)、GW flux、EM flux が有限振幅で長時間振動し、単純な plateau に入らない。

第1実験では上記のどれが現れるかを**結果として判定する**。事前に特定の結論を仮定しない。

---

# 2. スケールに対する実験上の立場

本研究モデルでは、無次元化した状態方程式と万能相互作用を用いるため、実験では絶対的な天体サイズそのものよりも、

\[
\frac{Q}{M},\quad
\frac{r}{M},\quad
\frac{T}{M},\quad
\Omega M
\]

などの無次元比を主要読出しとする。

第1実験は**天体系として解釈**して行う。

ただし、結果保存時には raw dimensionless data を残し、後に別スケールへ写像できる形式とする。

第1段階で「素粒子状態が得られた」とは解釈しない。

---

# 3. 第1実験では行わないこと

以下は禁止する。

1. \(m/n\) を横軸にした探索
2. 有理数候補を意図的に密集させた parameter scan
3. 電荷番号 \(n\) の整数周期探索
4. \(\Delta\tau=m\Delta\tau_0\) の物理的準位探索
5. 既知の谷位置に合わせた adaptive sampling
6. 結果を \(\cos^2(\pi m/n)\) に fit すること
7. smooth / interpolation により raw extrema を消すこと
8. 量子化を期待した post-selection

これらは連続系ベースライン完成後の別実験とする。

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
\epsilon_{\exp},\quad
m_{\exp,\max}
\]

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
Q_over_M_decimal = 0.3
Q_over_M_exact = 3/10
```

ただし第1実験では `3/10` を「有理比として解析」しない。

---

# 5. 状態変数として追加・確認すべき実験設定

親仕様の原則

\[
\boxed{\text{初期化するものは原則すべて状態}}
\]

に従い、本実験で初期設定する値も \(\psi\) に恒等状態として保持する。

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
r_{\rm common-horizon}
\]

これらは実験runの停止条件であり、初期化時に決める。

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

ただし固定絶対 tolerance ではなく、後述する convergence floor から判定閾値を作る。

---

# 6. 第1実験の初期物理条件

最初は複雑性を抑えるため、

\[
\boxed{\mathbf J_a^{(0)}=\mathbf J_b^{(0)}=0}
\]

とする。

spin は current state の独立レジスタとして発展させず、必要な場合は horizon geometry から readout する。

質量比の基準はまず

\[
\boxed{M_a=M_b}
\]

とする。

その後、質量比を変えたrunを追加する。

---

# 7. 電荷条件：連続古典系としてのrun family

正規化として

\[
q_a=\frac{Q_a}{M_a},
\qquad
q_b=\frac{Q_b}{M_b}
\]

を用いる。

regular nonspinning charged black hole として初期化するrunでは原則

\[
|q_a|\le1,\qquad |q_b|\le1
\]

を初期化ゲートとする。

第1段階では次のfamilyを使う。

## Family N: 無電荷基準

\[
q_a=q_b=0
\]

Einstein 真空連星に近い基準。

## Family R: 同符号・反発側

\[
q_a=q_b=q>0
\]

連続的に

\[
q=0.3,\ 0.6,\ 0.9
\]

を基準点とする。

数値安定性を確認後、

\[
q=0.95,\ 0.97,\ 0.99
\]

を追加し、古典的重力–電気反発均衡へ近づく挙動を見る。

### extremal anchor

\[
q=1
\]

は通常のgeneric binary runとは分け、Majumdar–Papapetrou 静的平衡解に対応する validation anchor として扱う。

初期化器が extremal horizon を安定に生成できない場合は、`SKIP_INITIALIZER_LIMIT` と記録し、近傍 \(q<1\) のみ走らせる。

## Family A: 逆符号・引力強化側

\[
q_a=+q,\qquad
q_b=-q
\]

\[
q=0.3,\ 0.6,\ 0.9
\]

を基準とする。

重力と電磁引力が同方向となり、同符号系との比較対象となる。

## Family D: 電荷非対称

\[
q_a\neq q_b
\]

例：

\[
(q_a,q_b)=(0.9,0.3),(0.9,0),(0.6,0.3)
\]

とする。

これは equal charge-to-mass 系では抑制される可能性がある電磁 dipole radiation を分離するための重要な比較群とする。

---

# 8. 初期軌道条件

第1runでは quasi-circular initial data を使用する。

初期化器は必ず次を report する。

\[
d_0
=
\left|
\mathbf x_a^{(0)}-\mathbf x_b^{(0)}
\right|
\]

equal-mass の場合は表示用軌道半径を

\[
\boxed{
r_{\rm orb}^{(0)}=\frac{d_0}{2}
}
\]

と定義する。

旧実験の `P0=50` との対応を使用する場合は、

\[
P_0\longleftrightarrow d_0
\]

の写像を initializer metadata に明示する。

**旧 `P0` から座標分離を暗黙に推定してはならない。**

第1基準runでは1つの十分広い初期分離を採用し、その後、

\[
d_0/M_{\rm tot}
\]

を複数値に変えて初期距離依存性を調べる。

---

# 9. Δτ の扱い

第1実験では \(\Delta\tau\) を**物理的離散準位として探索しない**。

\(\Delta\tau\) は状態に保持するが、第1実験での変更目的は純粋に numerical convergence test とする。

同一の物理初期条件について、

\[
h,\quad\frac h2,\quad\frac h4
\]

の3系列を走らせる。

ここで

\[
h=\Delta\tau
\]

である。

判定：

- 物理的readoutが \(h\to0\) で収束する
- 小さな構造が \(\Delta\tau\) に同期して移動するなら数値artifactを疑う
- 後の量子化実験に進む前に continuum baseline を確立する

ことを必須とする。

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

深い谷状のstructureが後に現れたとき、

- 物理
- timestep artifact
- spectral cutoff artifact

を区別するため、最初から収束データを残す。

---

# 11. 時間 readout

状態に

\[
T_\tau=e^{\kappa_\tau\tau}
\]

がある場合、図の横軸に使う物理時間readoutを

\[
\boxed{
T
\equiv
\mathcal R_T(\psi)
=
\frac{\log T_\tau}{\kappa_\tau}
}
\]

と定義する。

raw state としては

\[
T_\tau
\]

も保存する。

したがってCSVには

```text
T_state
T_read
```

を両方保存する。

以後の図の `T` は `T_read` を意味する。

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
- horizon-centroid separation
- areal-radius based separation surrogate

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

## 13.2 horizon absorption

各 horizon について

\[
F_{{\rm GW},a}^{H}(T),
\qquad
F_{{\rm GW},b}^{H}(T)
\]

を読み出し、

\[
\boxed{
F_{\rm GW}^{\rm abs}
=
F_{{\rm GW},a}^{H}
+
F_{{\rm GW},b}^{H}
}
\]

とする。

ここで `abs` は horizon へ流入する正のenergy flux magnitude とする。

---

# 14. EM readout

## 14.1 外向きEM

各 extraction sphere でPoynting fluxを積分し、

\[
F_{\rm EM}^{\rm out}(T;R_{{\rm ext},k})
\]

を得る。

可能なら electromagnetic Newman–Penrose radiative scalars も保存する。

## 14.2 horizon absorption

各 horizon で

\[
F_{{\rm EM},a}^{H}(T),
\qquad
F_{{\rm EM},b}^{H}(T)
\]

を読み、

\[
\boxed{
F_{\rm EM}^{\rm abs}
=
F_{{\rm EM},a}^{H}
+
F_{{\rm EM},b}^{H}
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

標準的なblack-hole energy balanceでは、

- 無限遠へ出ていくenergy
- horizonへ吸収されるenergy

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

# 18. horizon readout

各stepまたは指定readout cadenceで、

\[
M_a^{\rm read}(T),
\quad
M_b^{\rm read}(T)
\]

\[
Q_a^{\rm read}(T),
\quad
Q_b^{\rm read}(T)
\]

\[
\mathbf J_a^{\rm read}(T),
\quad
\mathbf J_b^{\rm read}(T)
\]

\[
A_{H,a}(T),\quad A_{H,b}(T)
\]

を保存する。

初期spinが0でも、

\[
J^{\rm read}(T)
\]

は保存する。

horizon absorption と

\[
\Delta M_H,\quad
\Delta J_H,\quad
\Delta Q_H
\]

の対応を監査できるようにする。

---

# 19. constraint / numerical readout

各stepで最低限、

\[
\|{\cal C}_{\rm GH}\|
\]

\[
\|D_iE^i\|
\]

\[
\|D_iB^i\|
\]

を保存する。

可能なら、

- \(L_2\)
- \(L_\infty\)

の両方を保存する。

さらに、

\[
\epsilon_{K\psi}
=
\|K(\psi)\psi-F_{\rm direct}(\psi)\|
\]

も実装監査として保存する。

---

# 20. 第1主図：ユーザー指定3段図

**Figure 01: Orbit radius and radiation–absorption residuals**

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

# 26. Figure 06：horizon量

最低限、

\[
M_a^{\rm read},M_b^{\rm read}
\]

\[
Q_a^{\rm read},Q_b^{\rm read}
\]

\[
|J_a^{\rm read}|,|J_b^{\rm read}|
\]

を共通時間軸で図化する。

放射吸収との対応を確認する。

---

# 27. Figure 07：Δτ convergence

同じ物理初期条件について、

\[
h,\quad h/2,\quad h/4
\]

を重ねる。

最低限：

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

重要：

第1実験ではこの図を**量子化図と解釈しない**。

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

# 30. Figure 10：charge continuum family

第1実験終了後、古典的parameter依存を見るため、

横軸：

\[
q=\frac{Q}{M}
\]

縦軸候補：

\[
r_{\rm final}/M
\]

\[
\langle|\dot r|\rangle_{\rm late}
\]

\[
\langle|\Delta F_{\rm GW}|\rangle_{\rm late}
\]

\[
\langle|\Delta F_{\rm EM}|\rangle_{\rm late}
\]

\[
T_{\rm merger}
\]

をscatterとして出す。

**線で滑らかに補間しない。**

全run点を表示する。

第1段階では \(m/n\) 軸へ変換しない。

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

# 32. 準定常判定は数値誤差床を基準にする

固定した任意値

\[
|\dot r|<10^{-k}
\]

だけで準定常と判定しない。

まず、

- \(\Delta\tau\) convergence
- \(N_{\rm cut}\) convergence
- extraction-radius convergence
- constraint residual

から numerical floor を推定する。

\[
\epsilon_{\rm num}
=
\max(
\epsilon_{\Delta\tau},
\epsilon_N,
\epsilon_{\rm ext},
\epsilon_{\rm constraint}
)
\]

とする。

late-time window で

\[
\left|\frac{d\log r}{dT}\right|
<
C_r\epsilon_{\rm num}
\]

が複数軌道にわたって続く場合を plateau candidate とする。

ただしこれは**readout分類**であり、発展へ帰還させない。

---

# 33. outcome分類

各runについて、最終的に次のラベルを出す。

## `INSPIRAL`

\[
\dot r<0
\]

が長時間継続し、plateau条件を満たさない。

## `MERGER`

common horizon または merger判定。

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

解像度収束が得られず物理解釈しない。

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

万能演算が自明状態を壊さないことを確認。

## Stage V1: single RN

\[
q=0,\ 0.3,\ 0.9
\]

でstationary test。

## Stage V2: static balance anchor

可能ならMajumdar–Papapetrouの

\[
q=1
\]

2-center static configuration。

これは重力–電気反発の古典的正確な均衡anchorとして使用する。

## Stage B0: neutral binary

\[
q_a=q_b=0
\]

## Stage B1: same-sign charged binary

\[
q_a=q_b=0.3
\]

\[
q_a=q_b=0.6
\]

\[
q_a=q_b=0.9
\]

## Stage B2: near-balance continuum

\[
q=0.95,\ 0.97,\ 0.99
\]

## Stage B3: opposite-sign

\[
(q_a,q_b)=(+0.3,-0.3),(+0.6,-0.6),(+0.9,-0.9)
\]

## Stage B4: asymmetric-charge

\[
(0.9,0.3),(0.9,0),(0.6,0.3)
\]

## Stage C: convergence

代表3caseについて、

- \(\Delta\tau\) 3水準
- \(N_{\rm cut}\) 3水準
- extraction radius複数

を走らせる。

---

# 36. 第1実験での「均衡」の意味

本実験で重力–電磁均衡を見るとき、単一のNewtonian force equalityだけを使わない。

少なくとも、

1. orbital radius drift
2. orbital angular frequency drift
3. GW outward flux
4. GW horizon absorption
5. EM outward flux
6. EM horizon absorption
7. horizon mass/spin/charge change
8. constraints

を同時に見る。

特に classical Einstein–Maxwell には extremal same-sign charged black holes が静的に均衡する Majumdar–Papapetrou class が存在するため、これは第1連続系実験の重要なanchorになる。

---

# 37. raw time-series file schema

最低限、1 row / readout step とする。

```text
run_id
step
T_state
T_read
dt_state
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
gw_flux_abs_a
gw_flux_abs_b
gw_flux_abs_total
gw_out_minus_abs
gw_out_plus_abs

em_flux_out
em_flux_abs_a
em_flux_abs_b
em_flux_abs_total
em_out_minus_abs
em_out_plus_abs

total_positive_drain

Ma_read
Mb_read
Qa_read
Qb_read
Jax_read
Jay_read
Jaz_read
Jbx_read
Jby_read
Jbz_read

horizon_area_a
horizon_area_b

gh_constraint_L2
gh_constraint_Linf
maxwell_divE_L2
maxwell_divE_Linf
maxwell_divB_L2
maxwell_divB_Linf

Kpsi_rhs_residual_L2
Kpsi_rhs_residual_Linf
```

---

# 38. waveform/raw-mode file schema

GW:

```text
T
R_ext
ell
m
Re_Psi4
Im_Psi4
```

EM:

```text
T
R_ext
ell
m
Re_EM_mode
Im_EM_mode
```

可能ならbinary format/HDF5を用い、CSVはsummary用とする。

---

# 39. run manifest

各runにつき、

```text
run_manifest.json
```

を作り、最低限、

```text
spec_version
code_version
state_schema_version
initializer_version
S_generator_version
CF4_version

all_initial_state_values
exact_parameter_forms
basis_definition
domain_definition
extraction_radii
dt
N_cut
solver_tolerances
stopping_conditions
readout_cadence
checkpoint_cadence
```

を保存する。

---

# 40. 推奨出力ファイル構成

```text
run_<RUNID>/
    run_manifest.json
    initial_state.*
    final_state.*
    timeseries_raw.*
    horizon_readout.*
    gw_modes.*
    em_modes.*
    energy_budget.*
    constraint_history.*

    figures/
        figure01_orbit_radius_gw_em_balance_linear.png
        figure01_orbit_radius_gw_em_balance_linear.svg
        figure01_orbit_radius_gw_em_balance_symlog.png
        figure01_orbit_radius_gw_em_balance_symlog.svg
        figure02_orbit_and_positive_energy_drain.png
        figure02_orbit_and_positive_energy_drain.svg
        figure03_gw_em_out_abs_decomposition_linear.png
        figure03_gw_em_out_abs_decomposition_log.png
        figure04_orbit_decay_rate.png
        figure05_phase_omega.png
        figure06_horizon_quantities.png
```

---

# 41. 図化の共通規則

1. 横軸は原則 `T_read`
2. raw data を基準にする
3. smoothingした場合は必ず別series
4. raw seriesを消さない
5. flux differenceはzero lineを描く
6. signed quantityに通常logを使わない
7. signed residualにはsymlogを使う
8. positive fluxにはlog版を追加可
9. 軸単位とnormalizationを図中に明記
10. run IDをcaptionかfilenameに含める
11. transient除外を行った場合、除外区間を明示
12. late-time averageだけを見せずfull historyも残す
13. interpolation curveだけの図は禁止
14. parameter sweepは全raw pointを表示する
15. y-rangeを手動で狭めて変動を誇張しない
16. 逆に極小構造を潰さないためsymlog補助図を作る

---

# 42. 第1実験完了条件

第1連続系実験は以下を満たしたとき完了とする。

## 数値側

\[
K(\psi)\psi=F_{\rm direct}(\psi)
\]

が設定精度内。

\[
\Delta\tau\to0
\]

で主要readoutが収束。

\[
N_{\rm cut}\uparrow
\]

で主要readoutが収束。

複数 extraction radius でradiative fluxが整合。

constraint が許容範囲。

## 物理readout側

少なくとも、

\[
r_{\rm orb}(T)
\]

\[
F_{\rm GW}^{\rm out},
F_{\rm GW}^{\rm abs},
\Delta F_{\rm GW},
D_{\rm GW}
\]

\[
F_{\rm EM}^{\rm out},
F_{\rm EM}^{\rm abs},
\Delta F_{\rm EM},
D_{\rm EM}
\]

\[
M_{a,b}^{\rm read},
Q_{a,b}^{\rm read},
J_{a,b}^{\rm read}
\]

が全runで取得可能。

## 結果側

各charge familyについて、

- inspiral
- merger
- quasi-stationary
- radiative quasi-stationary
- expanding
- oscillatory
- unresolved

のいずれかに、数値収束を伴って分類できる。

---

# 43. 第2段階へ進む条件

以下を満たすまでは \(m/n\) 探索へ進まない。

1. continuum run の収束確認
2. \(\Delta\tau\) を数値刻みとして小さくした極限の確認
3. charge continuum のclassical behaviorの把握
4. Majumdar–Papapetrou balance anchor の再現またはinitializer限界の明示
5. flux readoutのenergy accounting確認
6. raw data保存系の確立

その後初めて、

\[
m/n
\]

型の内部演算反復による深い谷が、continuum baseline からどのように分岐するかを別実験として調べる。

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

# 45. 第1実験で問う最終質問

第1段階で問うのは量子化そのものではない。

\[
\boxed{
\text{古典連続 Einstein–Maxwell 二体系は、万能相互作用 }S(\psi)\text{ の下で、どの条件なら放射し続け、どの条件なら有限半径の長寿命状態へ近づくか？}
}
\]

である。

その答えを完全なraw time seriesとして残す。

その上で後の段階において、

\[
\boxed{
\text{離散化された電荷・時間・反復比 }m/n\text{ が、この連続系のどこに深い選択谷を作るか}
}
\]

を比較できるようにする。

---

# 46. 参考になる標準物理上のアンカー

- Generalized harmonic numerical relativity では \(\Psi_4\) を用いたGW readoutが標準的に使用できる。
- black-hole horizon へのGW吸収は horizon energy / angular-momentum flux として扱われ、個々のBHの質量・spin変化に対応する。
- Einstein–Maxwell の Majumdar–Papapetrou 解は、extremal same-sign charge において gravitational attraction と electric repulsion が釣り合う静的multi-black-hole解であり、第1連続系実験の有力なexact anchorである。
- 本研究ではこれらを外部近似力として dynamics に注入せず、full field evolution のreadout / validationにのみ使用する。

参考：
1. SpECTRE General Relativity documentation, `gr::psi_4`.
2. Bernuzzi, Nagar & Zenginoğlu, Phys. Rev. D 86, 104038 (2012), horizon absorption.
3. Álvaro-Díaz & Morras, Phys. Rev. D 114, 044033 (2026), horizon absorption in eccentric/precessing BBH inspirals.
4. Lucietti, Annales Henri Poincaré 22, 2437–2450 (2021), Majumdar–Papapetrou multi-black-hole equilibrium.

---

# 47. 仕様上の最終固定事項

第1実験では、

\[
\boxed{
\text{continuous classical first}
}
\]

を守る。

量子化仮説を初期のparameter selectionへ持ち込まない。

一方で保存・図化は、

\[
\boxed{
\text{future narrow-valley detection ready}
}
\]

とする。

これにより、

1. 古典連続挙動を独立に確立できる
2. 数値離散化artifactを分離できる
3. 後の \(m/n\) 反復実験との差を明確にできる
4. 後から深い谷が見つかった場合にも再走行なしで初期runを再解析できる

ことを狙う。
