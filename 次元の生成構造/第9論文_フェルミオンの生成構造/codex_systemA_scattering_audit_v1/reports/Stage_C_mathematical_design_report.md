# Stage C 状態依存散乱API・Candidate 0〜3 数学設計報告

## 0. 本報告の位置

**[コード上の事実]** Stage A/Bにより、次を確定結果として引き継ぐ。

1. 現行散乱核は入力パリティに対して盲であり、型依存経路応答を生成しない。
2. 現行核は半周期移動と可換であり、二チャネル全体の偶奇セクターを保存する。
3. F×Bで各出力チャネルが混合するのは型変換ではなくチャネル輸送である。
4. `B_to_A_transfer` は経路交換量ではなく、初期Bスペクトルに対するA出力の余弦類似度である。
5. 現行設定ではA/Bの直交\(\eta\)モードにより経路干渉項が消える。
6. 現行のチャネル別正規化は基準実験では実質恒等だが、状態依存散乱ではrawチャネル強度差を消し得る。

**[モデル定義]** 本報告は数学設計だけを行う。Candidate 0〜3の数値実装、単体テスト、比較実験、既存System A／System Bへの組込みは行わない。

---

## 1. 共通状態空間とパリティ読出し

### 1.1 状態

**[モデル定義]** 二チャネル入力を、

\[
\Psi=
\begin{pmatrix}
a\\
b
\end{pmatrix},
\qquad
a,b\in\mathcal H_u\otimes\mathcal H_\eta
\]

とする。内積は離散実装を含めて正定値エルミート内積とし、

\[
\|a\|^2=\langle a,a\rangle
\]

とする。

### 1.2 搬送波を除いた相互作用表示

**[モデル定義]** 状態依存特徴量は、搬送波を含む表示ではなくカーネル座標で計算する。由来別搬送波除去写像を \(D\) とし、

\[
\widehat a=Da,
\qquad
\widehat b=Db
\]

とする。

\(D\) は、各\(\eta\)由来成分を射影し、その成分固有の \(q\) を除去して再結合する線形写像である。Stage Bで用いた、

```text
m=1 成分 → q=+1 を除去
m=2 成分 → q=-1 を除去
```

の一般化に相当する。\(D\) は対応する由来部分空間上でノルムを保存し、逆写像を持つことをAPIの前提条件とする。

カーネル座標の半周期移動を \(P\) とし、搬送波を含む元の表示へ持ち上げた演算子を、

\[
\widetilde P=D^{-1}PD
\]

とする。搬送波が奇数波数なら \(D\) と元表示の単純な半周期移動は符号因子のため一般には可換しない。以下の「半周期移動に対する性質」は、常にカーネル座標の \(P\)、またはそれと同値な元表示の \(\widetilde P\) に関する記述である。

**[未導出]** 搬送波・\(\eta\)由来メタデータが欠け、\(D\)を一意に構成できない入力については、Candidate 1〜3を黙って評価してはならない。APIは未定義として停止する。

### 1.3 半周期移動と射影

**[モデル定義]**

\[
(P\psi)(u,\eta)=\psi(u+\pi,\eta),
\qquad
P^2=I,
\qquad
P^\dagger=P.
\]

\[
\Pi_B=\frac{I+P}{2},
\qquad
\Pi_F=\frac{I-P}{2}.
\]

各カーネル状態について、

\[
\widehat a_B=\Pi_B\widehat a,
\quad
\widehat a_F=\Pi_F\widehat a,
\quad
\widehat b_B=\Pi_B\widehat b,
\quad
\widehat b_F=\Pi_F\widehat b
\]

とする。

単波の無次元半周期指標は、

\[
c_A=
\frac{\operatorname{Re}\langle\widehat a,P\widehat a\rangle}
{\|\widehat a\|^2},
\qquad
c_B=
\frac{\operatorname{Re}\langle\widehat b,P\widehat b\rangle}
{\|\widehat b\|^2}
\]

である。ここでは添字 \(A,B\) はチャネル名であり、\(p_B\) の \(B\) は偶数型を表すため、混同しない。

\[
-1\le c_A,c_B\le1.
\]

ゼロノルム入力では未定義とし、エラーにする。

---

## 2. 全候補に共通する有効角

### 2.1 基準角

**[モデル定義]** 既存の反射パラメータ \(R\in[0,1]\) との互換性を、

\[
\theta_0(R)=\arcsin\sqrt R
\]

で保持する。

\[
R=\sin^2\theta_0,
\qquad
T=\cos^2\theta_0.
\]

### 2.2 共通の範囲包絡

**[モデル定義]** Candidate 0〜3の違いを無次元構造関数 \(F_j\) だけへ集約する。

\[
\Delta\theta_j
=\rho(\theta_0)F_j(\widehat a,\widehat b),
\]

\[
\rho(\theta_0)
=\frac{2}{\pi}\theta_0
\left(\frac{\pi}{2}-\theta_0\right).
\]

各候補は、

\[
-1\le F_j\le1
\]

を満たすように定義する。

\[
\theta_{\mathrm{eff},j}
=\theta_0+\Delta\theta_j.
\]

\(\rho\) は、

\[
0\le\rho(\theta_0)
\le
\min\left(
\theta_0,\frac{\pi}{2}-\theta_0
\right)
\]

を満たすため、

\[
0\le\theta_{\mathrm{eff},j}\le\frac{\pi}{2}.
\]

したがって、既存の \(R\) 掃引の意味を保ったまま、

\[
R_{\mathrm{eff},j}
=\sin^2\theta_{\mathrm{eff},j},
\qquad
T_{\mathrm{eff},j}
=\cos^2\theta_{\mathrm{eff},j}
\]

を得る。

**[数学的帰結]** \(R=0,1\) では \(\rho=0\) となり、すべての候補が既存端点へ一致する。最大補正幅は \(\theta_0=\pi/4\) における \(\pi/8\) である。

**[未導出]** この包絡は、端点保存と範囲保証のための共通設計であり、力学から導出されたものではない。将来変更する場合も、Candidate間で共通に変更し、候補差を包絡へ混入させない。

### 2.3 共通係数

**[モデル定義]**

\[
t_{\mathrm{eff}}
=e^{i\theta_{\mathrm{eff}}}
\cos\theta_{\mathrm{eff}},
\]

\[
r_{\mathrm{eff}}
=-ie^{i\theta_{\mathrm{eff}}}
\sin\theta_{\mathrm{eff}}.
\]

全候補で、

\[
U(\theta_{\mathrm{eff}})
=
\begin{pmatrix}
r_{\mathrm{eff}}&t_{\mathrm{eff}}\\
t_{\mathrm{eff}}&r_{\mathrm{eff}}
\end{pmatrix}
\]

を用いる。Candidate間で係数式・更新式・正規化規則を変えてはならない。違いは \(\Delta\theta_j\) の定義だけである。

### 2.4 ユニタリ性

**[数学的帰結]** 任意の実数 \(\theta_{\mathrm{eff}}\) について、

\[
|r_{\mathrm{eff}}|^2+|t_{\mathrm{eff}}|^2=1,
\]

\[
r_{\mathrm{eff}}^*t_{\mathrm{eff}}
+t_{\mathrm{eff}}^*r_{\mathrm{eff}}=0.
\]

したがって、

\[
U^\dagger U=I.
\]

\(\theta_{\mathrm{eff}}\) が入力状態に依存しても、各入力対に対して評価された \(U\) はユニタリであるため、

\[
\|a_{\mathrm{raw}}\|^2+\|b_{\mathrm{raw}}\|^2
=
\|a\|^2+\|b\|^2
\]

を自動的に保存する。

Candidate 1〜3は入力から角度を決めるため写像全体としては非線形だが、入力ごとの二チャネル結合ノルムは保存する。

---

## 3. raw物理出力と反復状態の分離

### 3.1 経路配列

**[モデル定義]**

\[
A_{a\to a}=r_{\mathrm{eff}}a,
\qquad
A_{b\to a}=t_{\mathrm{eff}}b,
\]

\[
A_{b\to b}=r_{\mathrm{eff}}b,
\qquad
A_{a\to b}=t_{\mathrm{eff}}a.
\]

物理的な正本は、

\[
a_{\mathrm{physical,raw}}
=A_{a\to a}+A_{b\to a},
\]

\[
b_{\mathrm{physical,raw}}
=A_{b\to b}+A_{a\to b}
\]

である。

経路単独ノルムと干渉項を、

\[
N_{a\to a}=\|A_{a\to a}\|^2,
\quad
N_{b\to a}=\|A_{b\to a}\|^2,
\]

\[
I_a
=2\operatorname{Re}
\langle A_{a\to a},A_{b\to a}\rangle
\]

のように分離する。B出力も同様である。

\[
\|a_{\mathrm{physical,raw}}\|^2
=N_{a\to a}+N_{b\to a}+I_a.
\]

干渉項を反射・交換のどちらかへ配分してはならない。

### 3.2 反復用状態

**[モデル定義]** 反復計算が単位チャネルノルムを要求する場合だけ、

\[
a_{\mathrm{iteration}}
=\mathcal N_A(a_{\mathrm{physical,raw}}),
\qquad
b_{\mathrm{iteration}}
=\mathcal N_B(b_{\mathrm{physical,raw}})
\]

を別フィールドとして生成する。

```text
physical canonical state:
    a_physical_raw
    b_physical_raw

optional iteration state:
    a_iteration
    b_iteration
```

`a_iteration`,`b_iteration` はraw物理出力を上書きしない。

**[数学的帰結]** 状態依存散乱がrawチャネルノルム差を生じた場合、チャネル別正規化はその差を1対1へ置き換える。したがって候補判定は必ずraw経路ノルム・raw干渉項・rawチャネルノルムに基づける。

### 3.3 直交\(\eta\)モードと応答の位置

**[数学的帰結]** 入力が、

\[
\langle a,b\rangle=0
\]

を満たすとき、すべての候補で係数は複素スカラーなので、

\[
I_a
=2\operatorname{Re}
\left(
r_{\mathrm{eff}}^*t_{\mathrm{eff}}
\langle a,b\rangle
\right)
=0,
\]

\[
I_b=0.
\]

したがって、Candidate 1〜3が状態依存角を生成しても、直交\(\eta\)モードそのものから経路干渉が発生するわけではない。

非零の状態依存応答は、

\[
N_{a\to a}
=\sin^2\theta_{\mathrm{eff}}\|a\|^2,
\qquad
N_{b\to a}
=\cos^2\theta_{\mathrm{eff}}\|b\|^2
\]

のような経路単独ノルムの変化として現れる。以下でいう「直交\(\eta\)モード下で候補が非零」とは、干渉項ではなく \(\Delta\theta_j\neq0\) を意味する。

---

## 4. Candidate 0 — 現行基準

### 4.1 定義

**[モデル定義]**

\[
F_0=0,
\qquad
\Delta\theta_0=0,
\qquad
\theta_{\mathrm{eff},0}=\theta_0.
\]

### 4.2 性質

**[数学的帰結]**

- 範囲: \(\Delta\theta_0=0\)
- 次元: 無次元角
- ユニタリ性: 成立
- A/B交換対称性: 成立
- 同時半周期移動との可換性: 成立
- 状態依存性: なし
- 直交\(\eta\)モード下の補正: 常にゼロ

Candidate 0はStage Bの再現基準であり、新しい型依存応答を生成しない。

---

## 5. Candidate 1 — 自己半周期相関主導

### 5.1 交差相関

**[モデル定義]**

\[
c_{AB}^{\mathrm{sym}}
=
\frac{
\operatorname{Re}
\left(
\langle\widehat a,P\widehat b\rangle
+\langle\widehat b,P\widehat a\rangle
\right)
}
{2\|\widehat a\|\|\widehat b\|}.
\]

\[
-1\le c_{AB}^{\mathrm{sym}}\le1.
\]

### 5.2 有効角補正

**[モデル定義]**

\[
F_1
=
\alpha_s\frac{c_A+c_B}{2}
+\alpha_x c_{AB}^{\mathrm{sym}},
\]

\[
|\alpha_s|+|\alpha_x|\le1.
\]

\[
\Delta\theta_1=\rho(\theta_0)F_1.
\]

\(\alpha_s,\alpha_x\) は無次元であり、比較実験前に固定する。データを見て符号や大きさを変更してはならない。

### 5.3 範囲と次元

**[数学的帰結]**

\[
|F_1|\le1,
\qquad
|\Delta\theta_1|\le\rho(\theta_0).
\]

すべて無次元である。

### 5.4 対称性

**[数学的帰結]**

- \(A\leftrightarrow B\): \(F_1\) は不変
- 二波への共通全体位相: 不変
- 独立なA/B位相: 自己項は不変、交差項は相対位相に依存
- \((\widehat a,\widehat b)\mapsto(P\widehat a,P\widehat b)\): \(F_1\) は不変

したがってCandidate 1の非線形散乱写像は、カーネル座標の同時半周期移動と可換する。元表示では、同じ命題を \(\widetilde P=D^{-1}PD\) に対して述べる。

### 5.5 純セクターでの構造値

**[数学的帰結]** 自己項について、

| 入力 | \((c_A+c_B)/2\) |
|---|---:|
| F×F | \(-1\) |
| B×B | \(+1\) |
| F×B | \(0\) |
| 50:50 MIX×MIX | \(0\) |

したがって \(\alpha_s\neq0\) なら、F×FとB×Bは反対向きの角度補正を持つ。F×Bの自己項は相殺する。

### 5.6 全周期内積と\(\eta\)直交

**[数学的帰結]** 初期状態が、

\[
\widehat a(u,\eta)=f_A(u)e^{im_A\eta},
\qquad
\widehat b(u,\eta)=f_B(u)e^{im_B\eta}
\]

で \(m_A\neq m_B\) なら、

\[
\langle\widehat a,P\widehat b\rangle
\propto
\int_{-\pi}^{\pi}
e^{i(m_B-m_A)\eta}\,d\eta
=0.
\]

よって、

\[
c_{AB}^{\mathrm{sym}}=0.
\]

一方、自己相関は、

\[
\langle\widehat a,P\widehat a\rangle
\propto
\int_{-\pi}^{\pi}|e^{im_A\eta}|^2\,d\eta
\]

を含むためゼロにならない。

**[数学的帰結]** 現行の初期直交\(\eta\)モード下でCandidate 1の角度補正が非零になる必要十分条件は、

```text
0 < R < 1
and
alpha_s != 0
and
c_A + c_B != 0
```

である。交差項だけに依存する \(\alpha_s=0\) の設計は、初回衝突ではCandidate 0へ退化するため棄却する。

チャネル輸送後は各出力チャネルが複数の\(\eta\)由来成分を含み得るため、交差相関が一般状態で常にゼロとは限らない。しかし初回応答を生成できない交差項だけの設計を、直交性を破る機構として採用してはならない。

**[実装候補]** 直交\(\eta\)設定に対する最小基準は、

\[
\alpha_s=1,\qquad\alpha_x=0
\]

である。ただし採用値の決定は数値比較前の承認事項とする。

---

## 6. Candidate 2 — 相互作用作用素付き双線形パリティ相関

### 6.1 恒等作用素の退化

**[数学的帰結]** 通常の全周期内積で、

\[
C_{XY}
=
\langle\widehat a_X,\widehat b_Y\rangle
\]

とすると、二つの独立なゼロ化が起こる。

1. \(X\neq Y\) では直交射影により \(C_{BF}=C_{FB}=0\)
2. \(m_A\neq m_B\) では\(\eta\)直交により同型成分 \(C_{BB},C_{FF}\) もゼロ

したがって \(K_{\mathrm{int}}=I\) は現行初期設定で全行列ゼロとなる退化基準であり、初回の非自明Candidate 2としては使用できない。

チャネル輸送後は各チャネルが複数の\(\eta\)由来成分を含むため、同型ブロック \(C_{BB},C_{FF}\) が一般状態で常にゼロとは限らない。一方、\(K=I\) での異型ブロック \(C_{BF},C_{FB}\) はパリティ直交により常にゼロである。ここでは、輸送後に偶発的に活性化する可能性と、初期相互作用を生成できることを区別する。

### 6.2 相互作用作用素

**[モデル定義]** 有界自己共役作用素 \(K_{\mathrm{int}}\) を事前固定し、

\[
K_+
=\frac{K_{\mathrm{int}}+PK_{\mathrm{int}}P}{2},
\]

\[
K_-
=\frac{K_{\mathrm{int}}-PK_{\mathrm{int}}P}{2}
\]

と分解する。

\[
PK_+P=K_+,
\qquad
PK_-P=-K_-.
\]

\(K_+\) は同じパリティセクターを結び、\(K_-\) は異なるパリティセクターを結ぶ。

### 6.3 正規化相関

**[モデル定義]** \(K_+\neq0\) のとき、

\[
d_+
=
\frac{
\operatorname{Re}
\left(
\langle\widehat a_B,K_+\widehat b_B\rangle
-
\langle\widehat a_F,K_+\widehat b_F\rangle
\right)
}
{
\|\widehat a\|\|\widehat b\|\|K_+\|_{\mathrm{op}}
}.
\]

\(K_-=0\) でなければ、

\[
o_-
=
\frac{
\operatorname{Re}
\left(
\langle\widehat a_B,K_-\widehat b_F\rangle
+
\langle\widehat a_F,K_-\widehat b_B\rangle
\right)
}
{
\|\widehat a\|\|\widehat b\|\|K_-\|_{\mathrm{op}}
}.
\]

対応する作用素成分がゼロなら、その特徴量を0と定義し診断フラグを返す。

Cauchy–Schwarz不等式により、

\[
|d_+|\le1,
\qquad
|o_-|\le1.
\]

### 6.4 有効角補正

**[モデル定義]**

\[
F_2
=\beta_d d_+
+\beta_o o_-^2,
\]

\[
|\beta_d|+|\beta_o|\le1.
\]

\[
\Delta\theta_2=\rho(\theta_0)F_2.
\]

\(o_-^2\) を用いる理由は、異型交差相関の大きさを保持しながら、同時半周期移動で符号反転する \(o_-\) を不変量へ変えるためである。

### 6.5 範囲・次元・対称性

**[数学的帰結]**

- \(F_2\in[-1,1]\)
- \(\Delta\theta_2\) は無次元角
- \(K\) の物理次元は \(\|K\|_{\mathrm{op}}\) で除かれる
- \(A\leftrightarrow B\): 自己共役 \(K\) と実部対称化により不変
- 共通全体位相: 不変
- 独立なA/B位相: 一般に相対位相依存
- 同時半周期移動: \(d_+\) は不変、\(o_-\) は符号反転、\(o_-^2\) は不変

したがってCandidate 2も同時半周期移動と可換する。

### 6.6 初期\(\eta\)直交下で非零になる条件

**[数学的帰結]** \(u\) だけの局所窓、

\[
K_{\mathrm{int}}=W(u)\otimes I_\eta
\]

は\(\eta\)直交を破らないため、現行の \(m_A\neq m_B\) では全相関がゼロのままである。局所窓だけでは不十分である。

同型相関 \(d_+\) を非零にするには、少なくとも、

\[
K_{\eta,\mathrm{swap}}
=
|m_A\rangle\langle m_B|
+|m_B\rangle\langle m_A|
\]

のような自己共役\(\eta\)モード橋渡しが必要である。

異型相関 \(o_-\) を非零にするには、さらに \(u\) パリティを反転する作用が必要である。例えば衝突窓 \(W\) から、

\[
W_-(u)
=\frac{W(u)-W(u+\pi)}{2}
\]

を作り、

\[
K_-
=W_-(u)\otimes K_{\eta,\mathrm{swap}}
\]

とすれば、

\[
PK_-P=-K_-.
\]

**[数学的帰結]** Candidate 2の非零条件は次のように分かれる。

```text
same-parity response d_+:
    0 < R < 1 and beta_d != 0
    eta-mode bridge is required
    and the bridged spatial components must overlap

cross-parity response o_-:
    0 < R < 1 and beta_o != 0
    eta-mode bridge is required
    and a parity-odd interaction operator is required
    and the contact-weighted components must overlap
```

**[未導出]** \(K_{\mathrm{int}}\)、衝突窓 \(W\)、\(\beta_d,\beta_o\) の物理的起源と値は未導出である。結果を見て窓や結合を調整してはならず、数値実装前に全候補を固定する必要がある。

---

## 7. Candidate 3 — 積波の\(\mathbb Z_2\)パリティ構成

### 7.1 積波

**[モデル定義]** 同じ点ごとの格子・測度を持つカーネル状態について、

\[
g_{AB}(u,\eta)
=\widehat a(u,\eta)\widehat b(u,\eta)
\]

を作る。複素共役を入れない積である。

積波を、

\[
g_B=\Pi_B g_{AB},
\qquad
g_F=\Pi_F g_{AB}
\]

へ分解する。

### 7.2 積パリティ指標

**[モデル定義]**

\[
\chi_\times
=
\frac{\|g_B\|^2-\|g_F\|^2}
{\|g_B\|^2+\|g_F\|^2}.
\]

\[
-1\le\chi_\times\le1.
\]

\(g_{AB}=0\) がほとんど至る所で成立する場合は接触なしとみなし、

\[
\chi_\times=0,
\qquad
\Delta\theta_3=0
\]

とする。同時に `no_pointwise_overlap=true` を診断へ返す。

### 7.3 有効角補正

**[モデル定義]**

\[
F_3=\gamma\chi_\times,
\qquad
|\gamma|\le1,
\]

\[
\Delta\theta_3
=\rho(\theta_0)F_3.
\]

### 7.4 \(\mathbb Z_2\)合成則

**[数学的帰結]** 純パリティ状態では、

| 入力 | 積波 | \(\chi_\times\) |
|---|---|---:|
| F×F | 偶数型 | \(+1\) |
| B×B | 偶数型 | \(+1\) |
| B×F / F×B | 奇数型 | \(-1\) |

すなわち、

\[
F\times F\to B,
\qquad
B\times B\to B,
\qquad
B\times F\to F
\]

を一つの式で読出す。

Candidate 3はF×FとB×Bを区別せず、同型対と異型対を区別する候補である。

### 7.5 範囲・次元・対称性

**[数学的帰結]**

- \(\chi_\times,F_3\) は無次元
- 積波の場の次元は比で相殺される
- A/B交換: 点ごとの積が可換なので不変
- A/Bの独立な全体位相: 積波全体の位相となり、射影ノルム比では消える
- 同時半周期移動: \(g_{AB}\mapsto Pg_{AB}\) なので\(\chi_\times\)は不変

したがってCandidate 3も同時半周期移動と可換する。

### 7.6 全周期積分と\(\eta\)直交

**[数学的帰結]** 素朴な線形量、

\[
\int g_{AB}(u,\eta)\,du\,d\eta
\]

は、\(m_A+m_B\neq0\) のとき\(\eta\)積分でゼロになり得る。またB×Fの積波は奇数型なので、\(u\)全周期積分でもゼロになる。

したがって、積波の全周期積分をそのまま \(\Delta\theta_3\) に使う案は棄却する。

一方、

\[
\|g_B\|^2+\|g_F\|^2
=\|g_{AB}\|^2
\]

では\(\eta\)位相の絶対値二乗が残る。A/Bが直交\(\eta\)モードであっても、点ごとの積がゼロでない限り \(\chi_\times\) は定義され、純セクターでは \(\pm1\) となる。

**[数学的帰結]** Candidate 3が現行の直交\(\eta\)モード下で非零になる条件は、

```text
0 < R < 1
and
gamma != 0
and
a_hat * b_hat is not zero almost everywhere
and
the product parity contrast chi_cross is not zero
```

である。通常のA–B内積は使用しないため、\(\langle a,b\rangle=0\) だけではCandidate 3は消えない。

---

## 8. 全周期内積・\(\eta\)直交に対する事前判定

**[数学的帰結]**

| 候補 | 素朴形の問題 | 採用形 | 現行初期直交\(\eta\)下の非零条件 |
|---|---|---|---|
| Candidate 0 | 状態依存性なし | \(F_0=0\) | 常に補正ゼロ |
| Candidate 1 | A–B交差半周期相関はゼロ | 自己相関主項＋交差診断 | \(\alpha_s(c_A+c_B)\neq0\) |
| Candidate 2 | \(K=I\)では全ブロックがゼロ | \(\eta\)橋渡しを持つ固定 \(K_\pm\) | bridge、空間重なり、異型にはparity-odd作用 |
| Candidate 3 | 積波の線形全周期積分はゼロになり得る | 積波の偶奇射影ノルム比 | 点ごとの積が非零かつ\(\chi_\times\neq0\) |

**[棄却候補]**

1. Candidate 1を交差相関だけで構成する案
2. Candidate 2で \(K_{\mathrm{int}}=I\) を非自明候補として扱う案
3. Candidate 2で \(W(u)\otimes I_\eta\) だけを用いる案
4. Candidate 3で \(\int\widehat a\widehat b\) をそのまま角度補正にする案
5. 偶数・奇数ラベルによる条件分岐

これらは現行初期設定で自明にゼロ化するか、必要な異型交差応答を恒等的に失うか、仮説を直接埋め込むため採用しない。

---

## 9. 候補比較

**[数学的帰結]**

| 性質 | C0 | C1 | C2 | C3 |
|---|---|---|---|---|
| 共通 \(U(\theta_{\mathrm{eff}})\) | 使用 | 使用 | 使用 | 使用 |
| \(\Delta\theta\)の状態依存 | なし | 自己・交差半周期相関 | \(K\)付き双線形相関 | 積波パリティ |
| ユニタリ性 | 自動 | 自動 | 自動 | 自動 |
| A/B交換対称性 | あり | あり | あり | あり |
| 同時半周期移動と可換 | あり | あり | 採用形ではあり | あり |
| 全周期・\(\eta\)直交への耐性 | 該当なし | 自己項が残る | bridge必須 | 積射影ノルムが残る |
| F×FとB×Bの区別 | なし | 反対符号 | \(d_+\)で可能 | しない |
| 同型対と異型対の区別 | なし | F×Bでは自己項相殺 | \(o_-^2\)で可能 | 直接可能 |
| 新しい未導出構造 | なし | 結合係数 | \(K\)、窓、結合係数 | 積相互作用の物理解釈 |

**[未導出]** どの候補が物理的に採用されるべきかは、この数学設計だけでは決まらない。本報告は「非零になり得ること」と「保存・対称性を壊さないこと」を確認した段階であり、実在の散乱法則を導出したものではない。

---

## 10. 共通APIの数学仕様

### 10.1 入力

**[モデル定義]** 実装時の概念的入力を次とする。

```text
ScatteringRequest
    a: complex field array
    b: complex field array
    reflection_parameter: R in [0,1]
    coordinate_u
    half_shift_operator: P
    carrier_lineage_spec: specification for D
    model: Candidate 0 | 1 | 2 | 3
    fixed_model_parameters
    iteration_normalization: none | per_channel_l2 | pair_l2
```

Candidate 2では、事前固定した \(K_{\mathrm{int}}\) とその識別子・ハッシュを入力仕様へ追加する。

Candidate 3では、A/Bが同じ点ごとの格子・測度を持つことを必須とする。

### 10.2 出力

**[モデル定義]**

```text
StateDependentScatteringResult
    model
    theta_0
    delta_theta
    theta_eff
    t_eff
    r_eff
    T_eff
    R_eff

    a_physical_raw
    b_physical_raw

    path_a_to_a_amplitude
    path_b_to_a_amplitude
    path_b_to_b_amplitude
    path_a_to_b_amplitude

    path_a_to_a_norm
    path_b_to_a_norm
    path_b_to_b_norm
    path_a_to_b_norm

    interference_in_a
    interference_in_b
    output_a_raw_norm
    output_b_raw_norm

    a_iteration
    b_iteration
    iteration_normalization_policy
    iteration_scale_a
    iteration_scale_b

    parity_features_input
    parity_features_raw_output
    candidate_feature_values

    norm_before_pair
    norm_after_raw_pair
    norm_after_iteration_pair
    unitarity_residual
    half_shift_commutator_response

    diagnostics
    legacy_metrics
```

`a_iteration`,`b_iteration` は、正規化方針が `none` の場合はraw出力と同じ参照または明示的な `None` とする。

### 10.3 意味上の固定

**[モデル定義]**

```text
physical truth:
    a_physical_raw
    b_physical_raw
    raw path amplitudes
    raw path norms
    raw interference

iteration aid:
    a_iteration
    b_iteration

legacy diagnostic only:
    B_to_A_transfer
```

`B_to_A_transfer` を残す場合、必ず、

```text
spectral cosine similarity of A output to the initial B spectrum;
NOT a path-exchange norm
```

と記録する。

経路交換量の互換名は、

```text
exchange_b_to_a = path_b_to_a_norm
exchange_a_to_b = path_a_to_b_norm
```

に限る。干渉項を含めず、`B_to_A_transfer` とは別フィールドにする。

### 10.4 API不変条件

**[数学的帰結]** 全候補で次を要求する。

\[
\texttt{a\_physical\_raw}
=
\texttt{path\_a\_to\_a}
+
\texttt{path\_b\_to\_a},
\]

\[
\texttt{b\_physical\_raw}
=
\texttt{path\_b\_to\_b}
+
\texttt{path\_a\_to\_b},
\]

\[
\texttt{norm\_after\_raw\_pair}
=
\texttt{norm\_before\_pair}
\]

を数値許容差内で満たす。

候補間比較ではraw物理出力を使い、反復用正規化状態で判定してはならない。

---

## 11. 実装前に固定すべき未導出項目

**[未導出]**

1. Candidate 1の \(\alpha_s,\alpha_x\)
2. Candidate 2の \(K_{\mathrm{int}}\)、接触窓、\(\beta_d,\beta_o\)
3. Candidate 3の \(\gamma\)
4. 共通包絡 \(\rho\) を今後も維持するか
5. 相対チャネル位相を物理入力として扱うか
6. 反復状態の正規化方針
7. Candidate 2の\(\eta\)モード橋渡しを、どの力学から与えるか
8. Candidate 3の積相互作用を、どの物理量へ対応させるか

これらは数値結果を見て後付けで選ばず、実装前に固定する。

---

## 12. Stage C停止

**[コード上の事実]** 本段階で作成したのは本Markdown報告書だけである。Candidate 0〜3の数値モジュール、テスト、比較データ、図、既存System A／System Bへの組込みは作成していない。

Stage C数学設計完了。Candidate 0〜3の数値実装および既存System A／System B本体への組込みは未実施であり、人間の承認を待つ。
