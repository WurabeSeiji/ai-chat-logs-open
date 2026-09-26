# 第五思考実験：隠れた状態はなかったか
## ― 第四思考実験の二状態軌道生成を自己監査し、五つの永続状態・同期相互作用・作用列・近似次数を検証する ―

**著者:** 木原範昭 (Noriaki Kihara)  
**ORCID:** 0009-0004-6753-4020  
**Version DOI:** （取得後記入）  
**Concept DOI:** （取得後記入）  
**版:** v1.0  
**日付:** 2026-09-24  

---

## 要旨

第四思考実験 [3] では、二成分内部状態

$$
X_n=\begin{pmatrix}a_n\\b_n\end{pmatrix}
$$

に状態依存作用 $S_n$ を逐次作用させ、二次読み出し

$$
\Phi(X_n)=
\begin{pmatrix}
a_n^2-b_n^2\\
2a_nb_n
\end{pmatrix}
$$

から、近点移動と 2.5PN 放射反作用による縮退を含む重力軌道族を再構成できることを示した。しかし、その実装を改めて監査すると、$X_n$ の外側で処理位相 $\chi$、半直弦 $p$、離心率 $e$、方位角 $\phi$、時計 $t$ が保持・更新され、それらから $r$、$S_n$、$X_n$ が構成されていた。したがって「二状態と相互作用だけで軌道を生成した」という読み方は厳密ではなかった。

本稿では、この問題を第四思考実験の自己監査として扱う。監査規則は厳しい。次状態は現在の全状態だけを一度に読み込み、積・積和により一度に書き込む。同一処理中に先に更新された次状態を再利用しない。外部の処理番号や暗黙のレジスタを物理情報として参照しない。処理段階が必要なら、それ自体を位相状態として内部化する。固定定数は導出根拠が明確な場合に限り使用し、結果に合わせた自由な調整係数は認めない。読み出しに $\log$ を用いることは許すが、その値を状態更新へ戻さない。

この規則で Paper 4 の数式とプログラムを順に展開したところ、マクロな処理stepを越えて保持される状態として、現在の表現では少なくとも

$$
\boxed{(\chi,p,e,\phi,t)}
$$

の五つが必要であることが分かった。これは五状態が数学的最小であることの証明ではなく、現在の表現で隠れた永続状態を除去した結果である。さらに加法的な五状態を、乗法更新に適した

$$
\boxed{
Z=(U,P,E,H,Q)
=(e^{i\chi},p,e,e^{i\phi},e^{t/\tau_0})
}
$$

へ再符号化した。$\log$ は時計の観測時にのみ用い、状態遷移へはフィードバックしない。

Paper 4 の RK4 処理を隠れた逐次レジスタなしに実行するため、評価点、$k_1,\ldots,k_4$、中点、合成増分、および処理位相を明示的な内部状態として展開した 11 相の状態機械を構成した。500 個のランダム状態に対する一マクロstep照合では五状態すべてで差が 0 となり、12 条件の長時間照合でも最大差は $H$ で $3.38\times10^{-16}$、$Q$ で $4.00\times10^{-15}$ であった。保存済み Paper 4 C1 データとの照合でも $p,e,t,x,y$ は浮動小数点丸めの範囲で一致した。

また、放射反作用演算と $N=12$ 角度級数を $U^{\pm1}$、$P$、$E$ の明示的積和へ展開し、10 万ランダム点で元式と照合した。さらに角度級数を $N=12,16,20,24$ と上げると、最も厳しい $p=24.3,e=0.55$ の 3 軌道後の方位角誤差は $1.56\times10^{-7}$、$9.03\times10^{-10}$、$5.54\times10^{-12}$、$3.55\times10^{-14}$ と減少した。五つの永続状態を固定した Taylor 漸化でも 4、6、8 次への高次化が可能で、9 条件の 16→32 step/orbit の比較では $P,E,\phi$ の実測次数がそれぞれ設計次数に近い値を示した。

本稿は一般相対論の厳密解を状態相互作用から導出したものではない。また五状態の最小性や、複数の内部作用 $S_n$ が原理的に不可欠であることを証明していない。示したのは、背景時空や天体軌道を基本状態として保持せず、厳しい同期更新規則の下でも、既知の重力モデルと比較可能な軌道読み出しを有限状態の相互作用系として再構成できるという限定的な表現可能性である。

**キーワード:** 隠れた状態、状態相互作用、同期更新、乗法状態、状態依存作用、ポストニュートン近似、離散力学、数値監査、思考実験

---

# 0. 本稿の位置づけ

本稿は研究プロジェクト設計 [0] のもとでの第五思考実験である。

第二思考実験 [1] では、意味を定めない二つの値からなる状態に固定作用 $S$ を作用させ、二次位置読み出しと面積時計を介して逆二乗型運動を構成した。

第三思考実験 [2] では向きを逆転し、既知の Newtonian / PN 二体軌道から質量比・相対距離・相対時間を読み出せるかを調べた。

第四思考実験 [3] では、第三思考実験で外部計算していた軌道生成を再び内部相互作用へ戻すため、

$$
X_{n+1}=S_nX_n,
\qquad
S_n=\mathcal S(\mathcal R_n),
\qquad
\det S_n=1
$$

という状態依存作用を構成した。

しかし第五思考実験の問いは、第四思考実験の結果をさらに拡張することではない。

$$
\boxed{
\text{第四思考実験は、本当に明示した状態だけで閉じていたのか。}
}
$$

である。

本稿では札を次の意味で使う。

| 札 | 意味 |
|---|---|
| 【原則】 | 本稿で自己監査のために維持する計算原則 |
| 【設計仮定】 | 現在の表現・実装のために置く仮定。将来外してよい |
| 【構成】 | 既知式または明示状態から作用を構成する手順 |
| 【数値】 | プログラムによる数値照合で確認した結果 |
| 【照合】 | 既知物理・既存数学との比較。基礎原理からの導出ではない |
| 【未証明】 | 現在の実験では決められない事項 |

本稿の最も重要な方針は、**第四思考実験の成功を前提にしない**ことである。隠れた状態が見つかれば、それを失敗として隠さず、明示状態へ移す。そのため、本稿では途中で却下した実装や、監査中に発見した違反も再現可能な実験履歴として保存する。

---

# 1. 動機：二状態だけで本当に軌道を生成していたのか

第四思考実験の表面上の状態は

$$
X_n=(a_n,b_n)^{\mathsf T}
$$

であり、軌道は

$$
Y_n=\Phi(X_n)
$$

から読み出されていた。

しかし実際のプログラムでは、各 step で

$$
\chi_n,
\quad p_n,
\quad e_n,
\quad \phi_n,
\quad t_n
$$

を別に保持し、そこから

$$
r_n=\frac{p_n}{1+e_n\cos\chi_n},
\qquad
\theta_n=\frac{\phi_n}{2}
$$

を求め、さらに $S_n$ と $X_n$ を構成していた [3]。

したがって実際の情報流は

$$
(\chi,p,e,\phi,t)
\longrightarrow
(r,\theta)
\longrightarrow
S_n
\longrightarrow
X
\longrightarrow
Y
$$

であり、$X$ は上流状態から導かれた表現であった。

![図1　Paper 4 の見かけの二状態と、監査で明示された永続状態](figures/figure00_hidden_state_audit_structure.svg)

**図1.** 第四思考実験では二状態 $X=(a,b)$ が基本状態に見えていたが、実装を上流へ追跡すると $\chi,p,e,\phi,t$ が保持・更新されていた。本稿はこれらを隠れた永続状態として明示化する。

本稿の第一の問いは、

$$
\boxed{
\text{次状態を決めるために必要な情報を、すべて状態として数え直すと何状態必要か。}
}
$$

である。

第二の問いは、

$$
\boxed{
\text{それらの状態を一つの固定作用 }S\text{ だけで同期更新できるか。}
}
$$

である。

第三の問いは、Paper 4 が複数の有限次数近似から構成されていることに対して、

$$
\boxed{
\text{状態数や作用の種類を増やさず、近似次数だけを上げられるか。}
}
$$

である。

---

# 2. 本稿の目的と主張範囲

本稿の目的は一般相対論の新しい厳密解を導出することではない。

また、ここで構成する状態遷移則が自然界の基本法則であると主張するものでもない。

目的は、通常の物理計算では暗黙に利用しやすい

- 外部時計、
- 過去状態、
- 演算step番号、
- 更新途中の値、
- 暗黙のレジスタ、
- 結果に合わせた自由な補正係数、

を極力排除し、

$$
\boxed{
\text{現在の有限状態}
\rightarrow
\text{相互作用}
\rightarrow
\text{次の有限状態}
\rightarrow
\text{読み出し}
}
$$

という厳しい計算条件でも、既知の物理的軌道と比較可能な読み出しを構成できるかを検査することである。

したがって本稿が検証するのは

$$
\boxed{
\text{制約された状態相互作用系が、既知物理を表現できるか}
}
$$

という**表現可能性**である。

これは

$$
\boxed{
\text{制約された状態相互作用系が、自然界の基礎法則である}
}
$$

ことの証明ではない。

本稿で確認したことと未確認のことを先に分ける。

### 本稿で確認したこと

1. Paper 4 の二状態 $X$ だけでは実装は閉じていなかった。
2. 現在の表現では、マクロstepを越えて保持される状態として $\chi,p,e,\phi,t$ の五状態を明示すると Paper 4 の軌道計算を再構成できた。
3. Paper 4 の RK4 逐次計算を、処理段階を含む明示内部状態へ展開すると、同一 microstep 中に更新後値を再利用しない同期状態機械として再構成できた。
4. 角度級数と時間積分の次数は、五つの永続状態を増やさず高次化できた。

### 本稿では証明していないこと

1. 五状態が真の最小状態数であること。
2. 一つの固定作用 $S$ では原理的に不可能であり、作用列 $S_n$ が必ず必要であること。
3. 近似次数を無限に上げた極限が一般相対論の厳密解へ一致すること。
4. 本稿の相互作用則が自然界の基本法則であること。

---

# 3. 自己監査規則

## 3.1 現在状態だけから次状態を作る

状態集合を

$$
\Psi_n=(\psi_{1,n},\ldots,\psi_{N,n})
$$

とする。

一回の相互作用では、現在状態を一度に読み、次状態を一度に書く。

$$
\boxed{
\psi_{i,n+1}=F_i(\Psi_n)
}
$$

でなければならない。

同じ一回の相互作用中に先に計算した $\psi_{j,n+1}$ を別成分の計算に使うことは禁止する。

---

## 3.2 増分代入と積和を区別する

プログラム上の

```text
psi += delta
```

は、旧値と新値を暗黙に保持する逐次更新として読めるため、基本相互作用としては採用しない。

一方、現在状態だけから

$$
\psi_i'
=
\sum_{\alpha}
c_{i\alpha}
\prod_j \psi_j^{m_{i\alpha j}}
$$

を一度に計算して書き込む**積和**は許容する。

例えば

$$
P'=P+D_p
$$

は、逐次代入ではなく

$$
P'=1\cdot P+1\cdot D_p
$$

という二つの現在状態の積和として一度に計算する限り、本稿の規則に反しない。

---

## 3.3 演算回数を外部から参照しない

外部の loop counter によって作用を選んではならない。

処理段階の区別が必要なら、one-hot 位相状態

$$
q_n=(q_{0,n},\ldots,q_{m-1,n})
$$

を系内部に置き、巡回置換

$$
\boxed{
q_{n+1}=C_m q_n,
\qquad C_m^m=I
}
$$

で更新する。

作用の選択も数式上は

$$
\boxed{
\Psi_{n+1}
=
\sum_{j=0}^{m-1} q_{j,n}F_j(\Psi_n)
}
$$

と積和で書く。

実際の監査コードでは $q$ が厳密な one-hot であるため高速化として `argmax(q)` と分岐を使用するが、これは物理的な外部情報ではなく、上式の one-hot 選択を実装したものとして扱う。したがって、**数式上の基本定義は $q$ による積和選択**であり、`argmax` そのものを基本相互作用とはみなさない。

---

## 3.4 読み出しは状態へ書き戻さない

観測のための写像

$$
O=\mathcal O(\Psi)
$$

は許される。

しかし、特に

$$
\log
$$

を用いて加法座標へ戻し、旧来の加法計算を行って再び指数化する

$$
\Psi
\rightarrow
\log\Psi
\rightarrow
\text{additive update}
\rightarrow
\exp
\rightarrow
\Psi'
$$

という手順は基礎モデルとして認めない。

読み出しは

$$
\boxed{
\text{state}\rightarrow\text{observation}
}
$$

の一方向とする。

---

## 3.5 定数

固定定数との積和は、由来が明確なら許容する。

例えば RK4 の

$$
\frac16(1,2,2,1)
$$

や、級数

$$
c_{m+1}=\frac{2m+1}{m+1}c_m
$$

の係数は導出根拠が明確である。

一方、既知軌道に近づけるためだけの自由な fitting coefficient は用いない。

![図2　本稿で課す相互作用規則](figures/figure00b_interaction_constraints.svg)

**図2.** 現在状態だけを読み、積・積和により次状態を同期的に書く。外部 step 番号や更新途中の新値を物理情報として用いない。

---

# 4. 既存理論での位置づけ

本稿の計算規則は、組み合わせると厳しく見えるが、個々の考え方は理論物理・数理物理で特異なものではない。

まず、少数の変数へ射影した有効記述では、消去した自由度の影響が memory kernel として現れ得ることは Zwanzig の射影形式に代表される [4]。これは「現在状態だけで閉じないなら、過去をそのまま参照する前に欠落自由度を疑う」という本稿の監査方針の背景になる。ただし Zwanzig の結果から、本稿の五状態が一意に導かれるわけではない。

次に、時間発展を複数作用素の積として構成することは Trotter 積公式 [5] および Suzuki の系統的高次近似 [6] で標準的である。本稿の内部作用列はこれらと同一ではないが、**一つの連続発展を作用の積へ展開する**という表現自体は一般的である。

離散時刻の各状態を同一時刻の旧状態から同期的に更新する形式は、セルオートマトンを物理モデルとして扱う研究にも明確な先例がある [7]。本稿はセルオートマトンではないが、同期更新という計算原則の先例として参照する。

周期依存を拡張された状態空間へ持ち上げる考えは Floquet-Sambe 形式 [8] と構造的に類似する。本稿の one-hot 処理位相は量子 Floquet 状態ではないが、外部の周期ラベルを内部状態化するという点で対応する。

また、高次の時間積分を新しい物理自由度の追加と区別すべきことは、高次 symplectic integrator の構成 [9] からも明らかである。本稿の Taylor 漸化は Yoshida 法そのものではないが、**高次化が新しい物理状態の導入を意味しない**という数値解析上の背景を共有する。

したがって本稿固有の点は、これらの個別概念を新規発明したことではなく、**厳しい状態監査規則を同時に課し、前報の実装を自己監査したこと**にある。

---

# 5. Paper 4 の状態監査

## 5.1 表面状態 $X=(a,b)$ は基本状態ではなかった

Paper 4 では

$$
X_{n+1}=S_nX_n
$$

を主要な状態遷移として表示した [3]。

しかし実装では先に

$$
(\chi_n,p_n,e_n,\phi_n,t_n)
$$

を更新し、

$$
r_n=\frac{p_n}{1+e_n\cos\chi_n},
\qquad
\theta_n=\frac{\phi_n}{2}
$$

から

$$
X_n
=
\sqrt{r_n}
\begin{pmatrix}
\cos\theta_n\\
\sin\theta_n
\end{pmatrix}
$$

を構成していた。

したがって $X_n$ は上流状態の圧縮された**派生表現**であり、それだけでは次状態を決定できない。

---

## 5.2 五つの永続状態

監査の結果、現在の表現でマクロstepをまたいで必要な状態は

$$
\boxed{
S_n^{\rm macro}
=(\chi_n,p_n,e_n,\phi_n,t_n)
}
$$

であった。

ここで

- $\chi$：処理位相、
- $p$：半直弦、
- $e$：離心率、
- $\phi$：軌道方位角、
- $t$：時計読み出しの累積座標、

である。

表1に監査結果を示す。

### 表1　Paper 4 の状態監査

| 量 | Paper 4 での役割 | 永続状態か | 本稿での扱い |
|---|---|---:|---|
| $X=(a,b)$ | 内部2状態として表示 | No | $r,\phi$ から作る派生表現 |
| $\chi$ | 処理位相 | Yes | $U=e^{i\chi}$ へ再符号化 |
| $p$ | 軌道要素 | Yes | $P=p$ |
| $e$ | 軌道要素 | Yes | $E=e$ |
| $\phi$ | 方位角 | Yes | $H=e^{i\phi}$ へ再符号化 |
| $t$ | 時計 | Yes | $Q=e^{t/\tau_0}$ へ再符号化 |

五状態同期写像を Paper 4 の元計算と 12 条件で照合したところ、最悪値は

$$
|\Delta\chi|_{\max}=5.03\times10^{-12},
$$

$$
|\Delta p|_{\max}=4.19\times10^{-13},
$$

$$
|\Delta e|_{\max}=2.99\times10^{-14},
$$

$$
|\Delta\phi|_{\max}=1.95\times10^{-13},
$$

$$
|\Delta t|_{\max}=1.35\times10^{-8}
$$

であった。$t$ の絶対誤差が大きく見えるのは、長時間累積量に対する浮動小数点の加算順序差を含むためである。軌道読み出しでは

$$
|\Delta(x,y)|_{\max}=4.58\times10^{-10}
$$

であった。【数値】

![図3　五状態同期写像の再構成誤差](figures/figure02_five_state_reconstruction_errors.svg)

**図3.** Paper 4 の逐次処理と、五状態を明示した同期写像の差。これは五状態の最小性を証明する図ではなく、五状態で現在の実装を閉じられることの照合である。

---

## 5.3 「五状態」の意味

ここで「五状態」という言葉を厳密に限定する。

これは

$$
\boxed{
\text{マクロstepを越えて保持される永続状態が五つ}
}
$$

という意味である。

RK4 を隠れた逐次レジスタなしに実装するためには、$k_1,\ldots,k_4$、評価点、中点、合成増分、処理位相などを**明示内部状態**として置く必要がある。

したがって、全状態機械の実数成分数が五個だけという意味ではない。

また五が数学的最小であることも未証明である。

$$
\boxed{
\text{現在の状態表現では、五つの永続状態で閉じた}
}
$$

というのが本稿の正確な結論である。

---

# 6. 乗法状態への再符号化

加法座標

$$
(\chi,p,e,\phi,t)
$$

は監査には便利だが、本稿の「乗法を基本相互作用とする」という方針には必ずしも適していない。

そこで永続状態を

$$
\boxed{
Z=(U,P,E,H,Q)
}
$$

とし、

$$
U=e^{i\chi},
\qquad P=p,
\qquad E=e,
\qquad H=e^{i\phi},
\qquad Q=e^{t/\tau_0}
$$

とする。

このとき位相と時計は

$$
U'=Ue^{i\Delta\chi},
$$

$$
H'=He^{i\Delta\phi},
$$

$$
Q'=Qe^{\Delta t/\tau_0}
$$

という乗法更新になる。

一方、$P,E$ は現在の物理式において正の軌道要素として直接現れるので、本稿では無理に $\log P,\log E$ を状態として持たない。

【設計仮定】現在の検証プログラムは $U,H$ を複素数で簡潔に実装している。しかしこれは「基本状態が本質的に複素数でなければならない」という結論ではない。例えば

$$
U=(\cos\chi,\sin\chi)
$$

を実二成分とし、回転行列で更新することもできる。したがって、本稿の五状態という主張は**論理状態量の数**についてであり、最小実スカラー成分数を確定するものではない。

---

# 7. B方式：11相の同期内部状態機械

## 7.1 なぜ内部相が必要か

Paper 4 の放射反作用更新は RK4 であり、通常のコードでは

$$
k_1\rightarrow y_2\rightarrow k_2\rightarrow y_3\rightarrow k_3\rightarrow y_4\rightarrow k_4
$$

と逐次処理する。

これを一つの物理相互作用だとみなすと、処理途中の値が隠れたレジスタになる。

そこで各中間値を明示状態とし、11 回の microstep に分解した。

内部状態を模式的に

$$
\mathcal Z
=(Z,W,q)
$$

と書く。$Z=(U,P,E,H,Q)$ が永続状態、$W$ が $k_1,\ldots,k_4$、評価点、合成増分、中点などの work state、$q$ が 11 相 one-hot 状態である。

一 microstep は

$$
\boxed{
\mathcal Z_{m+1}
=
\sum_{j=0}^{10}q_{j,m}F_j(\mathcal Z_m)
}
$$

であり、同時に

$$
\boxed{
q_{m+1}=C_{11}q_m
}
$$

とする。

11 microstep 後に

$$
q_{m+11}=q_m
$$

となり、一つのマクロ状態更新が完了する。

---

## 7.2 実際の監査コード

以下は厳密監査に使用した中核部分である。`o` は更新前状態、`n` は次状態の書込先であり、各分岐の右辺は `o` または microstep 開始時の `U,P,E` だけを読む。

```python
def micro_fast(s, nu):
    j = int(np.argmax(s['q']))

    # old state -> new buffer
    o = {
        k: (v.copy() if isinstance(v, np.ndarray) else v)
        for k, v in s.items()
    }
    n = {
        k: (v.copy() if isinstance(v, np.ndarray) else v)
        for k, v in s.items()
    }

    U, P, E = o['U'], o['P'], o['E']

    if j == 0:
        n['k1'] = rr(U, P, E, nu)

    elif j == 1:
        n['Us'] = U * PH2
        n['Ps'] = P + 0.5 * H * o['k1'][0]
        n['Es'] = E + 0.5 * H * o['k1'][1]

    elif j == 2:
        n['k2'] = rr(o['Us'], o['Ps'], o['Es'], nu)

    elif j == 3:
        n['Us'] = U * PH2
        n['Ps'] = P + 0.5 * H * o['k2'][0]
        n['Es'] = E + 0.5 * H * o['k2'][1]

    elif j == 4:
        n['k3'] = rr(o['Us'], o['Ps'], o['Es'], nu)

    elif j == 5:
        n['Us'] = U * PH
        n['Ps'] = P + H * o['k3'][0]
        n['Es'] = E + H * o['k3'][1]

    elif j == 6:
        n['k4'] = rr(o['Us'], o['Ps'], o['Es'], nu)

    elif j == 7:
        n['d'] = (H / 6) * (
            o['k1'] + 2 * o['k2']
            + 2 * o['k3'] + o['k4']
        )

    elif j == 8:
        n['Um'] = U * PH2
        n['Pm'] = P + 0.5 * o['d'][0]
        n['Em'] = E + 0.5 * o['d'][1]

    elif j == 9:
        n['dphi'] = H * gN(
            o['Um'], o['Pm'], o['Em']
        )

    elif j == 10:
        n['U'] = U * PH
        n['P'] = P + o['d'][0]
        n['E'] = E + o['d'][1]

        n['H'] = o['H'] * complex(
            math.cos(o['dphi']),
            math.sin(o['dphi'])
        )

        n['Q'] = (
            o['Q']
            * math.exp(o['d'][2] / TAU0)
        )

    n['q'] = np.roll(o['q'], 1)
    return n
```

コード中の

```python
j = int(np.argmax(s['q']))
```

は one-hot 状態の非零成分を高速に選ぶ実装上の短縮である。数式上の相互作用定義は外部 loop counter ではなく

$$
\sum_{j=0}^{10}q_jF_j
$$

である。

一マクロstepは

```python
def macro_fast(s, nu):
    for _ in range(11):
        s = micro_fast(s, nu)
    return s
```

で検証した。ここでも `for` の回数を状態遷移の入力には使わない。作用選択に使われるのは $q$ だけである。

---

## 7.3 各内部相の数式

$h=\Delta\chi$、

$$
\Phi_{1/2}=e^{ih/2},
\qquad
\Phi_1=e^{ih}
$$

とする。

第1評価は

$$
K_1=R(U,P,E),
$$

第2評価点は

$$
U_2=U\Phi_{1/2},
$$

$$
P_2=P+\frac h2 K_{1p},
\qquad
E_2=E+\frac h2 K_{1e},
$$

第2評価は

$$
K_2=R(U_2,P_2,E_2).
$$

同様に

$$
U_3=U\Phi_{1/2},
\quad
P_3=P+\frac h2K_{2p},
\quad
E_3=E+\frac h2K_{2e},
$$

$$
K_3=R(U_3,P_3,E_3),
$$

$$
U_4=U\Phi_1,
\quad
P_4=P+hK_{3p},
\quad
E_4=E+hK_{3e},
$$

$$
K_4=R(U_4,P_4,E_4).
$$

合成は

$$
\boxed{
D
=\frac h6
(K_1+2K_2+2K_3+K_4)
}
$$

である。

角度評価点を

$$
U_m=U\Phi_{1/2},
\qquad
P_m=P+\frac12D_p,
\qquad
E_m=E+\frac12D_e
$$

とし、

$$
\Delta\phi=h\,g_N(U_m,P_m,E_m)
$$

から最終的に

$$
\boxed{
U'=U\Phi_1
}
$$

$$
\boxed{
P'=P+D_p,
\qquad
E'=E+D_e
}
$$

$$
\boxed{
H'=He^{i\Delta\phi}
}
$$

$$
\boxed{
Q'=Qe^{D_t/\tau_0}
}
$$

とする。

これは `P += Dp` のような逐次代入ではない。$P'$ と $E'$ は microstep 開始時の旧状態だけから一度に書かれる積和である。

---

# 8. 相互作用核の完全展開

## 8.1 $U$ から三角量を読む

$$
U=e^{i\chi}
$$

より

$$
C_\chi:=\cos\chi
=\frac{U+U^{-1}}2,
$$

$$
S_\chi:=\sin\chi
=\frac{U-U^{-1}}{2i}.
$$

したがって $\chi=\arg U$ を状態更新内で復元する必要はない。

---

## 8.2 2.5PN 放射反作用核

Paper 4 と同じ局所 2.5PN 放射反作用 [11,13] を用いる。

$$
d=1+EC_\chi,
\qquad
r=\frac Pd,
\qquad
\sqrt P=P^{1/2},
$$

$$
\dot r=\frac{ES_\chi}{\sqrt P},
\qquad
v_t=\frac d{\sqrt P},
$$

$$
v^2=\dot r^2+v_t^2,
\qquad
r^{-1}=\frac dP.
$$

放射反作用係数を

$$
C_{\rm RR}=\frac85\nu r^{-3},
$$

$$
A
=\dot r
\left(
18v^2+\frac23r^{-1}-25\dot r^2
\right),
$$

$$
B
=-\left(
6v^2-2r^{-1}-15\dot r^2
\right)
$$

とすると、

$$
\dot E_{\rm orb}
=C_{\rm RR}(A\dot r+Bv^2),
$$

$$
\dot p=2C_{\rm RR}Bp,
$$

$$
E_{\rm orb}=\frac{e^2-1}{2p},
$$

$$
\dot e
=\frac{p\dot E_{\rm orb}+E_{\rm orb}\dot p}{e}.
$$

また

$$
\frac{dt}{d\chi}
=\frac{r^2}{\sqrt p}.
$$

実際の次数検証コードの核は次である。

```python
def rr_numeric(U, P, E, nu=NU):
    Cc = 0.5 * (U + 1.0 / U)
    Ss = (U - 1.0 / U) / (2j)

    den = 1.0 + E * Cc
    r = P / den
    sqP = P**0.5

    rdot = E * Ss / sqP
    vt = den / sqP
    v2 = rdot * rdot + vt * vt
    invr = den / P

    Crr = 1.6 * nu * invr**3

    an = rdot * (
        18.0 * v2
        + (2.0 / 3.0) * invr
        - 25.0 * rdot * rdot
    )

    bv = -(
        6.0 * v2
        - 2.0 * invr
        - 15.0 * rdot * rdot
    )

    Edot = Crr * (an * rdot + bv * v2)
    pdot = 2.0 * Crr * bv * P

    Eorb = (E * E - 1.0) / (2.0 * P)
    edot = (P * Edot + Eorb * pdot) / E

    dt = r * r / sqP

    return pdot * dt, edot * dt, dt
```

この関数は $U,P,E$ だけを入力に取り、$dp/d\chi,de/d\chi,dt/d\chi$ を返す。

---

## 8.3 角度進行核

Schwarzschild 束縛軌道の Darwin パラメータ表示 [10] から

$$
\frac{d\phi}{d\chi}
=
\left[
1-\frac{2(3+e\cos\chi)}p
\right]^{-1/2}
$$

を用いる。

$$
u
:=
\frac{3+EC_\chi}{P}
$$

と置けば

$$
(1-2u)^{-1/2}
=
\sum_{m=0}^{\infty}c_mu^m,
$$

$$
c_0=1,
\qquad
c_{m+1}=\frac{2m+1}{m+1}c_m.
$$

$N$ 次で

$$
\boxed{
g_N(U,P,E)
=\sum_{m=0}^{N}c_m
\left[
\frac{3+\frac E2(U+U^{-1})}{P}
\right]^m
}
$$

である。

コードでは Horner 法で

```python
def gN_numeric(U, P, E, N):
    C = 0.5 * (U + 1.0 / U)
    u = (3.0 + E * C) / P

    c = COEFF[N]
    y = complex(c[-1])
    for a in c[-2::-1]:
        y = complex(a) + u * y
    return y
```

と評価する。

---

## 8.4 補助関数を消した完全展開の照合

途中の検討では `rr()` や `gN()` を残した段階を「完全展開」と誤って呼んだ。これは完全展開ではなかった。

最終監査では $dp/d\chi$ と $de/d\chi$ を $U^{\pm1}$ の有限 Laurent 多項式へ展開し、$g_{12}$ も係数を明示した上で、元式と 100,000 ランダム点で照合した。

最大誤差は

$$
|\Delta(dp/d\chi)|_{\max}
\approx1.25\times10^{-16},
$$

$$
|\Delta(de/d\chi)|_{\max}
\approx2.08\times10^{-17},
$$

$$
\left|\frac{\Delta(dt/d\chi)}{dt/d\chi}\right|_{\max}
\approx8.61\times10^{-16},
$$

$$
|\Delta g_{12}|_{\max}
\approx8.88\times10^{-16}
$$

であり、浮動小数点丸めの範囲で一致した。【数値】

---

# 9. 読み出し

## 9.1 軌道位置

状態 $U,P,E,H$ から

$$
C_\chi=\frac{U+U^{-1}}2
$$

を計算し、

$$
\boxed{
r=\frac{P}{1+EC_\chi}
}
$$

とする。

$$
H=e^{i\phi}
$$

なので

$$
\boxed{
x=r\operatorname{Re}H,
\qquad
y=r\operatorname{Im}H
}
$$

である。

---

## 9.2 時計

$$
Q=e^{t/\tau_0}
$$

なので観測時だけ

$$
\boxed{
t_{\rm obs}
=\tau_0\log |Q|
}
$$

と読む。

ここで得た $t_{\rm obs}$ を状態更新へ戻さない。

監査コード中の読み出し演算を関数として整理すると次である。これは監査で実際に用いた演算と同値である。

```python
def readout(U, P, E, H, Q, tau0):
    cos_chi = 0.5 * (U + 1.0 / U)
    r = (P / (1.0 + E * cos_chi)).real

    x = r * H.real
    y = r * H.imag

    # log is used only here, never in the state update
    t = tau0 * math.log(abs(Q))

    return x, y, r, t
```

したがって本稿の基本的な情報流は

$$
\boxed{
(U,P,E,H,Q)
\xrightarrow{\text{interaction}}
(U',P',E',H',Q')
}
$$

と

$$
\boxed{
(U,P,E,H,Q)
\xrightarrow{\text{readout}}
(x,y,r,t)
}
$$

に分離される。

---

# 10. B方式のプログラム監査結果

B方式について、単に最終軌道が似ているかではなく、状態機械としての規則を直接監査した。

検査項目は

1. 同一 microstep で更新後状態を読む箇所がないか、
2. one-hot 位相が壊れないか、
3. commit 相以外で永続五状態が変化しないか、
4. 11 microstep 後に位相が元へ戻るか、
5. 一マクロstepが元の Paper 4 RK4 step と一致するか、

である。

静的監査では更新後状態参照は 0 件であった。

500 個のランダム状態による一マクロstep比較では

$$
\boxed{
\max|\Delta(U,P,E,H,Q)|=(0,0,0,0,0)
}
$$

であった。

12 条件の長時間比較では

$$
\boxed{
\max|\Delta(U,P,E,H,Q)|
=
(0,0,0,3.38\times10^{-16},4.00\times10^{-15})
}
$$

であった。

また保存済み C1 CSV との比較では

$$
|\Delta p|_{\max}=4.97\times10^{-14},
$$

$$
|\Delta e|_{\max}=7.77\times10^{-16},
$$

$$
|\Delta t|_{\max}=1.94\times10^{-10},
$$

$$
|\Delta x|_{\max}=3.85\times10^{-12},
$$

$$
|\Delta y|_{\max}=5.11\times10^{-12}.
$$

![図4　12条件におけるB方式の厳密監査](figures/figure05_methodB_strict_heatmap.svg)

**図4.** 12条件での B方式と基準更新の差。0 の成分は表示下限へ切り下げている。非零差もほぼ浮動小数点丸めの範囲にある。

この結果は

$$
\boxed{
\text{Paper 4 の逐次処理を、隠れた更新途中値なしの明示状態機械へ再配置できる}
}
$$

ことを示す。

ただしこれは

$$
\boxed{
\text{この11相状態機械が自然界の基本作用である}
}
$$

ことを意味しない。

---

# 11. 近似次数を分離する

Paper 4 と本稿には複数の「次数」が存在する。これらを混同してはならない。

### 表2　近似次数の分類

| 層 | 現在または検証した次数 | 先頭の未保持項・誤差 | 意味 |
|---|---|---|---|
| 散逸物理 | leading 2.5PN | 高次 PN 項は未導入 | 物理モデルの近似次数 |
| 角度級数 | 基準 $N=12$、検証 $12,16,20,24$ | $O(u^{N+1})$ | Schwarzschild 角度因子の級数打切り |
| Paper 4 時間積分 | RK4 | 局所 $O(h^5)$、大域 $O(h^4)$ | 数値積分次数 |
| B 11相展開 | RK4 の代数的再配置 | 新しい打切り誤差なし | 内部状態機械の処理深度 |
| Taylor 次数検証 | $4,6,8$ | 実測収束を別途検査 | 同じ五状態での数値近似深度 |

特に

$$
\boxed{
\text{数値積分次数を上げること}
\neq
\text{PN物理次数を上げること}
}
$$

である。

また

$$
\boxed{
\text{近似次数を上げること}
\neq
\text{物理状態を追加すること}
}
$$

であるかを次に検証する。

---

# 12. 角度級数の高次化

角度因子の厳密式は

$$
g_\infty(u)=(1-2u)^{-1/2}
$$

であり、有限 $N$ では

$$
g_N(u)=\sum_{m=0}^{N}c_mu^m.
$$

状態 $Z=(U,P,E,H,Q)$ は固定し、$N$ だけを

$$
N=12,16,20,24
$$

へ増加させた。

最も厳しい試験条件の一つである

$$
p_0=24.3,
\qquad e_0=0.55
$$

では、3軌道後の方位角誤差が

| $N$ | $|\Delta\phi|$ |
|---:|---:|
| 12 | $1.56\times10^{-7}$ |
| 16 | $9.03\times10^{-10}$ |
| 20 | $5.54\times10^{-12}$ |
| 24 | $3.55\times10^{-14}$ |

と減少した。

![図5　角度級数の次数と誤差](figures/figure06_angular_order_convergence.svg)

**図5.** 9条件における $N=12,16,20,24$ の収束。状態数は変更せず、同じ級数生成則の打切り位置だけを変えている。

この結果から、少なくとも有限範囲では

$$
\boxed{
N\uparrow
\quad\text{while}\quad
\dim Z=5
}
$$

が成立する。

---

# 13. 五状態のまま数値積分次数を上げられるか

## 13.1 直接 Taylor 漸化

RK4 の段数を増やすと、内部 stage state も増える。

そこで「高次化のために新しい**永続物理状態**が必要か」を分離して調べるため、同じ五状態

$$
Z=(U,P,E,H,Q)
$$

に対して Taylor 係数を一 step 内だけで生成した。

連続ベクトル場を

$$
\frac{dZ}{d\chi}=F_N(Z)
$$

と書くと、

$$
Z(\chi+h)
=\sum_{k=0}^{p}A_k h^k+O(h^{p+1}),
$$

$$
A_0=Z_n,
$$

$$
\boxed{
(k+1)A_{k+1}
=[\xi^k]
F_N\!\left(\sum_{j=0}^{k}A_j\xi^j\right)
}
$$

とする。

コードの中核は次である。

```python
def taylor_step(z, h, order, N, tau0):
    # Same five persistent state variables.
    # No RK stage is promoted to a new persistent state.
    a = np.zeros((5, order + 1), complex)
    a[:, 0] = z

    for k in range(order):
        series = [TPS(a[j], order) for j in range(5)]
        U, P, E, H, Q = series

        dp, de, dt = rr_tps(U, P, E)
        g = gN_tps(U, P, E, N)

        F = [
            1j * U,
            dp,
            de,
            1j * g * H,
            (dt / tau0) * Q,
        ]

        for j in range(5):
            a[j, k + 1] = F[j].c[k] / (k + 1)

    powers = np.array(
        [h**k for k in range(order + 1)],
        complex
    )
    return a @ powers
```

ここで Taylor 係数 $A_k$ は一 step の計算係数であり、次のマクロstepへ保持しない。

したがって

$$
\boxed{
\text{persistent physical state dimension}=5
}
$$

を維持したまま $p=4,6,8$ を比較できる。

---

## 13.2 実測収束次数

9条件で step/orbit を 8、16、32、64 と変えた。

丸め誤差の影響が比較的小さい 16→32 step/orbit の収束率を 9条件平均で見ると、

### 表3　16→32 step/orbit の平均実測次数

| 設計次数 | $P$ | $E$ | $\phi$ |
|---:|---:|---:|---:|
| 4 | 3.998 | 3.993 | 3.928 |
| 6 | 5.997 | 6.010 | 5.921 |
| 8 | 7.918 | 8.072 | 7.916 |

となった。

4次・6次ではさらに 32→64 でもほぼ設計次数を保つ。8次では $P,E$ が機械精度へ近づき始めるため、より細かい step では参照解誤差と丸め誤差の影響が無視できなくなる。

![図6　同じ五状態での4・6・8次収束](figures/figure07_integration_order_convergence.svg)

**図6.** 同じベクトル場・同じ五つの永続状態に対して近似深度だけを変更した収束試験。

9条件の全体像を図7に示す。

![図7　9条件における高次化の比較](figures/figure08_integration_order_9case_summary.svg)

**図7.** 強場・中間・弱場および離心率の異なる条件でも、高次化に伴う誤差低下が確認できる。

したがって本稿で確認されたのは

$$
\boxed{
\text{近似次数の増加}
\not\Rightarrow
\text{新しい永続物理状態の追加}
}
$$

である。

ただし計算量や一時 Taylor 係数の個数は次数とともに増える。

---

# 14. 固定作用 $S$ と状態依存作用列について

本稿の監査では、現在の Paper 4 の表現をそのまま明示状態へ展開した場合、一つの固定作用だけでは処理を表せず、処理位相に応じた複数の内部作用を使用した。

数式上は

$$
\boxed{
\mathcal Z_{m+1}
=F_{\sigma_m}(\mathcal Z_m),
\qquad
\sigma_m\ \text{is encoded in }q_m
}
$$

である。

しかし、この結果から

$$
\boxed{
\text{自然界では必ず }S_n\text{ が必要}
}
$$

とは結論できない。

理由は三つある。

第一に、別の状態座標を選べば、一つの固定作用へ統合できる可能性がある。

第二に、内部11相は RK4 という数値アルゴリズムの処理構造を明示化したものであり、物理そのものの11段階性を証明するものではない。

第三に、Taylor 漸化では同じ五状態と一つのベクトル場 $F_N$ のまま高次化できた。

したがって本稿の結論は

$$
\boxed{
\text{現在の B 構成では作用列が必要だったが、その原理的必然性は未証明}
}
$$

である。

---

# 15. 本稿から分かったこと

## 15.1 「二状態モデル」への自己修正

Paper 4 の $X=(a,b)$ は有用な読み出し表現だったが、次状態を決める全情報を含む基本状態ではなかった。

これは Paper 4 の数値結果を無効にするものではない。

むしろ Paper 5 によって、その結果がどの状態情報に依存していたかが明示された。

---

## 15.2 隠れた状態を明示すると状態機械になる

通常の逐次数値コードでは、評価点や積分途中の値は「計算の都合」として扱われる。

しかし本稿の厳しい規則では、次の演算に使うなら情報を保持している以上、明示状態として扱う。

その結果、Paper 4 の更新は one-hot 位相を持つ有限状態機械として書き直せた。

---

## 15.3 読み出しと状態遷移を分離できる

特に時計について

$$
Q=e^{t/\tau_0}
$$

を状態とすれば、状態更新中に $t$ や $\log Q$ を必要としない。

$$
t=\tau_0\log|Q|
$$

は観測時だけに使える。

これは

$$
\boxed{
\text{内部状態の演算形式}
\neq
\text{観測者が読む座標形式}
}
$$

を具体的に示す。

---

## 15.4 高精度化と新しい物理を区別できる

角度 $N$ と数値積分次数を上げても、五つの永続状態は増えなかった。

したがって、少なくともこのモデルでは

$$
\boxed{
\text{高精度化}
\neq
\text{新しい物理自由度の追加}
}
$$

である。

これは今後、どの差が「不足している物理」で、どの差が「有限次数による数値誤差」かを区別するために重要である。

---

# 16. 限界と今後の検証

本稿には以下の限界がある。

1. **五状態の最小性は未証明である。**  
   別座標や保存量を利用すれば、より少ない永続状態で閉じる可能性がある。

2. **内部作用列の必然性は未証明である。**  
   B方式では11相が必要だったが、より自然な状態表現なら固定作用へ統合できる可能性がある。

3. **Paper 4 の物理近似をそのまま引き継ぐ。**  
   保存部分は Schwarzschild/test-particle 型、散逸は leading 2.5PN である [3,10-13]。完全な有限質量比一般相対論ではない。

4. **2.5PN より高い物理次数は検証していない。**  
   数値積分を8次へ上げても、物理モデル自身が2.5PNから高次化するわけではない。

5. **$U,H$ の複素表現は実装上の選択である。**  
   複素数を基本的実在と主張していない。実二成分回転表現との比較が今後必要である。

6. **内部 work state の最小性も未証明である。**  
   11相は RK4 をそのまま隠れ状態なしに展開した結果であり、最小内部状態機械とは限らない。

7. **厳密解への収束は証明していない。**  
   角度級数は既知の解析関数へ有限範囲で明確に収束したが、全物理モデルの無限次数極限が厳密 GR と一致することは別問題である。

今後の主要課題は、

$$
\boxed{
\text{五状態より小さい閉包が存在するか}
}
$$

および

$$
\boxed{
\text{内部作用列を一つの固定作用へ統合できるか}
}
$$

である。

---

# 17. 再現性と研究記録

本稿ではプログラム監査そのものが実験であるため、最終コードだけでなく途中経路も保存した。

Google Drive の研究フォルダ

```text
匿名頂点状態生成幾何-匿名内部観測者から不変な関係量の体系/
第五思考実験_乗法状態内部展開と近似次数検証_20260924/
```

に以下を保存している。

- `00_reference_from_paper4/` — Paper 4 の基準コードと C1 生データ
- `01_五状態同期写像の再構成/` — 五永続状態の監査
- `02_乗法状態表現と同期インライン化の経緯/` — 途中案・却下案を含む
- `03_A局所exp-log_vs_B内部状態展開/` — A/B比較
- `04_B方式_11相同期状態機械_厳密監査/` — B方式の静的・数値監査
- `05_演算ブロック完全展開_Laurent_N12/` — Laurent 展開と10万点照合
- `06_近似次数と状態_S不増加検証/` — $N=12,16,20,24$ および4,6,8次
- `figures/` — 論文用SVG
- `plot_scripts/` — 保存済みデータだけから図を再生成するコード
- `plot_data/` — 図表の元CSV

図の PNG は再生成可能なので Google Drive には保存せず、SVG と図化コードを保存した。

本稿の図化では物理計算を再実行していない。すべて既存 JSON / CSV の監査結果から生成している。

---

# 18. 結論

第四思考実験 [3] は、表面上は

$$
X_{n+1}=S_nX_n
$$

という二成分状態の相互作用から相対論的重力軌道を再構成する実験であった。

しかし本稿の自己監査により、実際の軌道生成には

$$
\chi,p,e,\phi,t
$$

という五つの永続状態が必要であり、$X=(a,b)$ はそれらから作られる派生表現であることが明らかになった。

したがって Paper 4 の結果をより正確に言い直すと、

$$
\boxed{
\text{二状態だけから軌道を生成したのではなく、}
\text{五つの永続状態から二状態表現を介して軌道を読み出していた。}
}
$$

となる。

次に、その五状態を明示し、RK4 の途中値も内部状態へ展開すると、Paper 4 の処理は隠れた更新途中値を使わない 11 相の同期状態機械として再構成できた。

500ランダム状態の一step照合では差が0、12条件の長時間照合でも最大差はほぼ機械精度であり、保存済み Paper 4 軌道も再現した。

また、角度級数を $N=12$ から $24$ へ上げても、数値積分次数を4から6、8へ上げても、永続状態

$$
Z=(U,P,E,H,Q)
$$

を増やす必要はなかった。

したがって、本稿で得た限定的な結論は

$$
\boxed{
\begin{array}{c}
\text{有限個の明示状態}\\
+\\
\text{現在状態だけを読む同期相互作用}\\
+\\
\text{状態へ戻さない読み出し}
\end{array}
\quad\Longrightarrow\quad
\text{Paper 3・4 と比較可能な重力軌道を再構成できる}
}
$$

である。

これは物理現象を解明したという主張ではない。

むしろ、背景時空、天体軌道、外部時計、隠れた逐次レジスタを基本状態として置かないという厳しい条件を課しても、既知の物理計算と同型の読み出しを構成できることを示した**モデルの表現能力の検証**である。

現段階で最も重要な未解決問題は、

$$
\boxed{
5\text{ が真の最小状態数か}
}
$$

と

$$
\boxed{
S_n\text{ が本質的に必要か、それとも一つの固定 }S\text{ へ統合できるか}
}
$$

である。

この二点を未証明のまま明示することが、次の思考実験の出発点になる。

---

# 参考文献

## 本研究系列

[0] 木原範昭, **状態と関係性から時空・粒子・相互作用の創発を探索する研究プロジェクト設計 ― 不確実性下の基礎物理モデル探索に対する予備設計と実装可能性検証の枠組み ―**, v1.0 (2026-09-20). Concept DOI: 10.5281/zenodo.22851944; Version DOI: 10.5281/zenodo.22851945.

[1] 木原範昭, **第二思考実験：クーロン型の逆二乗の力を、積の読み出しと面積の時計から導く ― 二つの値と線形の法則のまま、焦点・引力と斥力・中性を出し、二つの値では書けないものを確定する ―**, v1.2 (2026-09-21). Concept DOI: 10.5281/zenodo.22867336; Version DOI: 10.5281/zenodo.22876596.

[2] 木原範昭, **第三思考実験：a b 2体の天体の軌道から a b 2体の状態を読み出す ― 意味を定めない二つの値と G=c=M=1 の知識だけで、軌道と重力波の読み出しから質量比・相対距離・相対時間を得る ―**, v1.5 (2026-09-23). Concept DOI: 10.5281/zenodo.22909660; Version DOI: 10.5281/zenodo.22909661.

[3] 木原範昭, **第四思考実験：固定作用 $S$ から状態依存作用 $S_n$ へ ― 第三思考実験で外部計算していた相対論的重力軌道を、$X_{n+1}=S_nX_n$ の逐次相互作用から再生成できるか ―**, v1.0 (2026-09-23). Concept DOI: 10.5281/zenodo.22919787; Version DOI: 10.5281/zenodo.22919788.

## 外部文献：状態閉包・作用・同期更新・高次化

[4] R. Zwanzig, “Memory Effects in Irreversible Thermodynamics,” *Physical Review* **124**, 983–992 (1961). DOI: 10.1103/PhysRev.124.983.

[5] H. F. Trotter, “On the Product of Semi-Groups of Operators,” *Proceedings of the American Mathematical Society* **10**(4), 545–551 (1959). DOI: 10.1090/S0002-9939-1959-0108732-6.

[6] M. Suzuki, “Generalized Trotter's Formula and Systematic Approximants of Exponential Operators and Inner Derivations with Applications to Many-Body Problems,” *Communications in Mathematical Physics* **51**(2), 183–190 (1976). DOI: 10.1007/BF01609348.

[7] S. Wolfram, “Statistical Mechanics of Cellular Automata,” *Reviews of Modern Physics* **55**, 601–644 (1983). DOI: 10.1103/RevModPhys.55.601.

[8] H. Sambe, “Steady States and Quasienergies of a Quantum-Mechanical System in an Oscillating Field,” *Physical Review A* **7**, 2203–2213 (1973). DOI: 10.1103/PhysRevA.7.2203.

[9] H. Yoshida, “Construction of Higher Order Symplectic Integrators,” *Physics Letters A* **150**(5–7), 262–268 (1990). DOI: 10.1016/0375-9601(90)90092-3.

## 外部文献：Paper 4 から引き継ぐ重力モデル

[10] C. G. Darwin, “The Gravity Field of a Particle,” *Proceedings of the Royal Society of London A* **249**, 180–194 (1959). DOI: 10.1098/rspa.1959.0015.

[11] L. E. Kidder, “Coalescing Binary Systems of Compact Objects to (Post)$^{5/2}$-Newtonian Order. V. Spin Effects,” *Physical Review D* **52**, 821–847 (1995). DOI: 10.1103/PhysRevD.52.821.

[12] A. Pound and E. Poisson, “Osculating Orbits in Schwarzschild Spacetime, with an Application to Extreme Mass-Ratio Inspirals,” *Physical Review D* **77**, 044013 (2008). DOI: 10.1103/PhysRevD.77.044013.

[13] L. Blanchet, “Post-Newtonian Theory for Gravitational Waves,” *Living Reviews in Relativity* **27**, 4 (2024). DOI: 10.1007/s41114-024-00050-z.

---

# 付録A：論文本文で使用する図

本文で使用する SVG は、本稿と同じフォルダの `figures/` に保存した。

1. `figure00_hidden_state_audit_structure.svg`
2. `figure00b_interaction_constraints.svg`
3. `figure02_five_state_reconstruction_errors.svg`
4. `figure05_methodB_strict_heatmap.svg`
5. `figure06_angular_order_convergence.svg`
6. `figure07_integration_order_convergence.svg`
7. `figure08_integration_order_9case_summary.svg`

補足図として以下も保存している。

- `figure01_experiment_flow.svg`
- `figure03_full_inline_equivalence.svg`
- `figure04a_methodA_vs_B_internal.svg`
- `figure04b_methodA_vs_B_observables.svg`
- `figure05b_methodB_savedC1.svg`

---

# 付録B：主要な再現用ファイル

- `EXPERIMENT_INDEX.csv`
- `FIGURE_INDEX_ja.md`
- `FIGURE_TABLE_PLAN_ja.md`
- `REPRODUCIBILITY_POLICY_ja.md`
- `SHA256SUMS.txt`
- `FIGURE_SHA256SUMS.txt`
- `plot_scripts/generate_paper5_figures.py`
- `plot_scripts/generate_paper5_audit_summary.py`
- `plot_data/table01_state_audit_summary.csv`
- `plot_data/table02_approximation_orders.csv`
- `plot_data/table03_key_numerical_results.csv`

実験ごとのプログラム、導出式、結果 JSON/CSV、監査ログは各 `00`–`06` サブフォルダに保存している。
