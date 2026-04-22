# 対称テンソルの結晶に関連する $A_M^{(1)}$ オートマトン

**畑山泰佑\*, 日置健二†, 井上玲†, 国場敦夫\*, 高木太一郎‡, 時弘哲治§**

\* 東京大学大学院 総合文化研究科 物理学教室
† 東京大学大学院 理学系研究科 物理学教室
‡ 防衛大学校 数学物理学教室
§ 東京大学大学院 数理科学研究科

arXiv:math/9912209v1 [math.QA] (1999年12月27日)

---

## 概要

量子アフィン代数 $U'_q(A_M^{(1)})$ の対称テンソル表現の結晶に関連するソリトンセルオートマトンを導入する。これは、箱とキャリアの容量が任意かつ非一様である一般化箱玉系の結晶理論的定式化である。2つのソリトンの散乱行列は $U'_q(A_{M-1}^{(1)})$ の組合せ論的R行列と一致する。オートマトンの区分線形発展方程式は、非自律離散KP方程式の超離散極限として同定される。後者のソリトン解の超離散化を通じて、N-ソリトン解のクラスが得られる。

---

## 1. はじめに

高橋・薩摩 [TS] により発明された箱玉系は、ソリトンセルオートマトンの重要な例である。これは、ある規則の下で有限個のボールが一次元の箱の列に沿って移動する離散力学系である。その可積分性は [TTMS] において、超離散化と呼ばれる極限操作を通じてLotka-Volterra方程式の差分アナログ [H1] との関係を確立することにより証明された。

元の箱玉系はいくつかの方向に一般化されている。第一に、集合 $\{1, 2, \ldots, M\}$ からのインデックスで区別されるボールを導入できる。第二に、サイト $n$ の箱が $\theta_n$ 個までのボールを収容できるようにし、容量 $\theta_n$ を $n$ に依存させる。第三に、容量 $\kappa_t$ のキャリアを導入して時刻 $t$ の時間発展を再定義できる。キャリアは左から来て右へ進み、ある規則の下で箱からボールを拾い上げ別の箱に落とす。箱の列を通過する間に、連続的な積み降ろし過程がボールの運動、したがって系の時間発展を誘導する。

これらの一般化はパラメータ $(M, \theta_n, \kappa_t)$ で特徴付けられる。元の箱玉系 [TS] は $(M, \forall\theta_n, \forall\kappa_t) = (1, 1, \infty)$ に対応する。

本稿の目的は、一般の $(M, \theta_n, \kappa_t)$ の場合を研究することである。

---

## 2. 結晶からのオートマトン

### 2.1 $U'_q(A_M^{(1)})$ 結晶

$B_k$ を、$U_q(A_M)$ の $k$ 重対称テンソル表現に対応する $U'_q(A_M^{(1)})$ の古典的結晶とする。集合として、文字 $\{1, 2, \ldots, M+1\}$ 上の長さ $k$ の一行半標準タブローからなる：

$$B_k = \{ \boxed{m_1 \cdots m_k} \mid m_i \in \{1, \ldots, M+1\},\ m_1 \leq \cdots \leq m_k \}$$

内容の多重度による表記も用いる。すなわち $b = \boxed{m_1 \cdots m_k} \in B_k$ は $b = (x_1, x_2, \ldots, x_{M+1})$（$x_i = \#\{l \mid m_l = i\}$）とも表す。

柏原作用素 $\tilde{f}_i$, $\tilde{e}_i$（$i = 0, 1, \ldots, M$）の作用は $b = (x_1, \ldots, x_{M+1}) \in B_k$ に対して：

$$\tilde{e}_0 b = (x_1 - 1, x_2, \ldots, x_{M+1} + 1), \quad \tilde{f}_0 b = (x_1 + 1, x_2, \ldots, x_{M+1} - 1)$$
$$\tilde{e}_i b = (\ldots, x_i + 1, x_{i+1} - 1, \ldots), \quad \tilde{f}_i b = (\ldots, x_i - 1, x_{i+1} + 1, \ldots) \quad (i = 1, \ldots, M)$$

2つの結晶 $B$, $B'$ のテンソル積 $B \otimes B'$ を定義する。集合として $B \otimes B' = \{b_1 \otimes b_2 \mid b_1 \in B, b_2 \in B'\}$。柏原作用素の作用は：

$$\tilde{e}_i(b_1 \otimes b_2) = \begin{cases} \tilde{e}_i b_1 \otimes b_2 & \varphi_i(b_1) \geq \varepsilon_i(b_2) \text{ のとき} \\ b_1 \otimes \tilde{e}_i b_2 & \varphi_i(b_1) < \varepsilon_i(b_2) \text{ のとき} \end{cases}$$

$B' \otimes B$ と $B \otimes B'$ は標準的に同型であり、同型写像 $R : B' \otimes B \xrightarrow{\sim} B \otimes B'$ を**組合せ論的R行列**と呼ぶ。

### 2.2 同型写像

$R : B_k \otimes B_l \to B_l \otimes B_k$ を結晶グラフ全体を描くことなく得る明示的手続きを与える。

**命題2.3.** 同型写像 $R$ を得る規則：
1. $k \geq l$ のとき右列のドットを選び、左列でパートナーを見つける（巻き付き/巻き戻しの線）
2. 残りの未接続ドットについて繰り返す
3. 未対となった $(k-l)$ 個のドットを右にスライドさせる

### 2.3 オートマトン

正整数の2つの列 $\{\theta_n\}$（箱の容量）と $\{\kappa_t\}$（キャリアの容量）を考える。2次元格子の各点 $(t, n)$ で、上辺に $b_n^t \in B_{\theta_n}$、左辺に $v_n^t \in B_{\kappa_t}$ を配置する。組合せ論的R行列により：

$$R : v_n^t \otimes b_n^t \xrightarrow{\sim} b_n^{t+1} \otimes v_{n+1}^t$$

配列 $\{b_n^t\}$ から $\{b_n^{t+1}\}$ への発展は、$B_{\kappa_t}$ を補助空間とする $q = 0$ 行転送行列の作用として実現される。これを**$A_M^{(1)}$ オートマトン**と呼ぶ。

時間発展は可逆かつ可換：$T_\kappa T_{\kappa'} = T_{\kappa'} T_\kappa$。

**注意2.6.** 「空間」と「時間」の役割を交換して式(8)を見ることができる。

**例2.7.** $M = 3$, $\forall\theta_n = 1$, $\forall\kappa_t = \infty$ の場合：

```
···111142113111111111111···
···111111421311111111111···
···111111114231111111111···
···111111111124311111111···
···111111111111243111111···
···111111111111121143111···
```

典型的な2ソリトン散乱。振幅 $l$ のソリトンは速度 $l$ で右に移動。大きなソリトンが小さなものを追い越す。

### 2.4 箱玉系との等価性

$A_M^{(1)}$ オートマトンは一般化箱玉系として解釈できる。タブロー中の文字1を空箱、$2 \leq i \leq M+1$ の文字をインデックス $M+2-i$ のボールと解釈する。

$\forall\kappa_t = \infty$ のとき、$\{b_n^t\}$ の発展は [TTM] で研究された箱玉系と等価。力学系の規則：
1. すべてのボールを1回だけ動かす
2. インデックス1の最左のボールをその右の最近傍の空きのある箱へ移す
3. 同様に残りを繰り返す
4. インデックス1のボールがすべて移動したら、インデックス2に対して同じ手続きを行う
5. インデックス $M$ まで繰り返す

時間発展 $T_\infty$ は $T_\infty = \tilde{T}_M \cdots \tilde{T}_2 \tilde{T}_1$ と分解される。

---

## 3. 超離散ソリトンの散乱行列としての組合せ論的R行列

### 3.1 ソリトン

$B'_k$ を、$U'_q(A_{M-1}^{(1)})$ の $k$ 重対称テンソル表現に対応する古典的結晶とする：

$$B'_k = \{ \boxed{m_1 \cdots m_k} \mid m_i \in \{1, \ldots, M\},\ m_1 \leq \cdots \leq m_k \}$$

各 $k \in \mathbb{Z}_{\geq 1}$ に対して写像 $\iota_k$ を定義する：

$$\iota_k : B'_k \to (B_1)^{\otimes k}, \quad \boxed{m_1 \cdots m_k} \mapsto \boxed{m_k + 1} \otimes \cdots \otimes \boxed{m_1 + 1}$$

**漸近N-ソリトン状態**：ソリトン間の距離が十分大きい状態。各ソリトンは $B'_k$ の要素に関連付けられ、$k$ をソリトンの**振幅**と呼ぶ。

**命題3.3.** 振幅 $l$ の1-ソリトン $\mathbf{p} = \iota_l^{(L_0, L_1)}(b)$ の時間発展 $T_\kappa(\mathbf{p})$ は再び1-ソリトンであり、同じ $b \in B'_l$ に関連付けられる（内部状態が保存される）。速度は：

$$x(T_\kappa(\mathbf{p})) - x(\mathbf{p}) = \begin{cases} \kappa & \kappa < y(\mathbf{p}) \text{ のとき} \\ \min(\kappa, l) + \max(\theta_{n+k} - l, 0) & \kappa \geq y(\mathbf{p}) \text{ のとき} \end{cases}$$

### 3.2–3.3 2-ソリトン散乱

2つのソリトンの散乱行列は $U'_q(A_{M-1}^{(1)})$ の組合せ論的R行列 $R' : B'_l \otimes B'_k \to B'_k \otimes B'_l$ と一致する。

**定理3.10.** 振幅 $l, k$（$l > k$）のソリトンの散乱において：
- 散乱前の内部状態 $b \otimes c \in B'_l \otimes B'_k$ が
- 散乱後に $\tilde{c} \otimes \tilde{b} \in B'_k \otimes B'_l$ に変わる
- この写像は組合せ論的R行列 $R'$ に等しい

### 3.4 保存量

$A_M^{(1)}$ オートマトンはエネルギー関数から定義される保存量を持つ。テンソル積のエネルギー関数 $H_{B_l B_k}$ は、巻き戻し線の数の $(-1)$ 倍に等しい。

---

## 4. 区分線形発展方程式と非自律離散KP方程式

### 4.1 組合せ論的R行列の区分線形方程式

**命題4.1.** $R : B_\kappa \otimes B_\theta \xrightarrow{\sim} B_\theta \otimes B_\kappa$ について、$R(v \otimes b) = b' \otimes v'$ を $v = (v_0, \ldots, v_{M+1})$, $b = (b_0, \ldots, b_{M+1})$ 等として区分線形方程式で表す：

$$\begin{cases} v'_i = v_i - W_i + W_{i-1} \\ b'_i = b_i + W_i - W_{i-1} \end{cases}$$

$$W_i = \min(v_0 + \cdots + v_i,\ b_{i+1} + \cdots + b_{M+1}), \quad W_0 = 0, \quad W_{M+1} = \min(\kappa, \theta)$$

### 4.2 発展方程式と非自律離散KP方程式

$A_M^{(1)}$ オートマトンの発展方程式を $\tau$ 関数で表現し、非自律離散KP（ndKP）方程式の超離散極限として同定する。

**定理4.5.** ndKP方程式

$$\tau_n^{t+1} \tau_{n+1}^t - \tau_{n+1}^{t+1} \tau_n^t = (\delta_n - \delta_{n+1}) \tau_n^{t+1} \tau_{n+1}^{t+1}$$

の超離散極限が $A_M^{(1)}$ オートマトンの発展方程式を与える。

### 4.3 N-ソリトン解

ndKP方程式の $\tau$ 関数のソリトン解の超離散化を通じて、N-ソリトン解のクラスが得られる。各オートマトンのソリトンは、ndKPの $M$ 個のソリトンが超離散極限で合体して得られる。

---

## 5. まとめ

一般の $(M, \theta_n, \kappa_t)$ に対する $A_M^{(1)}$ オートマトンを結晶理論により定式化した。散乱行列が $U'_q(A_{M-1}^{(1)})$ の組合せ論的R行列と一致することを証明した。発展方程式はndKP方程式の超離散極限として同定され、$\tau$ 関数を通じてN-ソリトン解が構成された。

---

## 参考文献（主要なもの）

- [TS] D. Takahashi, J. Satsuma, J. Phys. Soc. Jpn. **59** (1990) 3514. — 箱玉系の原論文
- [TTMS] T. Tokihiro, D. Takahashi, J. Matsukidaira, J. Satsuma, Phys. Rev. Lett. **76** (1996) 3247. — 超離散化の確立
- [KMN1] S.-J. Kang, M. Kashiwara, K. C. Misra, T. Miwa, T. Nakashima, A. Nakayashiki, Int. J. Mod. Phys. A **7** (1992) 449. — 完全結晶
- [NY] A. Nakayashiki, Y. Yamada, Selecta Mathematica, New Series **3** (1997) 547. — 区分線形方程式
- [K] M. Kashiwara, Duke Math. J. **63** (1991) 465. — 結晶基底
- [FOY] G. Hatayama, K. Hikami, R. Inoue, A. Kuniba, T. Takagi, T. Tokihiro, Contemp. Math. **297** (2002) 151. — 結晶理論的アプローチ
- [HKT] G. Hatayama, A. Kuniba, T. Takagi, Nucl. Phys. B **577** (2000) 619. — 一般化BBS
- [TTM] T. Tokihiro, D. Takahashi, J. Matsukidaira, J. Math. Phys. **40** (1999) 2–. — 非一様BBS
