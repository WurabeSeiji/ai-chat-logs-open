# 90度自己無撞着床の固有振動スペクトル：incidence低ランク因数分解による解析導出

作成日: 2026-09-10

## 1. 目的

辺空間の生成子 `K` はサイズ

$$
M\times M,\qquad M=N(N-1)/2
$$

だが、90度系列の全固有スペクトルは偶数Nで2族、奇数Nで3族にしかならない。
本メモでは、この事実を数値フィットではなく頂点 incidence 行列から解析的に導出する。

## 2. 一般の incidence 因数分解

完全グラフ K_N の unsigned incidence 行列を

$$
B\in\mathbb{R}^{N\times M},
$$

$$
B_{ve}=1 \quad (v\in e),\qquad 0\quad(\text{otherwise})
$$

とする。

辺 e の位相を $\theta_e$ とし、

$$
c_e=\cos\theta_e,\qquad s_e=\sin\theta_e
$$

を対角行列 $D_c,D_s$ に並べる。

完全グラフの辺隣接では、対角を除けば

$$
(B^TB)_{ef}=1
$$

は e,f が頂点を共有するとき、0 は共有しないときである。

したがって現行生成子

$$
K_{ef}=A_{ef}\sin(\theta_f-\theta_e)
$$

は、対角では正弦が0なのでそのまま

$$
K=D_c B^T B D_s-D_s B^T B D_c
$$

と書ける。

ここで

$$
X=
\begin{pmatrix}
BD_c\\
BD_s
\end{pmatrix}\in\mathbb{R}^{2N\times M},
\qquad
J=
\begin{pmatrix}
0&I\\
-I&0
\end{pmatrix}
$$

とすれば

$$
\boxed{K=X^T J X}.
$$

これは90度床に限らない、現行 `sin(位相差)` 生成子の厳密因数分解である。

従って辺空間 M が $O(N^2)$ であっても、非零スペクトルは最大でも $2N$ 次元の頂点拡張空間で決まる。

## 3. 90度床での簡約

90度床では

$$
\theta_e\in\{0,\pi/2,\pi,3\pi/2\}.
$$

各辺は実軸群または虚軸群のどちらか一方に属し、各辺で $c_es_e=0$。

ラベル parity は

$$
q_e\bmod2=(a+b)\bmod2
$$

なので、頂点を偶奇二群 A,B に分けると

- 実軸側（even label）: 同 parity 辺 = AA と BB
- 虚軸側（odd label）: 異 parity 辺 = AB

となる。

実軸辺の incidence を $B_R$、虚軸辺の incidence を $B_I$ とすると

$$
XX^T=
\begin{pmatrix}
G_R&0\\
0&G_I
\end{pmatrix},
\qquad
G_R=B_RB_R^T,\quad G_I=B_IB_I^T.
$$

非零固有値は `AB/BA` の標準事実より $JXX^T$ で求められ、その二乗は

$$
(JXX^T)^2=
-\begin{pmatrix}
G_IG_R&0\\
0&G_RG_I
\end{pmatrix}.
$$

したがって

$$
\boxed{\omega^2\text{ は }G_IG_R\text{ の非零固有値}}
$$

であり、M×M問題が N×N問題にまで落ちる。

## 4. G_R と G_I の閉形式

$$
p=|A|,\qquad q=|B|,\qquad p+q=N.
$$

### 同 parity 辺（AA,BB）

A側は完全グラフ $K_p$、B側は完全グラフ $K_q$ なので unsigned incidence Gram は

$$
G_R=
\begin{pmatrix}
(p-2)I_p+J_p&0\\
0&(q-2)I_q+J_q
\end{pmatrix}.
$$

### 異 parity 辺（AB）

AB辺は完全二部グラフ $K_{p,q}$ なので

$$
G_I=
\begin{pmatrix}
qI_p&J_{p\times q}\\
J_{q\times p}&pI_q
\end{pmatrix}.
$$

ここで $J$ は全成分1の行列。

## 5. 不変部分空間による完全対角化

頂点空間は三つに分解できる。

### A零和部分空間

A上で成分和0、B上0。次元 $p-1$。

この空間では

$$
G_R=(p-2)I,\qquad G_I=qI
$$

なので

$$
\boxed{\omega_A^2=q(p-2)},\qquad \text{multiplicity }p-1.
$$

### B零和部分空間

B上で成分和0、A上0。次元 $q-1$。

$$
\boxed{\omega_B^2=p(q-2)},\qquad \text{multiplicity }q-1.
$$

### A/B定数部分空間

A上一定値 $\alpha$、B上一定値 $\beta$ の2次元空間。

この空間で $G_IG_R$ は

$$
2
\begin{pmatrix}
q(p-1)&q(q-1)\\
p(p-1)&p(q-1)
\end{pmatrix}
$$

となり rank 1。

固有値は

$$
0,
\qquad
4pq-2N.
$$

したがって全体モードは

$$
\boxed{\omega_G^2=4pq-2N},\qquad \text{multiplicity }1.
$$

もう1方向は零モードである。

## 6. 偶数N

偶数Nでは

$$
p=q=N/2.
$$

A零和族とB零和族が完全縮退し、

$$
\omega_A^2=\omega_B^2
=\frac N2\left(\frac N2-2\right)
=\boxed{\frac{N(N-4)}4}.
$$

多重度は

$$
(p-1)+(q-1)=N-2.
$$

全体モードは

$$
\omega_G^2=4(N/2)^2-2N
=\boxed{N(N-2)}.
$$

従って正固有値は2族のみ。

## 7. 奇数N

奇数Nでは

$$
p=\frac{N+1}{2},\qquad q=\frac{N-1}{2}.
$$

A零和族:

$$
\omega_A^2=q(p-2)
=\boxed{\frac{(N-1)(N-3)}4},
$$

多重度

$$
p-1=\boxed{\frac{N-1}{2}}.
$$

B零和族:

$$
\omega_B^2=p(q-2)
=\boxed{\frac{(N+1)(N-5)}4},
$$

多重度

$$
q-1=\boxed{\frac{N-3}{2}}.
$$

全体モード:

$$
\omega_G^2=4pq-2N
=N^2-1-2N
=\boxed{N^2-2N-1},
$$

多重度1。

これで、N=7以上の奇数系列で数値的に見つかった3族の固有値・多重度がすべて解析的に導出される。

## 8. 奇偶差の本質

奇偶差は「奇数だから特殊」という抽象的なものではなく、

$$
\boxed{|A|-|B|=N\bmod2}
$$

という頂点 parity 二群の個数差そのもの。

- 偶数N: $p=q$ → A/B零和モードが縮退 → 低周波1族 + global 1族 = 2族
- 奇数N: $p=q+1$ → A/B零和モードの縮退が解ける → 低周波2族 + global 1族 = 3族

したがって、先に観測した奇数Nの二サブ格子モードは、この完全グラフの parity 分割から必然的に生じる。

## 9. ランク

正固有値の多重度和は常に

$$
(p-1)+(q-1)+1=N-1.
$$

反対称Kでは ± 固有値が対になるので

$$
\boxed{\operatorname{rank}K=2(N-1)}
$$

（一般系列 N>=5。N=3,4 は小N特殊縮退あり）。

従って辺空間の零モード多重度は

$$
M-2(N-1)
=\boxed{\frac{(N-1)(N-4)}2}.
$$

これは「M=O(N^2) の関係空間に対し、動的非零モードは O(N)」という低ランク性の解析的起源を与える。

## 10. 物理的読みの範囲

ここで証明したのは90度自己無撞着床の線形生成子スペクトルの構造である。

FEM的には、

- A零和族
- B零和族
- A/B一様 global 族

という固有振動族に分かれる。

奇数Nでは A/B の1頂点不均衡が低周波縮退を分裂させる。偶数Nでは縮退が復活する。

これが非線形時間発展で観測した奇偶差の直接原因かどうかは別途検証が必要であり、現段階では「線形スペクトル側に奇偶差の明確な構造的原因がある」とまで主張する。
