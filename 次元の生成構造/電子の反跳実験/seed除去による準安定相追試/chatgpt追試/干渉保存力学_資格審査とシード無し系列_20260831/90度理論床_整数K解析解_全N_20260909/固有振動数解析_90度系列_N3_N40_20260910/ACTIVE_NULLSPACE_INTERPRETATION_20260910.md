# 90度系列の実効振動自由度と零モード分解

作成日: 2026-09-10

## 1. 結論

90度系列の整数生成子は、辺を

- real sector: q_e even
- imaginary sector: q_e odd

に並べ替えると

$$
K=
\begin{pmatrix}
0 & C\\
-C^T & 0
\end{pmatrix}
$$

となる。

さらに C は unsigned incidence 行列を使って厳密に

$$
\boxed{
C=D_R\,B_R^T B_I\,D_I
}
$$

と因数分解できる。

ここで

- $B_R$: 同 parity 辺 $AA\cup BB$ の incidence 行列
- $B_I$: 異 parity 辺 $AB$ の incidence 行列
- $D_R,D_I$: 0/π, π/2/3π/2 の符号行列

である。

N=6..40 の全ケースでこの因数分解は浮動小数点誤差すらなく完全一致した。

---

## 2. rank(C)=N-1

$B_I$ は完全二部グラフ $K_{p,q}$ の unsigned incidence 行列である。

$$
p=\left\lceil\frac N2\right\rceil,\qquad
q=\left\lfloor\frac N2\right\rfloor.
$$

完全二部グラフは連結なので、符号変換で通常の oriented incidence 行列へ写せる。
従って

$$
\operatorname{rank}B_I=N-1.
$$

N>=6 では $B_R$ は $K_p\cup K_q$ の unsigned incidence で、各成分が非二部グラフとなるため

$$
\operatorname{rank}B_R=N.
$$

実測でも N=6..40 で

$$
\boxed{\operatorname{rank}C=N-1}
$$

が全て成立した。

従って

$$
\boxed{\operatorname{rank}K=2(N-1)}.
$$

---

## 3. これは N-1 個の独立振動モード

実反対称行列は直交変換により

$$
K\sim
\bigoplus_{j=1}^{r}
\begin{pmatrix}
0&\omega_j\\
-\omega_j&0
\end{pmatrix}
\oplus 0
$$

と標準化できる。

ここで

$$
r=\frac{\operatorname{rank}K}{2}=N-1.
$$

したがって 2(N-1) は独立モード数ではなく、

$$
\boxed{N-1\text{ 個の2次元回転モード}}
$$

を意味する。

各正固有値 $\omega_j$ が1個の独立な正常振動に対応する。

これは「N頂点のうち1自由度が全体ゲージとして消え、N-1個の相対自由度が残る」構造に非常に近い。

ただし現段階では、その1自由度を物理的な並進・時間ゲージと同一視してはいけない。数学的には $B_I$ の1次元左零空間、すなわち二部グラフのA/B符号ベクトルに由来する。

---

## 4. imaginary-sector の零モードは cycle space

imaginary sector は $AB$ 辺、すなわち完全二部グラフ $K_{p,q}$ の辺空間である。

$$
\dim E_I=pq.
$$

そして

$$
\ker C=\ker B_I.
$$

したがって

$$
\dim\ker C
=
pq-(N-1)
=
(p-1)(q-1).
$$

これは連結グラフの cycle-space dimension

$$
E-V+1
$$

そのものである。

従って imaginary-sector の巨大な零モード群は

$$
\boxed{\text{完全二部グラフ上の閉じた循環フロー}}
$$

と解釈できる。

これらは各頂点への net incidence をゼロにするため、生成子の active vertex space からは見えない。

---

## 5. real-sector の零モード

real sector は

$$
AA\cup BB
$$

の同 parity 辺で、

$$
\dim E_R=
\binom p2+\binom q2.
$$

その零モード次元は

$$
\dim\ker C^T
=
\binom p2+\binom q2-(N-1).
$$

N>=6 では $B_R$ が full row rank N なので、

- $\ker B_R$ に属する内部循環成分
- $B_R y$ が $B_I^T$ の1次元零方向へ写る追加1方向

に分けられる。

つまり real-sector の零モードも、基本的には「頂点へ見えない内部辺循環」とA/Bバランス方向から成る。

---

## 6. 全零モード次元

合計すると

$$
\dim\ker K
=
M-2(N-1)
$$

であり、

$$
\boxed{
\dim\ker K
=
\frac{(N-1)(N-4)}2
}
$$

となる。

N=40では

$$
M=780,\qquad
\operatorname{rank}K=78,\qquad
\dim\ker K=702.
$$

つまり 780本の関係波のうち、生成子が直接見る実効回転空間は78次元であり、702次元は零モードである。

---

## 7. 物理的示唆

この結果から直ちに「空間次元がN-1」とは言えない。

しかし少なくとも、

$$
\boxed{
O(N^2)\text{ の関係辺}
\rightarrow
N-1\text{ 個の正常振動モード}
}
$$

という非常に強い縮約が自己無撞着90度床から自然に出ている。

さらに零モードの大部分は cycle-space であり、辺レベルでは存在するが頂点レベルの active dynamics には現れない。

これは

- microscopic relation space
- observable/active mode space

を分ける数学的機構の候補になる。

---

## 8. 次の課題

最も重要なのは、N-1個の正常モードを

1. global mode
2. A/B sublattice mode
3. zero-cycle sector

へ整理したうえで、時間発展データの実際の振幅がどのモードに乗っているかを投影すること。

静的スペクトルだけでなく、既存の N=3..40 の時間発展 trajectory を正常モード基底へ射影すれば、

- どの固有振動が実際に励起されるか
- 奇数Nで低二族間のエネルギー交換が起きるか
- 偶数Nで一族へ固定されるか

を直接検証できる。
