# 複素ガウス Bob-star における κ の厳密分布

## 1. 定義

Bob-star の複素関係を

$$
z_1,\ldots,z_m\in\mathbb C
$$

とし、各成分は独立な円対称複素ガウスとする。

$$
\Re z_j,\Im z_j\sim\mathcal N(0,1).
$$

スケール不変量

$$
\kappa
=
\frac{\left|\sum_{j=1}^{m}z_j^2\right|}
{\sum_{j=1}^{m}|z_j|^2}
$$

の分布を求める。

絶対スケールは消えるので、

$$
w_j=\frac{z_j}{\sqrt{\sum_k|z_k|^2}}
$$

と置けば、$w=(w_1,\ldots,w_m)$ は複素単位球面 $S^{2m-1}$ 上で一様であり、

$$
\kappa=\left|\sum_j w_j^2\right|.
$$

---

## 2. 2×m 実ガウス行列への変換

各 $z_j=x_j+i y_j$ と書き、

$$
X=
\begin{pmatrix}
x_1&\cdots&x_m\\
y_1&\cdots&y_m
\end{pmatrix}
$$

を考える。

2×2 Gram 行列

$$
W=XX^{\mathsf T}
=
\begin{pmatrix}
A&C\\
C&B
\end{pmatrix}
$$

を置くと、

$$
A=\sum_jx_j^2,
\qquad
B=\sum_jy_j^2,
\qquad
C=\sum_jx_jy_j.
$$

一方、

$$
\sum_j z_j^2=(A-B)+2iC,
$$

$$
\sum_j|z_j|^2=A+B.
$$

したがって

$$
\kappa^2
=
\frac{(A-B)^2+4C^2}{(A+B)^2}.
$$

$W$ の固有値を $\lambda_1\ge\lambda_2\ge0$ とすれば、

$$
(\lambda_1-\lambda_2)^2=(A-B)^2+4C^2,
$$

$$
\lambda_1+\lambda_2=A+B.
$$

よって

$$
\boxed{
\kappa
=
\frac{\lambda_1-\lambda_2}{\lambda_1+\lambda_2}
}.
$$

つまり κ は 2×m real Wishart 行列の固有値非対称度そのものである。

---

## 3. κ の厳密密度

$W=XX^{\mathsf T}$ は自由度 $m$ の 2×2 real Wishart 行列であり、固有値の同時密度は定数因子を除いて

$$
f(\lambda_1,\lambda_2)
\propto
\exp\!\left[-\frac{\lambda_1+\lambda_2}{2}\right]
(\lambda_1\lambda_2)^{(m-3)/2}
(\lambda_1-\lambda_2)
$$

である。

変数を

$$
s=\lambda_1+\lambda_2,
\qquad
q=\frac{\lambda_1-\lambda_2}{\lambda_1+\lambda_2}
$$

へ変換する。

$$
\lambda_1=\frac{s}{2}(1+q),
\qquad
\lambda_2=\frac{s}{2}(1-q),
$$

$$
0\le q\le1,
$$

かつ Jacobian は

$$
\left|\frac{\partial(\lambda_1,\lambda_2)}{\partial(s,q)}\right|
=\frac{s}{2}.
$$

$q=\kappa$ であり、$s$ と $q$ の因子は分離する。

$q$ に依存する部分だけを残すと

$$
f_\kappa(q)
\propto
q(1-q^2)^{(m-3)/2}.
$$

規格化すると、$m\ge2$ で

$$
\boxed{
f_\kappa(q)
=(m-1)q(1-q^2)^{(m-3)/2},
\qquad
0\le q\le1
}.
$$

従って累積分布は

$$
\boxed{
F_\kappa(q)
=
1-(1-q^2)^{(m-1)/2}
}.
$$

さらに

$$
Y=\kappa^2
$$

と置けば、

$$
\boxed{
Y\sim\operatorname{Beta}\!\left(1,\frac{m-1}{2}\right)
}.
$$

$m=1$ は例外的に

$$
\kappa\equiv1
$$

である。

---

## 4. 平均値の厳密式

$m\ge2$ では

$$
E[\kappa]
=(m-1)\int_0^1 q^2(1-q^2)^{(m-3)/2}\,dq.
$$

$y=q^2$ と置けば

$$
E[\kappa]
=
\frac{m-1}{2}
B\!\left(\frac32,\frac{m-1}{2}\right).
$$

よって

$$
\boxed{
E[\kappa]
=
\frac{\sqrt\pi}{2}
\frac{\Gamma\!\left(\frac{m+1}{2}\right)}
{\Gamma\!\left(\frac{m+2}{2}\right)}
}.
$$

特に $m=2$ では

$$
\boxed{
E[\kappa]=\frac{\pi}{4}
}.
$$

---

## 5. 中央値の厳密式

中央値 $q_{1/2}$ は

$$
F_\kappa(q_{1/2})=\frac12
$$

を満たす。

従って

$$
(1-q_{1/2}^2)^{(m-1)/2}=\frac12,
$$

すなわち

$$
\boxed{
q_{1/2}
=
\sqrt{1-2^{-2/(m-1)}}
}.
$$

特に $m=2$ では

$$
\boxed{
\operatorname{median}(\kappa)
=\frac{\sqrt3}{2}
}.
$$

---

## 6. m=2 の幾何学的証明

$m=2$ では、正規化された $w=(w_1,w_2)\in S^3$ に対して Hopf 写像

$$
X=2\Re(w_1\bar w_2),
$$

$$
Y=2\Im(w_1\bar w_2),
$$

$$
Z=|w_1|^2-|w_2|^2
$$

を考える。

$(X,Y,Z)$ は $S^2$ 上で一様であり、

$$
X^2+Y^2+Z^2=1.
$$

また

$$
\kappa^2
=|w_1^2+w_2^2|^2
=1-Y^2.
$$

$S^2$ 一様分布では任意の直交座標成分 $Y$ は

$$
Y\sim U[-1,1]
$$

なので、

$$
\kappa=\sqrt{1-Y^2}.
$$

従って

$$
E[\kappa]
=
\frac12\int_{-1}^{1}\sqrt{1-y^2}\,dy
=
\frac\pi4.
$$

また

$$
P(\kappa\le q)
=
P\left(|Y|\ge\sqrt{1-q^2}\right)
=
1-\sqrt{1-q^2},
$$

中央値は

$$
1-\sqrt{1-q^2}=\frac12
$$

から

$$
q=\frac{\sqrt3}{2}.
$$

これは Wishart 導出の $m=2$ 特殊例である。

---

## 7. 全モーメント

任意の $p>-2$ に対して

$$
E[\kappa^p]
=
\frac{m-1}{2}
B\!\left(1+\frac p2,\frac{m-1}{2}\right).
$$

従って

$$
\boxed{
E[\kappa^p]
=
\frac{\Gamma\!\left(\frac{m+1}{2}\right)
\Gamma\!\left(1+\frac p2\right)}
{\Gamma\!\left(\frac{m+1+p}{2}\right)}
}.
$$

特に

$$
\boxed{
E[\kappa^2]=\frac{2}{m+1}
}.
$$

---

## 8. 大 m 漸近

厳密 CDF

$$
F_\kappa(q)=1-(1-q^2)^{(m-1)/2}
$$

で

$$
q=\frac{x}{\sqrt m}
$$

と置くと、

$$
(1-q^2)^{(m-1)/2}
\longrightarrow
\exp\left(-\frac{x^2}{2}\right).
$$

従って

$$
\boxed{
\sqrt m\,\kappa
\xrightarrow{d}
\operatorname{Rayleigh}(1)
}.
$$

特に平均は

$$
\boxed{
E[\kappa]
\sim
\sqrt{\frac{\pi}{2m}}
}.
$$

中央値は

$$
\boxed{
\operatorname{median}(\kappa)
\sim
\sqrt{\frac{2\ln2}{m}}
}.
$$

---

## 9. 低位 N 実験との比較

$m=N-1$ とする。

| N | m | 実験平均 κ | 厳密平均 κ | 実験中央値 κ | 厳密中央値 κ |
|---:|---:|---:|---:|---:|---:|
| 2 | 1 | 1.000000 | 1 | 1.000000 | 1 |
| 3 | 2 | 0.784638 | $\pi/4=0.785398$ | 0.866433 | $\sqrt3/2=0.866025$ |
| 4 | 3 | 0.664529 | $2/3=0.666667$ | 0.709719 | $1/\sqrt2=0.707107$ |
| 5 | 4 | 0.579951 | $3\pi/16=0.589049$ | 0.599575 | 0.608309 |
| 6 | 5 | 0.530283 | $8/15=0.533333$ | 0.538024 | 0.541196 |

全て有限標本誤差の範囲で厳密分布と整合する。

---

## 10. 結論

複素ガウス Bob-star に対する κ は、単なる数値的秩序量ではなく、任意の $m\ge2$ について厳密分布を持つ。

$$
\boxed{
\kappa^2
\sim
\operatorname{Beta}\!\left(1,\frac{m-1}{2}\right)
}
$$

この一式から、

- $m=2$ の $E[\kappa]=\pi/4$
- $m=2$ の中央値 $\sqrt3/2$
- 任意 $m$ の平均・中央値・全モーメント
- 大 $m$ の Rayleigh 漸近

がすべて統一的に導出される。

したがって N=2〜40 の統計実験では、generic 複素ガウス系の理論基準曲線が完全に既知となり、自己無撞着な実力学系との差を直接測定できる。
