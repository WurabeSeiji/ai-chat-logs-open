# 8A 数値結果要約

最小共面位相差結合

$$
K(\delta)=\begin{pmatrix}
1 & \cos\delta & \cos2\delta\\
\cos\delta & 1 & \cos\delta\\
\cos2\delta & \cos\delta & 1
\end{pmatrix}
$$

を用いて、n=1 の `aa,ab,bb` に対する GW/EM 放射振幅残差を走査した。

全ての評価半径 `r=50,40,30,20` で最小点は

$$
\boxed{\delta=\frac{2\pi}{3}}
$$

（対称解として $4\pi/3$ も存在）に現れた。

$\delta=2\pi/3$ では

$$
\cos\delta=\cos2\delta=-\frac12,
$$

したがって

$$
K=\begin{pmatrix}
1 & -1/2 & -1/2\\
-1/2 & 1 & -1/2\\
-1/2 & -1/2 & 1
\end{pmatrix},
$$

であり、

$$
K\begin{pmatrix}1\\1\\1\end{pmatrix}=0.
$$

数値走査でも GW と EM の両残差は double precision の丸め誤差水準まで低下した。

この結果の意味は限定的である。今回仮定した K が、`aa,ab,bb` の合成振幅を打ち消す離散位相 sector を持つことを確認しただけである。K の物理的妥当性、位相ロックの動力学的安定性、Paper 7 の局所状態写像へ交差項を feedback した場合の挙動は未検証であり、次段階の課題である。
