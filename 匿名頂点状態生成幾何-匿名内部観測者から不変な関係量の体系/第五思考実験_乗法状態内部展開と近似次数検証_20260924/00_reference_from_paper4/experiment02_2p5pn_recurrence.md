# 実験02: 2.5PN放射反作用を S(n) に追加する最小漸化式

基準の S(n) は保存的で det S(n)=1 とする。
放射反作用は S の determinant を壊すのではなく、軌道要素 p_n, e_n のゆっくりしたドリフトとして入れる。

$$
X_{n+1}=S(n;p_n,e_n)X_n.
$$

論文3と同じ 2.5PN 放射反作用項を

$$
\mathbf a_{RR}=C(A\mathbf n+B\mathbf v)
$$

$$
C=\frac{8}{5}\nu r^{-3},
$$

$$
A=\dot r\left(18v^2+\frac{2}{3r}-25\dot r^2\right),
$$

$$
B=-\left(6v^2-\frac{2}{r}-15\dot r^2\right)
$$

とする。

Newtonian osculating elements では

$$
E=\frac{v^2}{2}-\frac1r,\qquad h=|\mathbf r\times\mathbf v|,
$$

$$
p=h^2,\qquad e^2=1+2Ep.
$$

したがって

$$
\dot E=\mathbf v\cdot\mathbf a_{RR}
=C(A\dot r+Bv^2),
$$

$$
\dot h=CBh,
$$

$$
\boxed{\dot p=2CBp},
$$

$$
\boxed{\dot e=\frac{p\dot E+E\dot p}{e}}.
$$

また

$$
r=\frac{p}{1+e\cos\chi},\qquad
\frac{dt}{d\chi}=\frac{r^2}{\sqrt p}
$$

なので

$$
\frac{dp}{d\chi}=\dot p\frac{dt}{d\chi},\qquad
\frac{de}{d\chi}=\dot e\frac{dt}{d\chi}.
$$

本実験ではこの2式を RK4 で更新し、更新後の p_n,e_n を同じ S(n) 生成式に渡す。
