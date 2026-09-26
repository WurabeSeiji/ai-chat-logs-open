# 導出式: RRと角度$N=12$の演算ブロック展開

ここで「完全展開」とは、これまで残していた `R()` / `gN()` という**名前付き演算ブロックを消す**という意味に限定する。加算、逆数、半整数冪、指数関数まで原始演算から除去したという意味ではない。

$U=e^{i\chi}$ とすると

$$
\cos\chi=\frac{U+U^{-1}}2,\qquad
\sin\chi=\frac{U-U^{-1}}{2i}.
$$

2.5PN局所RR式を代数展開すると

$$
\frac{dp}{d\chi}=-\frac{\nu}{P^{3/2}}\Bigg[
6E^3(U^3+U^{-3})+20E^2(U^2+U^{-2})
+\left(\frac{18}{5}E^3+\frac{112}{5}E\right)(U+U^{-1})
+\frac{56}{5}E^2+\frac{64}{5}\Bigg].
$$

$$
\begin{aligned}
\frac{de}{d\chi}=-\frac{\nu}{P^{5/2}}\Bigg[&
\frac54E^4(U^5+U^{-5})+\frac{55}{6}E^3(U^4+U^{-4})\\
&+\left(\frac{117}{20}E^4+\frac{316}{15}E^2\right)(U^3+U^{-3})\\
&+\left(22E^3+\frac{56}{3}E\right)(U^2+U^{-2})\\
&+\left(\frac52E^4+\frac{404}{15}E^2+\frac{32}{5}\right)(U+U^{-1})\\
&+\frac{121}{15}E^3+\frac{304}{15}E\Bigg].
\end{aligned}
$$

$$
\frac{dt}{d\chi}=
\frac{P^{3/2}}{\left[1+\frac E2(U+U^{-1})\right]^2}.
$$

角度級数は

$$
(1-2u)^{-1/2}=\sum_{m=0}^{\infty}c_m u^m,
\qquad c_{m+1}=\frac{2m+1}{m+1}c_m.
$$

$N=12$ の係数は

$$
1,1,\frac32,\frac52,\frac{35}{8},\frac{63}{8},\frac{231}{16},\frac{429}{16},
\frac{6435}{128},\frac{12155}{128},\frac{46189}{256},\frac{88179}{256},\frac{676039}{1024}.
$$
