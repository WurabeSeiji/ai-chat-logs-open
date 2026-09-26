# 導出式: 五状態同期写像

Paper 4 の逐次コードで永続していた量を

$$
X_n=(\chi_n,p_n,e_n,\phi_n,t_n)^T
$$

とする。RK4の中間量をすべて現在状態の関数として代入すれば、

$$
X_{n+1}=F_h(X_n)
$$

と1つの同期写像にできる。

放射反作用の局所率を

$$
R(\chi,p,e;\nu)=\left(\frac{dp}{d\chi},\frac{de}{d\chi},\frac{dt}{d\chi}\right)^T
$$

とすると、RK4は

$$
k_1=R(\chi,p,e),
$$
$$
k_2=R\left(\chi+\frac h2,p+\frac h2k_{1p},e+\frac h2k_{1e}\right),
$$
$$
k_3=R\left(\chi+\frac h2,p+\frac h2k_{2p},e+\frac h2k_{2e}\right),
$$
$$
k_4=R(\chi+h,p+hk_{3p},e+hk_{3e}),
$$
$$
d=\frac h6(k_1+2k_2+2k_3+k_4).
$$

したがって

$$
p_{n+1}=p_n+d_p,\quad e_{n+1}=e_n+d_e,\quad t_{n+1}=t_n+d_t.
$$

角度は中点

$$
p_m=p+\frac{d_p}{2},\qquad e_m=e+\frac{d_e}{2},\qquad \chi_m=\chi+\frac h2
$$

だけから

$$
\phi_{n+1}=\phi_n+h\,g_N(\chi_m,p_m,e_m)
$$

と書ける。よって右辺に更新後状態を残す必要はない。
