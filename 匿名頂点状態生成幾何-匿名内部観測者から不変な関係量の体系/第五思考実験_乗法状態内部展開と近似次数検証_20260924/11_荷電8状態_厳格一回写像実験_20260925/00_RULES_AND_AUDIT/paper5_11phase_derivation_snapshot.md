# 導出式: 11相同期状態機械

内部状態は、永続5状態に加えて

$$
(k_1,k_2,k_3,k_4),\ (U_s,P_s,E_s),\ d,\ (U_m,P_m,E_m),\ \Delta\phi,\ q
$$

を持つ。

各microstepでは**旧microstateだけ**を読み、新しい全レジスタを同時に生成する。
例としてstage2は

$$
U_s'=Ue^{ih/2},\qquad
P_s'=P+\frac h2 k_{1p},\qquad
E_s'=E+\frac h2 k_{1e}.
$$

次のmicrostepで初めて

$$
k_2'=R(U_s,P_s,E_s)
$$

を評価する。

commit相も

$$
U'=Ue^{ih},\quad P'=P+d_p,\quad E'=E+d_e,
$$
$$
H'=He^{i\Delta\phi},\quad Q'=Qe^{d_t/\tau_0}
$$

を**旧状態から直接**計算する。
