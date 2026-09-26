# 荷電二体系・解析近似基準解 v1
## 計算結果と軌道・GW・EM 放射の分析

## 1. 基準条件

$$
m_A=m_B=0.5,\quad
\lambda_A=+0.30,\quad
\lambda_B=-0.30,
$$

$$
Z=1-\lambda_A\lambda_B=1.09,\quad
r:50M\rightarrow20M.
$$

本 v1 は leading-order Einstein–Maxwell adiabatic quasi-circular inspiral の解析近似基準であり、一般二体厳密解でも 1PN/2PN 精度の軌道でもない。

## 2. 主要数値

| 量 | 結果 |
|---|---:|
| 最終時刻 | $t_f=169032.317607\,M$ |
| 累積軌道位相 | $\phi_f=767.088999402$ rad |
| 軌道回数 | $122.086006046$ cycles |
| 初期角周波数 | $0.002952964612\,M^{-1}$ |
| 最終角周波数 | $0.011672617530\,M^{-1}$ |
| 初期 $P_{\rm EM}/P_{\rm GW}$ | $1.720183486$ |
| 最終 $P_{\rm EM}/P_{\rm GW}$ | $0.688073394$ |
| 累積 GW 放出エネルギー | $0.002028847173\,M$ |
| 累積 EM 放出エネルギー | $0.002058652837\,M$ |
| 累積総放出エネルギー | $0.004087500009\,M$ |
| 軌道エネルギー差からの期待値 | $0.004087500000\,M$ |
| 最大 energy-balance residual | $9.14\times10^{-12}M$ |

## 3. 軌道分析

`figure01_relative_orbit_spiral.svg/png` は $r=50M$ から $20M$ まで約 122.09 周する滑らかな inward spiral を示す。quasi-circular 条件なので近点・遠点振動や periapsis precession はこの v1 の検証対象ではない。

equal mass のため COM 系では

$$
|\mathbf x_A|=|\mathbf x_B|=r/2
$$

であり、`figure02_two_body_com_orbits.svg/png` では両物体が COM の反対側を保つ。

半径減衰は

$$
-\dot r=\frac{A}{r^2}+\frac{B}{r^3}
$$

なので、`figure03_radius_vs_time.svg/png` では後半ほど inspiral が加速する。

## 4. 重力波

leading mass-quadrupole power:

$$
P_{\rm GW}=\frac{32}{5}\frac{\mu^2Z^3}{r^5}.
$$

数値は

$$
P_{\rm GW}(50M)=1.65763712\times10^{-9},
$$

$$
P_{\rm GW}(20M)=1.61878625\times10^{-7},
$$

で、97.65625 倍へ増加する。

dominant GW frequency は

$$
f_{\rm GW}=2f_{\rm orb}=\frac{\omega}{\pi}.
$$

初期は $9.39957830\times10^{-4}M^{-1}$、最終は $3.71550956\times10^{-3}M^{-1}$ で、約 3.953 倍へ chirp する。

## 5. 電磁波

leading electric-dipole power:

$$
P_{\rm EM}=\frac23\frac{\mu^2(\Delta\lambda)^2Z^2}{r^4}.
$$

数値は

$$
P_{\rm EM}(50M)=2.85144\times10^{-9},
$$

$$
P_{\rm EM}(20M)=1.11384375\times10^{-7},
$$

で、39.0625 倍へ増加する。

dipole frequency は

$$
f_{\rm EM}=f_{\rm orb},
$$

したがって全 raw data で

$$
f_{\rm GW}=2f_{\rm EM}
$$

となる。

## 6. EM と GW の競合

解析的に

$$
\frac{P_{\rm EM}}{P_{\rm GW}}
=
\frac{5}{48}\frac{(\Delta\lambda)^2}{Z}r.
$$

よって

$$
r_{\rm cross}=\frac{48Z}{5(\Delta\lambda)^2}
=29.0666667M.
$$

raw data では $r\simeq29.06675M$, $t\simeq149537.149M$, 約 96.004 orbit で $P_{\rm EM}\simeq P_{\rm GW}$。

- $r>29.07M$: EM dipole 優勢
- $r<29.07M$: GW quadrupole 優勢

理由は

$$
P_{\rm EM}\propto r^{-4},\qquad
P_{\rm GW}\propto r^{-5}
$$

であり、接近すると GW の増加が一段速い。

## 7. 累積放出エネルギー

最終:

$$
E_{\rm GW}=0.002028847173M,
$$

$$
E_{\rm EM}=0.002058652837M.
$$

割合は

$$
\text{GW}\simeq49.64\%,\qquad
\text{EM}\simeq50.36\%.
$$

初期には EM power が GW の約 1.72 倍だが、後半で GW が急増するため累積値はほぼ半々になる。

`figure05_cumulative_radiated_energy.svg/png` では累積 GW+EM 放射と orbital binding-energy loss が数値精度内で一致する。

## 8. 独立検算

`validate_charged_binary_analytic_reference_v1.py` の結果:

- $dt/dr$ 閉形式の有限差分照合: 最大相対誤差 $9.80\times10^{-10}$
- $d\phi/dr$ 閉形式の有限差分照合: $9.68\times10^{-11}$
- energy balance から復元した $\dot r$: $5.55\times10^{-16}$
- Kepler identity: $4.44\times10^{-16}$
- $f_{\rm GW}=2f_{\rm orb}$: 保存精度内で誤差 0
- $f_{\rm EM}=f_{\rm orb}$: 保存精度内で誤差 0
- 最大 cumulative energy residual: $9.14\times10^{-12}$

## 9. charged 1PN 診断

trajectory 生成には使わず、Zhang et al. の charged 1PN Kepler law を truncation diagnostic としてのみ評価した。

今回

$$
\epsilon_{\rm PN}=Z/r=0.0218\rightarrow0.0545.
$$

固定 $\omega$ での 1PN radius correction は約

$$
-1.84\%\rightarrow-4.60\%,
$$

固定 $r$ での angular-frequency correction は約

$$
-3.05\%\rightarrow-7.82\%.
$$

したがって本 v1 を「1PN 精度の正解」と扱ってはならない。この診断は LO benchmark の truncation の大きさを測るだけであり、LO 軌道へ戻していない。

## 10. 後続実験との比較原則

後続の匿名状態相互作用系では、同じ初期条件から独立に生成した結果を

- $r(t)$
- $\phi(t)$
- $(x,y)$
- $\omega(t)$
- $P_{\rm GW}(t)$
- $P_{\rm EM}(t)$
- radiation frequency
- cumulative energy loss

で比較する。

比較方向は常に

$$
\boxed{
\text{anonymous generator output}
\longrightarrow
\text{reference comparison}
}
$$

であり、

$$
\boxed{
\text{reference error}
\not\longrightarrow
\text{generator update}
}
$$

を厳守する。

## 11. 限界

本 v1 は eccentric orbit、periapsis precession、generic 1PN/2PN conservative dynamics、NLO radiation、spin、tidal/finite-size、strong-field merger、full nonlinear Einstein–Maxwell numerical relativityを含まない。

必要なら v1 は固定保存したまま、別の v2 で高次 PN benchmark を追加する。
