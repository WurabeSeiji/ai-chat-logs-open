# RWZ spin-2 波方程式 再実行・分析 v2

木原範昭 研究系列 / 2026-09-29

## 1. 目的

重力波側について、quadrupole 放射近似ではなく Schwarzschild 背景上の実際の spin-2 master equation である Regge-Wheeler / Zerilli-Moncrief 方程式を直接時間発展し、同一 master field から無限遠放射と horizon 吸収を同時に読み出せることを確認する。

本実験は full nonlinear binary Einstein solution そのものではない。背景時空は exact Schwarzschild、波動部分はその上の exact linear spin-2 perturbation equation である。したがって、今回確認したのは「spin-2 放射・吸収 sector への対応」であり、一般二体非線形重力場の完成ではない。

## 2. 使用方程式

Schwarzschild tortoise coordinate r_* 上で

\[
\left[-\partial_t^2+\partial_{r_*}^2-V_{\ell}(r)\right]\Psi_{\ell m}=0.
\]

odd parity には Regge-Wheeler potential、even parity には Zerilli-Moncrief potential を用いた。

Regge-Wheeler:

\[
V_{\rm RW}=f\left[\frac{\ell(\ell+1)}{r^2}-\frac{6M}{r^3}\right],\qquad
f=1-\frac{2M}{r}.
\]

Zerilli:

\[
V_{\rm Z}=
\frac{2f\left[n^2(n+1)r^3+3n^2Mr^2+9nM^2r+9M^3\right]}
{r^3(nr+3M)^2},
\qquad n=\frac{(\ell-1)(\ell+2)}2.
\]

今回は最低 gravitational radiative multipole \(\ell=2\) を使用した。

## 3. フラックス

Moncrief master field の規格化に対応して

\[
C_\ell=\frac{(\ell+2)!}{64\pi(\ell-2)!}
\]

を用い、局所エネルギー流を

\[
S=-C_\ell\,\partial_t\Psi\,\partial_{r_*}\Psi
\]

として左右境界から積分した。

左境界を horizon 側、右境界を infinity 側として

\[
E_H=\int F_Hdt,\qquad E_\infty=\int F_\infty dt
\]

を同じ master field から得る。

## 4. 再実行条件

\[
M=1,\quad \ell=2,\quad r_*\in[-120,220],
\]

\[
N_x=6801,\quad \Delta r_*=0.05,\quad
\Delta t=0.0225,\quad t_{\max}=420.
\]

初期波束は

\[
\Psi(0,r_*)=10^{-3}\exp\left[-\frac{(r_*-5)^2}{2\,7^2}\right],
\qquad \partial_t\Psi(0,r_*)=0.
\]

境界は horizon 側を ingoing、infinity 側を outgoing Sommerfeld 条件とした。

## 5. Regge-Wheeler odd \(\ell=2\) 再実行結果

\[
E_0=8.113272173738477\times10^{-8},
\]

\[
E_H=2.7961289218009255\times10^{-8},
\]

\[
E_\infty=5.3169762084753656\times10^{-8}.
\]

したがって

\[
\frac{E_H}{E_0}=0.344636400939636,
\]

\[
\frac{E_\infty}{E_0}=0.6553430101464698.
\]

残留エネルギー比は

\[
5.429102740277031\times10^{-6}
\]

で、エネルギー閉包は

\[
\frac{E_{\rm final}+E_H+E_\infty}{E_0}
=0.9999848401888461.
\]

## 6. Zerilli-Moncrief even \(\ell=2\) 再実行結果

\[
E_0=8.003711626513542\times10^{-8},
\]

\[
E_H=2.5351047493045664\times10^{-8},
\]

\[
E_\infty=5.468444288102532\times10^{-8}.
\]

したがって

\[
\frac{E_H}{E_0}=0.31674114056116626,
\]

\[
\frac{E_\infty}{E_0}=0.683238545225375.
\]

閉包は

\[
0.9999851085592671.
\]

## 7. 放射 = 吸収バランス探索

同じ方程式・同じ potential のまま、初期 Gaussian の中心 \(x_0\) のみを変え、

\[
E_\infty-E_H=0
\]

となる点を二分探索した。

### odd / Regge-Wheeler

\[
x_0\simeq2.39013671875.
\]

\[
\frac{E_H}{E_0}=0.49997667352320063,
\qquad
\frac{E_\infty}{E_0}=0.49998528792115576.
\]

差は

\[
\frac{E_\infty-E_H}{E_0}
=8.614397955131192\times10^{-6}.
\]

### even / Zerilli-Moncrief

\[
x_0\simeq1.89697265625.
\]

\[
\frac{E_H}{E_0}=0.5001262773373583,
\qquad
\frac{E_\infty}{E_0}=0.4998357001699545.
\]

差は

\[
\frac{E_\infty-E_H}{E_0}
=-2.905771674037716\times10^{-4}.
\]

探索分解能を上げればさらに詰められるが、この時点で同じ spin-2 master equation の中に「外向き放射と horizon 吸収がほぼ等量になる状態」が存在することを数値的に確認した。

## 8. 本研究系への意味

重要なのは、放射と吸収を別々の経験式で追加していないことである。

\[
\boxed{
\Psi_{\ell m}
\rightarrow
\{E_H,E_\infty\}
}
\]

という一つの spin-2 状態から、同時に両フラックスが生成される。

したがって本研究で要求していた

\[
\boxed{
\text{spin-2 state}
\rightarrow
\text{emission + absorption}
}
\]

という構造には対応できた。

さらに、ほぼ

\[
E_H=E_\infty
\]

となる初期状態も存在するため、「放射ゼロ」ではなく

\[
\boxed{
\text{radiation}\neq0,\qquad
\text{absorption}\neq0,\qquad
\text{balance}\simeq0
}
\]

という動的バランス状態を RWZ sector 内で表現できる。

## 9. aa / ab / bb との接続

前段で確認した二次関係状態

\[
\mathcal Q=\{a^2,2ab,b^2\}
\]

はいずれも共通回転で spin-weight 2 型に変換される。

今回の RWZ master field は物理的 spin-2 gravitational perturbation を伝播させるので、次の接続課題は

\[
\{a^2,2ab,b^2\}
\longrightarrow
\Psi^{\rm RW/Z}_{2m}
\]

を同じ有限履歴型万能相互作用から生成し、各 sector の和ではなく full master field の一状態として進めることである。

## 10. 現在の到達点と未対応部分

対応済み:

- exact Schwarzschild background
- Regge-Wheeler odd spin-2 propagation
- Zerilli-Moncrief even spin-2 propagation
- horizon absorption
- infinity radiation
- energy closure
- emission/absorption balance point search

未対応:

- Kerr background の Teukolsky equation
- spin-weighted spheroidal harmonics
- Kerr horizon superradiance
- full nonlinear dynamical Einstein-Maxwell two-body evolution
- finite-history harmonic persistent state と RWZ/Teukolsky solver の直接結合

したがって「GR重力波に未対応」という段階は脱したが、万能相互作用全体が完成したわけではない。次段階は Kerr/Teukolsky と、aa/ab/bb finite-history state からの source/readout の直結である。
