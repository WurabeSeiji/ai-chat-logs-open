# 荷電二体系・解析近似基準解 v1
## 解析式と初期条件

**目的:** 後続の匿名状態相互作用系とは完全に独立な「答え合わせ側」の基準軌道を構成する。  
**重要:** 本ファイルの式・読み出し値・誤差は、匿名状態相互作用系の生成則へ一切フィードバックしない。

---

## 1. 採用した近似

本 v1 は、Einstein–Maxwell 二体系の一般厳密解ではない。次の条件を同時に課した、**leading-order adiabatic quasi-circular inspiral** である。

- 二つの非回転点粒子。
- 幾何単位・自然単位:
  $$
  G=M=c=1,\qquad \epsilon_0=\frac{1}{4\pi}.
  $$
- 総質量:
  $$
  M=m_A+m_B=1.
  $$
- 保存力学は Newtonian 重力 + Coulomb の leading order。
- 軌道は準円軌道で、放射反作用による半径変化は一周の運動より十分遅いとする（adiabatic）。
- 重力波は leading mass-quadrupole 放射。
- 電磁波は leading electric-dipole 放射。
- 放射損失は energy-balance で軌道エネルギーへ結びつける。
- 軌道の時間発展は時間刻みで数値積分せず、$t(r)$ と $\phi(r)$ を閉形式で評価する。
- charged 1PN Kepler 則は近似誤差の診断にのみ用い、v1 軌道の生成へ戻さない。

これは一般の Einstein–Maxwell 二体問題を解いたものではなく、後続の反復相互作用系を検証するための「条件を明示した解析近似基準」である。

---

## 2. 文献上の基礎

主参照:

Zi-Han Zhang, Tan Liu, Shuai Zhang, Zong-Kuan Guo,  
**Post-Newtonian dynamics of charged compact binaries**, arXiv:2604.10654 (2026).

この論文は Einstein–Maxwell 作用から荷電二体系を 1PN で展開し、準円軌道に対して GW/EM の flux balance を扱う。

本 v1 で直接参照した定義は

$$
Q=\frac{q_Aq_B}{Gm_Am_B},
\qquad
Z=1-Q,
$$

および 1PN 診断用の

$$
\zeta_0
=
\frac{q_A^2}{m_AM}
+
\frac{q_B^2}{m_BM}
$$

である（本ファイルでは $G=M=c=1$）。

文献の 1PN circular Kepler law は診断専用として

$$
\omega_{\rm 1PN}^2
=
\frac{1}{r^3}
\left[
Z+\frac{F_1}{r}
\right],
$$

$$
F_1
=
\eta-3
-\frac12(4\eta+1)Q^2
+\left(\frac92+\eta\right)Q
-\zeta_0,
$$

と用いた。

---

## 3. 質量・電荷の無次元化

charge-to-mass ratio を

$$
\lambda_A=\frac{q_A}{m_A},
\qquad
\lambda_B=\frac{q_B}{m_B}
$$

と置く（$G=1$）。

換算すると

$$
q_A=\lambda_A m_A,
\qquad
q_B=\lambda_B m_B,
$$

なので

$$
Q=\lambda_A\lambda_B,
\qquad
Z=1-\lambda_A\lambda_B.
$$

換算質量は

$$
\mu=\frac{m_Am_B}{M},
\qquad
\eta=\frac{m_Am_B}{M^2}.
$$

$M=1$ では

$$
\mu=\eta.
$$

---

## 4. leading-order の準円軌道

重力 + Coulomb の有効中心引力から、

$$
\boxed{
\omega^2=\frac{Z}{r^3}
}
$$

とする。

したがって

$$
v_{\rm rel}=r\omega=\sqrt{\frac{Z}{r}},
$$

leading PN parameter は

$$
\epsilon_{\rm PN}\sim v^2=\frac{Z}{r}.
$$

軌道束縛エネルギーは

$$
\boxed{
E_{\rm orb}
=
-\frac{\mu Z}{2r}
}
$$

である。

---

## 5. 重力波放射

leading mass-quadrupole power は

$$
P_{\rm GW}
=
\frac{32}{5}\mu^2 r^4\omega^6.
$$

$\omega^2=Z/r^3$ を代入すると

$$
\boxed{
P_{\rm GW}
=
\frac{32}{5}
\frac{\mu^2Z^3}{r^5}
}
$$

となる。

dominant quadrupole radiation frequency は

$$
\boxed{
f_{\rm GW}=2f_{\rm orb}
=\frac{\omega}{\pi}
}.
$$

---

## 6. 電磁波放射

COM 系での electric dipole moment は leading order で

$$
\mathbf d
=
\mu(\lambda_A-\lambda_B)\mathbf r.
$$

ここで

$$
\Delta\lambda=\lambda_A-\lambda_B.
$$

円運動では

$$
|\ddot{\mathbf d}|
=
\mu|\Delta\lambda|r\omega^2.
$$

$\epsilon_0=1/(4\pi),c=1$ の Larmor dipole formula

$$
P_{\rm EM}
=
\frac23|\ddot{\mathbf d}|^2
$$

より、

$$
P_{\rm EM}
=
\frac23
\mu^2(\Delta\lambda)^2r^2\omega^4.
$$

したがって

$$
\boxed{
P_{\rm EM}
=
\frac23
\frac{\mu^2(\Delta\lambda)^2Z^2}{r^4}
}
$$

となる。

leading electric-dipole の放射周波数は軌道周波数なので

$$
\boxed{
f_{\rm EM}=f_{\rm orb}
=\frac{\omega}{2\pi}
}.
$$

---

## 7. energy balance と半径減衰

放射損失を

$$
\frac{dE_{\rm orb}}{dt}
=
-\left(P_{\rm GW}+P_{\rm EM}\right)
$$

とする。

$$
\frac{dE_{\rm orb}}{dr}
=
\frac{\mu Z}{2r^2}
$$

なので

$$
\frac{dr}{dt}
=
-\frac{P_{\rm GW}+P_{\rm EM}}
{dE_{\rm orb}/dr}.
$$

整理すると

$$
\boxed{
\frac{dr}{dt}
=
-\frac{A}{r^2}
-\frac{B}{r^3}
=
-\frac{Ar+B}{r^3}
}
$$

ただし

$$
\boxed{
A=
\frac43\mu(\Delta\lambda)^2Z
}
$$

は EM dipole 項、

$$
\boxed{
B=
\frac{64}{5}\mu Z^2
}
$$

は GW quadrupole 項である。

---

## 8. $t(r)$ の閉形式

$$
\frac{dt}{dr}
=
-\frac{r^3}{Ar+B}.
$$

次の原始関数を定義する:

$$
F_t(r)
=
\frac{r^3}{3A}
-\frac{Br^2}{2A^2}
+\frac{B^2r}{A^3}
-\frac{B^3}{A^4}\ln(Ar+B).
$$

これは

$$
\frac{dF_t}{dr}
=
\frac{r^3}{Ar+B}
$$

を満たす。

初期半径 $r_0$ で $t=0$ とすれば

$$
\boxed{
t(r)=F_t(r_0)-F_t(r)
}.
$$

したがって軌道の時間は time-step ODE 積分ではなく解析式から直接得る。

---

## 9. $\phi(r)$ の閉形式

$$
\frac{d\phi}{dt}=\omega=\sqrt{\frac{Z}{r^3}}
$$

と $dt/dr$ から

$$
\frac{d\phi}{dr}
=
-\frac{\sqrt{Z}\,r^{3/2}}{Ar+B}.
$$

$u=\sqrt r$ を用いると原始関数は

$$
F_\phi(r)
=
\sqrt Z
\left[
\frac{2r^{3/2}}{3A}
-\frac{2B\sqrt r}{A^2}
+\frac{2B^{3/2}}{A^{5/2}}
\tan^{-1}\sqrt{\frac{Ar}{B}}
\right].
$$

よって $\phi(r_0)=0$ として

$$
\boxed{
\phi(r)=F_\phi(r_0)-F_\phi(r)
}.
$$

相対軌道は

$$
\boxed{
x_{\rm rel}=r\cos\phi,
\qquad
y_{\rm rel}=r\sin\phi
}.
$$

COM 座標は

$$
\mathbf x_A=\frac{m_B}{M}\mathbf r,
\qquad
\mathbf x_B=-\frac{m_A}{M}\mathbf r.
$$

---

## 10. GW と EM の競合

比を取ると

$$
\frac{P_{\rm EM}}{P_{\rm GW}}
=
\frac{5}{48}
\frac{(\Delta\lambda)^2}{Z}
r.
$$

したがって

$$
\boxed{
r_{\rm cross}
=
\frac{48Z}{5(\Delta\lambda)^2}
}
$$

で

$$
P_{\rm EM}=P_{\rm GW}
$$

となる。

EM dipole は $r^{-4}$、GW quadrupole は $r^{-5}$ なので、同一条件では遠方ほど EM が相対的に強く、内側へ入るほど GW が追い越す。

---

## 11. v1 の初期条件

後続の匿名相互作用実験で同じ無次元条件を使えるよう、基準条件を固定した。

$$
M=1,
\qquad
m_A=m_B=0.5,
$$

$$
\eta=\mu=0.25,
$$

$$
\lambda_A=+0.30,
\qquad
\lambda_B=-0.30,
$$

$$
q_A=+0.15,
\qquad
q_B=-0.15,
$$

$$
Q=-0.09,
\qquad
Z=1.09,
$$

$$
\Delta\lambda=0.60,
\qquad
\zeta_0=0.09.
$$

軌道条件:

$$
r_0=50M,
\qquad
r_f=20M,
\qquad
\phi_0=0.
$$

raw data は半径方向に 40001 点保存する。

この設定を選んだ理由:

1. $Z>0$ で安定な leading-order attractive circular benchmark を構成できる。
2. $\lambda_A\ne\lambda_B$ のため electric dipole radiation が非零。
3. equal mass + equal-and-opposite charge により比較条件を対称化できる。
4. $r=50M\to20M$ では
   $$
   \epsilon_{\rm PN}=Z/r=0.0218\to0.0545
   $$
   で、leading-order benchmark として truncation を監視できる。
5. この範囲内で EM/GW power crossover が発生するため、両放射チャネルを同一実験で検証できる。

---

## 12. 1PN 診断は生成式に混ぜない

1PN charged Kepler law を用いた補正量を raw CSV に保存するが、

$$
\boxed{
\text{1PN diagnostic}
\not\rightarrow
\text{LO trajectory generation}
}
$$

である。

これは「v1 がどの程度 leading-order 近似なのか」を測るための答え合わせであり、解析式を途中で都合よく補正するためではない。

後続の匿名相互作用系に対しても同様に、

$$
\boxed{
\text{reference branch}
\not\rightarrow
\text{generator branch}
}
$$

を厳守する。

---

## 13. 出力の役割

- `solve_charged_binary_analytic_reference_v1.py`  
  上記閉形式を評価し raw CSV と metadata を生成する。
- `charged_binary_analytic_reference_v1_raw.csv`  
  一次データ。図化プログラムはこのファイルだけを読む。
- `charged_binary_analytic_reference_v1_metadata.json`  
  初期条件、解析係数、式、主要結果、raw CSV SHA-256。
- `validate_charged_binary_analytic_reference_v1.py`  
  微分恒等式、energy balance、Kepler identity を独立検算する。
- `plot_charged_binary_analytic_reference_v1.py`  
  raw CSV のみから SVG/PNG を生成する。軌道を再計算しない。

