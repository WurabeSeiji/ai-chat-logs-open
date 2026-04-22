# 2+1次元ソリトンセルオートマトン

**森脇慎一*, 永井敦†, 薩摩順吉†, 時弘哲治†, 鳥居正城†, 高橋大輔‡, 松木平淳太‡**

\* 日本モトローラ株式会社
† 東京大学大学院 数理科学研究科
‡ 龍谷大学 応用数学情報学科

London Mathematical Society Lecture Note Series 255,
*Symmetries and Integrability of Difference Equations*, pp. 335–342 (1999)

---

## 概要

離散化された一般化戸田方程式から、超離散化と呼ばれる操作を通じて2+1次元ソリトンセルオートマトンを導出する。そのソリトン解と時間発展についても議論する。

---

## 1. はじめに

近年、離散ソリトン系が大きな注目を集めている。中でも、従属変数も独立変数も離散値をとる「超離散」ソリトン方程式が活発に研究されている。最も重要な超離散ソリトン系の一つは、いわゆる「ソリトンセルオートマトン」（SCA）である [1]。これは1（空間）+1（時間）次元で2値（0と1）である。時刻 $t$ における $j$ 番目のセルの値 $u_j^t$ の時間発展は次式で与えられる：

$$u_j^{t+1} = \begin{cases} 1 & u_j^t = 0 \text{ かつ } \sum_{i=-\infty}^{j-1} u_i^t > \sum_{i=-\infty}^{j-1} u_i^{t+1} \text{ のとき} \\ 0 & \text{それ以外} \end{cases} \tag{1}$$

式(1)の注目すべき性質は、あらゆる状態がソリトンのみから構成され、KdVソリトンと同様に相互作用することである。さらに、豊かな組合せ論的構造と無限個の保存量を持つ [2]。最近、SCAとKdV方程式の可積分離散化の一つであるLotka-Volterra方程式との直接的関係が明らかにされた [3, 4]。超離散化の鍵は次の公式である：

$$\lim_{\varepsilon \to +0} \varepsilon \log(1 + e^{X/\varepsilon}) = F(X) = \max[0, X]. \tag{2}$$

本稿の目的は、2+1次元戸田方程式の超離散化から導出される2+1次元SCAを提示することである。

---

## 2. 一般化戸田方程式の離散アナログと2+1次元SCA

広田 [5] が提案した以下の差分方程式から出発する：

$$[Z_1 \exp(D_1) + Z_2 \exp(D_2) + Z_3 \exp(D_3)] f \cdot f = 0, \tag{3}$$

ここで $Z_i$（$i = 1, 2, 3$）は任意のパラメータ、$D_i$ は未知関数 $f$ の変数に関する広田微分を表す。式(3)は一般化戸田方程式の離散アナログ（広田-三輪方程式）と呼ばれる。式(3)の適切な極限をとることで多くのソリトン方程式が得られる [5]。

本稿では、式(3)の特別な場合を考える：

$$f(t-1, x, y) f(t+1, x, y) - \delta^2 f(t, x-1, y) f(t, x+1, y) - (1 - \delta^2) f(t, x, y+1) f(t, x, y-1) = 0. \tag{5}$$

式(5)は独立・従属変数の変換を通じて、2+1次元戸田方程式の離散アナログに帰着する [7, 8]。

### 超離散極限

従属変数変換 $f(t, x, y) = \exp[S(t, x, y)]$ を施すと、

$$\exp[(\Delta_t^2 - \Delta_y^2) S] = (1 - \delta^2)\left(1 + \frac{\delta^2}{1 - \delta^2} \exp[(\Delta_x^2 - \Delta_y^2) S]\right) \tag{13}$$

が得られる。対数をとり $(\Delta_x^2 - \Delta_y^2)$ を作用させると、$u(t, x, y) = (\Delta_x^2 - \Delta_y^2) S(t, x, y)$ として、

$$(\Delta_t^2 - \Delta_y^2) u = (\Delta_x^2 - \Delta_y^2) \log\!\left(1 + \frac{\delta^2}{1 - \delta^2} \exp[u]\right). \tag{15}$$

$u = v_\varepsilon / \varepsilon$、$\delta^2/(1 - \delta^2) = e^{-\theta_0/\varepsilon}$ と置き、$\varepsilon$ の小さな極限をとると、**2+1次元SCA**の支配方程式が得られる：

$$(\Delta_t^2 - \Delta_y^2) v(t, x, y) = (\Delta_x^2 - \Delta_y^2) F(v(t, x, y) - \theta_0), \tag{18}$$

$$F(X) = \max[0, X]. \tag{19}$$

---

## 3. 2+1次元SCAのソリトン解

### 1-ソリトン解

双線形方程式(5)は1-ソリトン解 $f = 1 + e^\eta$（$\eta = px + qy + \omega t$）を持ち、分散関係を満たす。超離散極限をとると、新変数 $P = \varepsilon p$, $Q = \varepsilon q$, $\Omega = \varepsilon \omega$, $K = Px + Qy + \Omega t$ として：

$$v(t, x, y) = F(K + P) + F(K - P) - F(K + Q) - F(K - Q). \tag{25}$$

分散関係は：

$$|\Omega| = \max[|P|, |Q| + \theta_0] - \max[0, \theta_0]. \tag{26}$$

### 2-ソリトン解

2-ソリトン解は $f = 1 + e^{\eta_1} + e^{\eta_2} + e^{\eta_1 + \eta_2 + \theta_{12}}$ で与えられる。超離散極限をとると：

$$v(t, x, y) = \max[0, K_1 + P_1, K_2 + P_2, K_1 + K_2 + P_1 + P_2 + \Theta_{12}]$$
$$+ \max[0, K_1 - P_1, K_2 - P_2, K_1 + K_2 - P_1 - P_2 + \Theta_{12}]$$
$$- \max[0, K_1 + Q_1, K_2 + Q_2, K_1 + K_2 + Q_1 + Q_2 + \Theta_{12}]$$
$$- \max[0, K_1 - Q_1, K_2 - Q_2, K_1 + K_2 - Q_1 - Q_2 + \Theta_{12}] \tag{32}$$

位相シフト $\Theta_{12}$ は式(34)で決定される。

### N-ソリトン解

同様の極限操作によりN-ソリトン解も得られる：

$$v = \max_{\mu=0,1} \left[\sum_i \mu_i(K_i + P_i) + \sum_{i < j} \mu_i \mu_j \Theta_{ij}\right] + \max_{\mu=0,1} \left[\sum_i \mu_i(K_i - P_i) + \sum_{i < j} \mu_i \mu_j \Theta_{ij}\right]$$
$$- \max_{\mu=0,1} \left[\sum_i \mu_i(K_i + Q_i) + \sum_{i < j} \mu_i \mu_j \Theta_{ij}\right] - \max_{\mu=0,1} \left[\sum_i \mu_i(K_i - Q_i) + \sum_{i < j} \mu_i \mu_j \Theta_{ij}\right] \tag{35}$$

分散関係：$|\Omega_i| = \max[|P_i|, |Q_i| + \theta_0] - \max[0, \theta_0]$

---

## 4. 結語

2+1次元戸田方程式の超離散極限をとることにより、2+1次元SCAが得られることを示した。そのソリトン解も発見した。有理解、分子解、準周期解など他の種類の解の超離散版を構成することが今後の課題である。

---

## 付録A. 1-ソリトンおよび2-ソリトンの時間発展

パラメータ $P = 5$, $Q = 1$, $\theta_0 = 2$ の1-ソリトン解を図1に示す。パラメータ $P_1 = 6$, $Q_1 = 1$, $P_2 = 6$, $Q_2 = 5$, $\theta_0 = 2$ の2-ソリトン解を図2, 3に示す。

2-ソリトン解の時間発展（$t = -4$ から $t = 7$）では、2つのソリトンが衝突し位相シフトを伴って通過する様子が観察される。

---

## 参考文献

[1] D. Takahashi, J. Satsuma, J. Phys. Soc. Jpn. **59** (1990) 3514.

[2] M. Torii, D. Takahashi, J. Satsuma, Physica D **92** (1996) 209.

[3] D. Takahashi, J. Matsukidaira, Phys. Lett. A **209** (1995) 184.

[4] T. Tokihiro, D. Takahashi, J. Matsukidaira, J. Satsuma, Phys. Rev. Lett. **76** (1996) 3247.

[5] R. Hirota, J. Phys. Soc. Jpn. **45** (1981) 3785.

[6] R. Hirota, *ソリトンの数理*（岩波書店, 1992）.

[7] R. Hirota, M. Ito, F. Kako, Prog. Theor. Phys. Suppl. **94** (1988) 42.

[8] R. Hirota, S. Tsujimoto, T. Imai, RIMS講究録 **822** (1992) 144.
