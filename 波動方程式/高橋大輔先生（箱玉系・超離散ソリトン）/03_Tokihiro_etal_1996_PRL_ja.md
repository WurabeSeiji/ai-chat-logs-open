# ソリトン方程式から可積分セルオートマトンへ ── 極限手続きによる導出

**T. 時弘哲治¹, 高橋大輔², 松木平淳太², 薩摩順吉¹**

¹ 東京大学大学院 数理科学研究科, 東京 153
² 龍谷大学 応用数学情報学科, 瀬田, 大津 520-21

(1995年12月26日受理)

Physical Review Letters, Vol. 76, No. 18, pp. 3247–3250 (1996年4月29日)

---

## 概要

セルオートマトンと可積分非線形波動方程式との直接的な関係を示す。また、セルオートマトンのN-ソリトン公式を提示する。最後に、このような可積分セルオートマトンとそのN-ソリトン解を構成する一般的方法を提案する。

PACS番号: 03.40.Kf, 02.10.Jf, 52.35.Sb

---

## 本文

複雑な物理現象を調べるためには、現象の本質的な特徴を示す単純なモデルを採用する必要がある。セルオートマトン（CA）はそのような単純なモデルとして機能し、物理学、化学、生物学、計算機科学において広く研究されてきた [1]。CAはその時間発展において多様な構造を呈する。約10年前、Park, Steiglitz, Thurstionはその中からソリトンの概念を抽出した [2]。彼らは、フィルタ型のCAにおいて、非ゼロのセル値のあるパターンがしばしば一定の有限速度で伝播し、衝突後にもその恒等性を保持することを発見した。これらの振る舞いは、Korteweg–de Vries（KdV）方程式のような非線形波動方程式の孤立波解（ソリトン）と極めて類似している。ソリトン構造を持つ種々のCAが研究され、系の多くの特徴が明らかにされてきたが [3–5]、我々の知る限り、これらのCAと非線形波動方程式との直接的な関係はこれまで報告されていなかった。

数年前、著者のうち2名（高橋、薩摩）が新しい型のフィルタCAを提案した [6]。このCAは1（空間）+1（時間）次元で、2値（0と1）である。時刻 $t$ における $j$ 番目のセルの値 $u_j^t$ は次のように与えられる：

$$u_j^{t+1} = \begin{cases} 1 & \text{if } u_j^t = 0 \text{ かつ } \sum_{i=-\infty}^{j-1} u_i^t > \sum_{i=-\infty}^{j-1} u_i^{t+1}, \\ 0 & \text{それ以外.} \end{cases} \tag{1}$$

ここで $|j| \gg 1$ に対して $u_j^t = 0$ を仮定する。発展規則を別の言い方で述べることもできる [7]。時刻 $t$ において、0と1からなる無限列がある。1の個数は有限である。時刻 $t+1$ の状態を決定する規則は以下の通りである：(1) すべての1を1回だけ動かす。(2) 最も左の1をその最も近い右の0と交換する。(3) 残りの1のうち最も左のものをその最も近い右の0と交換する。(4) すべての1が移動するまでこの手続きを繰り返す。

このCAの特異な性質は、あらゆる状態がソリトンのみから構成され、KdVソリトンと同じ様式で相互作用することである（図1）。さらに、無限個の保存量を持つ [8]。したがって、可積分非線形波動方程式のアナログと考えられていた。本論文の目的は、CA [式(1)] と可積分非線形波動方程式との関係を明らかにすることにより、CA（または離散値をとる差分-差分方程式）の一つのクラスと可積分非線形波動方程式との直接的関係を示すことである。これにより、このような可積分CAを構成する一般的方法を提示できる。N-ソリトン解の明示的な形も提示する。

**図1.** 式(1)の時間発展の例。1のパターン3つ（1111, 11, 1）が衝突後に位相シフトを伴いつつ形を保持する。

---

## KdV方程式からCAへの道

最もよく知られた非線形波動方程式の一つはKdV方程式である：

$$\frac{\partial}{\partial t} a(x,t) = \frac{\partial^3}{\partial x^3} a(x,t) + a(x,t) \frac{\partial}{\partial x} a(x,t). \tag{2}$$

KdV方程式はLax表現を許容し、ハミルトン構造を持ち、逆散乱変換により厳密に解ける意味で可積分である [9]。

KdV方程式の可積分離散化（微分差分方程式）はLotka-Volterra方程式である [10]：

$$\frac{d}{dt} b_j(t) = b_j(t) \left[ b_{j+1}(t) - b_{j-1}(t) \right]. \tag{3}$$

これは生物種の生存競争のモデルやプラズマ中のLangmuir振動の薄い構造に現れる。$b_j(t) = 1 + (1/6)\varepsilon^2 a((j + 2t)\varepsilon, \varepsilon^3 t/3)$ と置き、$\varepsilon \to 0$ の極限をとると、式(3)はKdV方程式に変換される。Lotka-Volterra方程式も可積分であり、無限個の保存量を持つ [10,11]。

Lotka-Volterra方程式に関連する可積分な差分-差分方程式が広田・辻本により提案された [12]：

$$\frac{c_j^{t+1}}{c_j^t} = \frac{1 + \delta c_{j-1}^t}{1 + \delta c_{j+1}^{t+1}} \quad (j, t \in \mathbf{Z}). \tag{4}$$

この方程式もN-ソリトン解と無限個の保存量を持つ。$c_j^t = b_j(-\delta t)$ と置き $\delta \to 0$ とすれば、Lotka-Volterra方程式が式(4)の連続極限であることは容易に確認できる。$c_j^t$ を $\exp(d_j^t)$ で置き換えると、

$$d_j^{t+1} - d_j^t = \ln\!\left(\frac{1 + \delta \exp(d_{j-1}^t)}{1 + \delta \exp(d_{j+1}^{t+1})}\right). \tag{5}$$

ここで、可積分CAを得るための鍵となる重要な極限操作を行う。正のパラメータ $\varepsilon \equiv -(\ln \delta)^{-1}$（すなわち $\delta = e^{-1/\varepsilon}$）を導入し、$d_j^t = e_j^t / \varepsilon$ と置く。すると、

$$F(X) \equiv \lim_{\varepsilon \to +0} \varepsilon \ln(1 + e^{X/\varepsilon}) = \max[0, X] \tag{6}$$

という事実に注意すれば、$\varepsilon \to +0$ の極限で式(5)から

$$e_j^{t+1} - e_j^t = -F(e_{j+1}^{t+1} - 1) + F(e_{j-1}^t - 1)$$

が得られ、$f_x^{-y+x} \equiv e_y^x$ と置けば、

$$f_{j+1}^{t+1} - f_j^t = -F(f_{j+1}^t - 1) + F(f_j^{t+1} - 1) \equiv -(\Delta_j - \Delta_t) F(f_j^t - 1) \tag{7}$$

となる。ここで $\Delta_j X_j^t \equiv X_{j+1}^t - X_j^t$、$\Delta_t X_j^t \equiv X_j^{t+1} - X_j^t$ である。

この方程式は $f_j^t$ の値を整数に制限すればフィルタ型CAを記述することに注意されたい。

---

## 式(7)と式(1)の同値性

式(1)はmin関数またはmax関数を用いて次のように表現できる：

$$u_j^{t+1} = \sum_{i=-\infty}^{j-1} u_i^t - \sum_{i=-\infty}^{j-1} u_i^{t+1} - \max\!\left(0,\; \sum_{i=-\infty}^{j-1} u_i^t - \sum_{i=-\infty}^{j-1} u_i^{t+1} + u_j^t - 1\right). \tag{8}$$

$S_j^t = \sum_{i=-\infty}^{j} u_i^t$ を導入すると、

$$S_{j+1}^{t+1} - S_j^t = -F(S_{j+1}^t - S_j^{t+1} - 1) \tag{9}$$

が得られる。これは $f_j^t = S_{j+1}^t - S_j^{t+1} = (\Delta_j - \Delta_t) S_j^t$ と置けば式(7)と同値である。こうして、CA (1)とKdV方程式およびLotka-Volterra方程式との関係が示された。

**図2.** KdV方程式 [式(2)] から可積分セルオートマトン [式(1)] への道筋。

```
KdV方程式 ... (2)
  連続極限 ↑↓ 離散化
Lotka-Volterra方程式 ... (3)
  連続極限 ↑↓ 離散化
式(3)の離散アナログ ... (4) or (5)
  ↓ （非解析的）極限操作
離散方程式（拡張されたCA）... (6)
  ↕ 変数変換
セルオートマトン ... (1)
```

---

## N-ソリトン解

上で明らかにしたように、式(1)は差分-差分方程式(4)の（非解析的）極限であるから、式(4)のN-ソリトン解が式(1)のN-ソリトン解になることが期待される。これは実際にその通りである。式(4)のN-ソリトン解は、広田の双線形恒等式を通じてCasorati行列式またはGram型行列式により与えられる [13]。$\varepsilon \to +0$ の極限をとれば、式(4)のN-ソリトン解が式(1)の解を次のように与えることが証明できる：

$$u_j^t = \rho_j^t - \rho_j^{t+1} - \rho_{j-1}^t + \rho_{j-1}^{t+1}, \tag{10}$$

ここで

$$\rho_j^t = \max_{\mu_i = 0,1}\!\left[\sum_{i=1}^N \mu_i \eta_i - \sum_{i > l} \mu_i \mu_l A_{il}\right], \tag{11}$$

$\eta_i = k_i j - \omega_i t + \eta_i^0$、$A_{il} = 2 \min(\omega_i, \omega_l)$、および

$$k_i = \text{sgn}(\omega_i)\, \min(1, |\omega_i|). \tag{12}$$

$\eta_i^0$（$i = 1, 2, \ldots, N$）は任意の整数、$\omega_i$（$i = 1, 2, \ldots, N$）は任意の（ただし全て正または全て負の）整数であり、$\max_{\mu_i=0,1}[X(\{\mu_i\})]$ は各 $\mu_i$ を0または1で置き換えて得られる $2^N$ 個の $X(\{\mu_i\})$ の値の最大値を表す。

例えば、3-ソリトン解は次のように与えられる：

$$\rho_j^t = \max[0,\, \eta_1,\, \eta_2,\, \eta_3,\, \eta_1 + \eta_2 - A_{12},\, \eta_2 + \eta_3 - A_{23},\, \eta_3 + \eta_1 - A_{31},\, \eta_1 + \eta_2 + \eta_3 - A_{12} - A_{23} - A_{31}].$$

$\rho_j^t$ は次の式を満たすことに注意すべきである：

$$\rho_{j+1}^{t+1} + \rho_j^{t-1} = \max\!\left[\rho_{j+1}^t + \rho_j^t,\; \rho_j^{t-1} + \rho_{j+1}^t - 1\right], \tag{13}$$

これは式(8)と(10)から得られ、双線形恒等式のアナログと考えることができる。

---

## 一般的な可積分CAの構成法

CA方程式(1)と可積分非線形波動方程式との関係が見出されれば、他の非線形可積分方程式から他の種類の可積分CAを構成することは直截的である。佐藤理論 [14] の観点から、可積分CAの一般的構成法を概説する。

1980年代初頭、佐藤幹夫はソリトンの統一理論を確立した。彼は、あらゆる可積分微分方程式が普遍Grassmann多様体（UGM）上の力学系と見なせることを示した。非線形方程式の解はUGMの点に対応し、τ関数と呼ばれる。UGMのPlücker関係式を用いると、τ関数に対する広田の双線形恒等式が得られる。伊達-神保-柏原-三輪は、場の理論と頂点作用素の方法により佐藤理論を無限次元Lie代数との関連に発展させた [15]。τ関数はフェルミオン場の演算子の真空期待値として表される。

通常のフェルミオン生成・消滅演算子を用いると、N-ソリトン解は次のように表される [16]：

$$\tau(\mathbf{t}) = \langle\text{vac}|\prod_{i=1}^{N}[1 + c_i \psi(p_i, \mathbf{t})\psi^*(q_i, \mathbf{t})]|\text{vac}\rangle \tag{14}$$

可積分CAを構成するために、$\mathbf{t} = \{t, j_1, j_2, \ldots, j_m\}$ と置き、条件

$$e^{\xi(\mathbf{t},p) - \xi(\mathbf{t},q)} = e^{(-\omega t + k_1 j_1 + \cdots + k_m j_m)/\varepsilon}$$

を課す。ここで $\omega$ と $k_j$ のうち2つは任意の整数である。この条件は分散関係 $\omega = \omega(k_1, \ldots, k_m : \varepsilon)$ を与える。同時に $p$ と $q$ は漸近形 $p = e^{P/\varepsilon} + \cdots$、$q = e^{Q/\varepsilon} + \cdots$ を持つ。$\rho(\mathbf{t}) \equiv \lim_{\varepsilon \to +0} \varepsilon \ln[\tau(\mathbf{t})]$ と置けば、$\rho(\mathbf{t})$ は式(13)に類似した方程式を満たし、ここから $m$ 次元CAとそのN-ソリトン解が得られる。CAはこのように構成されるため、τ関数の幾何学的・代数的性質を自然に継承する。

この方法により、一般化Lotka-Volterra系、1次元・2次元戸田格子、KP方程式に関連する他の種類の可積分CAも構成できる。解と保存量の詳細な解析とともに別の論文で報告する。

---

## 謝辞

広田良吾教授の有益なコメントに感謝する。著者の一人（時弘）はLotka-Volterra系に関する伊藤義教授の助言に感謝する。筧三郎、森脇慎一、永井敦、鳥居正城との議論にも謝意を表する。本研究は文部省科学研究費の一部により支援された。

---

## 参考文献

[1] 例えば、S. Wolfram 編、*Theory and Applications of Cellular Automata* (World Scientific, Singapore, 1986).

[2] K. Park, K. Steiglitz, W. P. Thurston, Physica D **19**, 423 (1986).

[3] T. S. Papatheodrou, M. J. Ablowitz, Y. G. Saridakis, Stud. Appl. Math. **79**, 173 (1988).

[4] T. S. Papatheodrou, A. S. Focas, Stud. Appl. Math. **80**, 165 (1989).

[5] M. Bruschi, P. M. Santini, O. Ragnisco, Phys. Lett. A **169**, 151 (1992); M. Bruschi, P. M. Santini, Physica D **70**, 185 (1994).

[6] D. Takahashi, J. Satsuma, J. Phys. Soc. Jpn. **59**, 3514 (1990).

[7] D. Takahashi, Proceedings of NOLTA '93 (1993), p. 555.

[8] M. Torii, D. Takahashi, J. Satsuma（出版予定）.

[9] 例えば、M. J. Ablowitz, H. Segur, *Solitons and the Inverse Scattering Transform* (SIAM, Philadelphia, 1985).

[10] 例えば、O. I. Bogoyavlenskii, Russ. Math. Surveys **46**:3, 1 (1991) およびその参考文献.

[11] Y. Itoh, Prog. Theor. Phys. **78**, 507 (1987).

[12] R. Hirota, S. Tsujimoto, J. Phys. Soc. Jpn. **64**, 3125 (1995).

[13] Y. Ohta, R. Hirota, S. Tsujimoto, T. Imai, J. Phys. Soc. Jpn. **62**, 1872 (1993).

[14] M. Sato, Y. Sato, in *Nonlinear PDEs in Applied Sciences*, eds. H. Fujita et al. (Kinokuniya/North-Holland, Tokyo, 1983); E. Date et al., Physica D **4**, 343 (1982); Y. Ohta et al., Prog. Theor. Phys. Suppl. **94**, 210 (1988).

[15] E. Date, M. Jimbo, M. Kashiwara, T. Miwa, Proceedings of RIMS Symposium (World Scientific, Singapore, 1983).

[16] E. Date, M. Jimbo, T. Miwa, J. Phys. Soc. Jpn. **51**, 4125 (1982).
