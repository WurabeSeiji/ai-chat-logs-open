# Paper 4 charged extension v1 — 実装・初期検証・課題

作成日: 2026-09-25

## 1. 目的

第四実験の状態機械・読み出し構造をできるだけ維持したまま、離散電荷

$$
n_A,n_B\in\mathbb Z
$$

と連続スケール

$$
\Gamma=\frac{k e_0^2}{G M^2}
$$

だけを追加した荷電拡張の最初の数値実装である。

相対運動の Newton 型中心結合は

$$
\alpha=1-\frac{\Gamma n_A n_B}{\nu},\qquad
\nu=\frac{m_A m_B}{(m_A+m_B)^2}
$$

とし、

$$
\left|\frac{F_C}{F_G}\right|
=\frac{\Gamma |n_A n_B|}{\nu}
$$

を制御量とした。

第四実験の

$$
(\chi,p,e,\phi,t)\rightarrow X\rightarrow Y
$$

および

$$
S_n=B_{n+1}B_n^{-1},\qquad \det S_n=1
$$

はそのまま維持した。

## 2. 追加した読み出し

### 2.1 GW

質量四重極の二階微分から face-on の正規化波形

$$
h_+,\ h_\times
$$

を生成する。波形係数は

$$
\frac{\nu}{D}
$$

なので、既知または共通の距離スケールを除けば質量側の交換不変量 $\nu$ を読む。

### 2.2 EM

電気双極子係数

$$
\beta_q=n_A x_B-n_B x_A,
\qquad x_A=\frac{m_A}M,\quad x_B=\frac{m_B}M
$$

を用い、正規化遠方電場を

$$
(E_x,E_y)
\propto
\frac{\sqrt\Gamma\,\beta_q}D\,\ddot{\mathbf r}
$$

として生成する。

したがって EM 波形から

$$
\Gamma\beta_q^2
$$

を匿名量として読める。

### 2.3 軌道・時計から Coulomb 積を読む

中心結合から

$$
\Gamma n_A n_B=\nu(1-\alpha)
$$

である。

一周周期から $\alpha$ を読む方法と、第四実験の局所時計

$$
\frac{dt}{d\chi}=\frac{r^2}{\sqrt{\alpha p}}
$$

から局所的に

$$
\boxed{
\alpha
=\frac{r^4}{p(dt/d\chi)^2}
}
$$

を読む方法を両方実装した。

## 3. 放射反作用

v1 では、荷電した軌道でも使えるよう、Paper 4 の Kidder 2.5PN 局所項をそのまま流用するのではなく、Newton 型 osculating orbit 上で先頭多重極フラックスを直接評価した。

- GW: 質量四重極放射
- EM: 電気双極子放射

をそれぞれエネルギー・角運動量流束へ変換し、

$$
p=\frac{h^2}\alpha,
\qquad
E_{orb}=-\frac{\alpha(1-e^2)}{2p}
$$

から $\dot p,\dot e$ を得ている。

固定軌道での解析式との照合では、数値平均と解析平均は機械精度で一致した。

- uncharged GW: Peters の平均四重極パワーと相対差 $\sim2.2\times10^{-16}$
- charged GW, EM: 解析平均式と相対差 $\sim2\text{--}3\times10^{-16}$

## 4. 推奨した最初の試験条件

離散電荷は

$$
(n_A,n_B)=(+1,-1)
$$

質量比は

$$
m_A/m_B=1,\qquad \nu=1/4
$$

とした。

Coulomb/重力比を

$$
\left|F_C/F_G\right|=100
$$

とすると

$$
\Gamma=25,
\qquad
\alpha=101.
$$

同じ Paper-4 C1 の $p=24.3$ をそのまま使うと速度が超光速になるため、軌道スケールをまず $\alpha$ 倍し、さらに放射反作用を十分緩やかにするため 10 倍した。

したがって

$$
\text{scale}=1010,
\qquad
p_0=24543,
\qquad
e_0=0.435889894354.
$$

## 5. 実測結果 — 推奨ケース

3 radial periods、4000 steps/orbit、Newtonian angular mode で実行した。

|量|実測|
|---|---:|
|$\max(v/c)$|0.092278372|
|$p_{final}/p_0$|0.993419255|
|$e_{final}/e_0$|0.995059824|
|初周平均 EM/GW 放射パワー比|6960.575819|
|$\max|\det S_n-1|$|5.551e-16|
|state-map error max|5.745e-12|

放射反作用を含めても 3周で

$$
\Delta p/p_0\approx-0.6581\%,
$$

$$
\Delta e/e_0\approx-0.4940\%
$$

に抑えられた。

## 6. 読み出し結果

### 質量側

GW 波形係数から

$$
\nu_{fit}=0.250000000000
$$

を得た。真値は $0.25$。

### 電荷双極子側

EM 波形から

$$
(\Gamma\beta_q^2)_{fit}
=25.000000000000
$$

真値は

$$
\Gamma\beta_q^2=25.000000000000.
$$

### Coulomb 積側

一周周期だけで読むと

$$
\alpha_{period}=101.449722254
$$

で、真値 $101$ から放射反作用の影響を受ける。

一方、局所時計からは

$$
\boxed{\alpha_{local}=100.999939562}
$$

となり、

$$
(\Gamma n_A n_B)_{local}
=-24.999984891
$$

で、真値 $-25$ をほぼそのまま回収した。

既知の $\Gamma=25$ を用いて整数候補を列挙すると最良候補は

$$
(n_A,n_B)=(-1,+1),\ (+1,-1)
$$

であり、交換以外は排除された。

## 7. 数値収束

推奨条件を1周だけ、2000 / 4000 / 8000 steps/orbit で比較した。

8000 step を参照した 4000 step の差は

- $|\Delta p|=2.547e-11$
- $|\Delta e|=1.776e-15$
- $|\Delta t|=2.887e-08$

で、今回見えている secular 変化は積分刻みのアーティファクトではない。

# 8. 抽出された問題・課題

## 課題 A — C1 の数値スケールをそのまま使えない

Coulomb 比を100にして $p=24.3$ のままだと、初期近点速度は概算で

$$
v_p\sim\sqrt{101}\,v_{p,C1}>c
$$

となる。

したがって「第四実験へ電荷だけ追加」は**状態機械としては可能**だが、同じ数値半径をそのまま使うことはできない。

今回は $p,r$ を $\alpha$ 倍、さらに10倍して低速化した。

## 課題 B — EM dipole radiation が非常に強い

軌道を $\alpha$ 倍するだけのケースでは、3周で

$$
p/p_0=0.779026,
\qquad
e/e_0=0.829680
$$

まで変化した。

この場合、一周周期から Coulomb 結合を読むと

$$
\alpha_{period}=116.485
$$

となり、真値101から大きくずれる。

これは数値誤差ではなく、v1 の leading EM dipole RR が実際に強いことによる。

**対策:** 軌道をさらに大きくして adiabatic 化するか、局所時計読み出しを使う。

## 課題 C — 「余分な状態」が読み出しに意味を持つことが実際に現れた

今回、全周周期のような粗い読み出しは放射反作用でバイアスされたが、$t$ と $p,r,\chi$ を使う局所時計読み出しは

$$
\alpha_{local}\simeq101
$$

を回収した。

したがって、力学核だけを見ると冗長に見える時計状態 $Q=e^{t/\tau_0}$ が、**電荷結合と質量結合を匿名に区別する局所計量として有効になる**可能性が、具体的な数値例で出た。

これは現在の状態削減議論に対して重要な反例候補である。

## 課題 D — $\Gamma$ が未知なら電荷量子数とスケールが縮退する

等質量で $(+1,-1)$ の場合、

$$
\Gamma n_A n_B=-\Gamma,
\qquad
\Gamma\beta_q^2=\Gamma.
$$

しかし $\Gamma$ 自身を未知にすると、例えば

$$
(n_A,n_B)=(+1,-1),\ \Gamma=25
$$

と

$$
(n_A,n_B)=(+2,-2),\ \Gamma=25/4
$$

は同じ二つの不変量を作れる。

したがって絶対的な整数電荷を読むには、

- $e_0$ を普遍単位として既知とする、または
- 全質量 $M$ を独立に読む、または
- 追加の高次 EM 読み出しを使う

必要がある。

## 課題 E — $\beta_q=0$ セクター

$$
\frac{q_A}{m_A}=\frac{q_B}{m_B}
$$

なら electric dipole radiation は消える。

この場合、Coulomb 相互作用は存在しても今回の EM dipole readout はゼロになるため、電気四重極・磁気双極子などの高次放射を追加しないと一般の電荷セクターを読めない。

## 課題 F — conservative charged GR はまだ未導出

`paper4` angular mode は既存の Schwarzschild/Darwin 角度因子をそのまま残せるが、Coulomb 優勢時には Einstein-Maxwell/Lorentz-force から導出したものではない。

したがって今回の推奨ケースは conservative 部を `newtonian` とし、まず Coulomb + GW/EM readout の構造だけを検証した。

次段階では

- 弱重力 + Coulomb の摂動論
- charged test-particle in curved spacetime
- 必要なら Einstein-Maxwell 二体系

のどのレベルまで採用するかを分離する必要がある。

## 課題 G — GW の EM-field stress-energy を無視している

v1 の GW は質点質量四重極だけで生成している。電磁場自身の stress-energy も重力源なので、精密な Einstein-Maxwell 波形では追加寄与がある。

ただし推奨ケースの近点で Coulomb potential energy / $Mc^2$ は概算

$$
\frac{\Gamma}{r_p}\sim1.5\times10^{-3}
$$

なので、最初の構造試験として質点四重極を使う余地はある。

## 課題 H — 実在の素電荷 $e$ と直接同一視すると古典実験としては極端

$\Gamma=25$ を本物の elementary charge $e$ と結びつけると、全質量は約

$$
M\approx3.72\times10^{-10}\,\mathrm{kg}
$$

となる。今回の dimensionless 軌道をそのまま物理長へ写すと古典点粒子モデルとして不自然な極小長さになる。

したがって現段階では

$$
\boxed{\text{構造検証用の無次元思考実験}}
$$

として扱い、具体的な実在粒子への写像は別課題とするべきである。

# 9. 現時点の最も重要な結果

第一の荷電拡張だけで、第四実験の状態構造をほぼそのまま用いながら、

$$
\boxed{\nu}
$$

を GW から、

$$
\boxed{\Gamma\beta_q^2}
$$

を EM 波から、

$$
\boxed{\Gamma n_A n_B}
$$

を軌道 + 局所時計から、別々に読み出せた。

したがって、質量と電荷の区別を読むためには「軌道を作る最小状態」だけではなく、局所時計や位相などの**読み出し計量を増やす状態**が意味を持つ、という今回の仮説を実際に検査できる形になった。
