# 曲率半径の自己整合的決定と Tangherlini 蒸発スケーリングの回復（論文5：動力学的接続）

**著者**：木原 範昭（WF System Co., Ltd.）
**ORCID**：[0009-0004-6753-4020](https://orcid.org/0009-0004-6753-4020)
**日付**：2026年4月
**DOI**：[10.5281/zenodo.19837596](https://doi.org/10.5281/zenodo.19837596)

---

## 要旨

中心投影構成は4次元主観座標図の自然な幾何学的スケールとして曲率半径 $R$ を導入する。全エネルギー $E_{\mathrm{total}}$ を持つ物質の存在下で、$R$ が Friedmann 様関係 $E_{\mathrm{total}} \propto R^2$ から自己整合的に決定され、4+1次元 Schwarzschild–Tangherlini ブラックホール特有の質量–半径スケーリング $R \propto M^{1/2}$ に至ることを論じる。地平面 $r = R$ での表面重力から計算される Hawking 温度は $T_H \propto 1/\sqrt{M}$ で、Stefan–Boltzmann 様の蒸発論証は寿命 $\propto M^{3/2}$ を与え、Tangherlini 結果と一致する。当初プログラム要旨で指摘された蒸発レート指数の「不一致」は、3+1 と 4+1 規約の混同であったことを §5.1 で明らかにする。

---

## §1. 序論

論文2 の中心投影構成は、曲率半径 $R$ を自由幾何学的パラメータとして含む4次元計量を生む。この構成をブラックホール熱力学に接続するには、4次元理論の物質内容（特にブラックホールの質量 $M$）から $R$ を決定する必要がある。

本論文は $M$ から $R$ を固定する自己整合論証を提示する。論証は Friedmann 方程式の構造を持つ：曲がった幾何学に関連する全エネルギーが物質の総質量エネルギーと同定され、結果の方程式が $R$ を決定する。

主要結果：

1. **質量–半径関係**：$R \propto M^{1/2}$。
2. **Hawking 温度**：$T_H \propto 1/\sqrt{M}$。
3. **寿命**：$\tau \propto M^{3/2}$（4+1次元）。
4. **蒸発レート指数**：$dM/dt \propto -M^{n}$ で $n = -1/2$（4+1次元、Tangherlini と一致）。

最初の3つは Tangherlini 4+1次元 Schwarzschild 結果と定量的に一致する。プログラム当初要旨で言及された「不一致」は3+1 と 4+1 規約の混同であり、4+1次元では一致する。

論文構成。第2節で自己整合関係 $E_{\mathrm{total}} \propto R^2$ を定式化。第3節で質量–半径スケーリングを導出。第4節で Hawking 温度と寿命を計算。第5節で蒸発レート指数の解析と当初の混同の解明。第6節で結論。

---

## §2. 全エネルギーと曲率半径

### §2.1 Friedmann 類似

論文2 の中心投影計量は（質量パラメータを伴う静的球座標で）形
$$ds^2 = -f(r)\, dt^2 + f(r)^{-1} dr^2 + r^2\, d\Omega_2^2, \qquad f(r) = 1 - \frac{2GM}{r} - \frac{r^2}{R^2}$$
を持つ。これは宇宙項 $\Lambda = 3/R^2$ を持つ Schwarzschild–de Sitter 計量。

宇宙項 $\Lambda$ を伴う閉宇宙論に適用された4次元 Friedmann 方程式：
$$H^2 + \frac{k}{a^2} = \frac{8\pi G}{3} \rho + \frac{\Lambda}{3}$$
ここで $H = \dot a/a$ は Hubble パラメータ、$k$ は空間曲率、$\rho$ はエネルギー密度、$a$ はスケール因子。

宇宙項 $\Lambda = 3/R^2$ を持つ閉（正曲率）宇宙論で、サイズ $R$ の領域内の物質全エネルギーと曲率の関係は
$$E_{\mathrm{total}} \sim \frac{R^2}{G}.$$
全エネルギーは $R^2$ にスケールし、比例定数は Newton 定数で設定される。

### §2.2 $E_{\mathrm{total}}$ と $M c^2$ の同定

4次元理論の質量 $M$ のブラックホールで、全エネルギーは
$$E_{\mathrm{total}} = M c^2.$$

§2.1 と組み合わせ：
$$M c^2 \propto \frac{R^2}{G} \quad \Longrightarrow \quad R^2 \propto G M / c^2.$$

自然単位（$c = 1$, $G = 1$）で $R^2 \propto M$、すなわち
$$\boxed{\;R \propto M^{1/2}.\;}$$

### §2.3 Tangherlini 地平面半径との比較

4+1次元 Schwarzschild–Tangherlini 地平面：
$$r_h = \left(\frac{8 G_5 M}{3\pi}\right)^{1/2}.$$

スケーリング $r_h \propto M^{1/2}$ は中心投影スケーリング $R \propto M^{1/2}$ と同一。これは論文4 でエントロピー比 $\Delta(R)/S_{BH}$ を計算するために用いられた同定 $R = r_h$ を正当化する。

比例定数は4次元 $G$ に対する $G_5$ の規格化に依存する因子で Tangherlini 値と異なる。**後で検討**：定数の正確な調整。

---

## §3. 質量–半径関係

### §3.1 関係の詳細

§2.2 から、自然単位で $R = c_1 M^{1/2}$（無次元定数 $c_1$）。

姉妹論文（論文4）はこの $R$ を Tangherlini 地平面半径 $r_h$ と同定する。Tangherlini 関係 $r_h^2 = 8 G_5 M/(3\pi)$ により
$$R^2 = \frac{8 G_5 M}{3\pi}$$
で、その単位で $c_1^2 = 8 G_5/(3\pi)$ が固定される。

### §3.2 妥当性領域

関係 $R \propto M^{1/2}$ は以下の領域で成立：
- 質量が Planck スケールに比して大（半古典解析が妥当）
- 質量が任意の宇宙論的スケールに比して小（局所近似が妥当）
- 物質が $R$ に比して本質的に点状

標準的天体物理または理論的領域のブラックホールでは3条件すべて満足。

---

## §4. Hawking 温度と寿命

### §4.1 表面重力と温度

中心投影計量の地平面（支配的地平面で $f(r) = 0$）は $r = R$（質量項が宇宙項に比して小さい極限で）。この地平面での表面重力：
$$\kappa = \frac{1}{2} |f'(r)|_{r=R} = \frac{R}{R^2} = \frac{1}{R}.$$

（より精密に、事象地平面 $r_h \sim R$ での Tangherlini 類似は $\kappa = 1/r_h$。）

Hawking 温度：
$$T_H = \frac{\kappa}{2\pi} = \frac{1}{2\pi R} \propto \frac{1}{R} \propto M^{-1/2}.$$

これは $r_h \leftrightarrow R$ で Tangherlini 結果 $T_H = 1/(2\pi r_h)$ と一致。

### §4.2 Stefan–Boltzmann による寿命

Stefan–Boltzmann 様の輻射光度：
$$\frac{dE}{dt} \sim \sigma A_3 T_H^4$$
ここで $A_3 = 2\pi^2 R^3$ は3球地平面の面積、$T_H \sim 1/R$、$\sigma$ は Stefan-Boltzmann 様定数。

代入：
$$\frac{dM}{dt} \sim -\sigma \cdot R^3 \cdot \frac{1}{R^4} = -\frac{\sigma}{R}.$$

$R \propto M^{1/2}$ により：
$$\frac{dM}{dt} \propto -M^{-1/2}.$$

積分：
$$\int M^{1/2}\, dM \propto -t \quad \Longrightarrow \quad M^{3/2} \propto -t \quad \Longrightarrow \quad \tau \propto M^{3/2}.$$

これは4+1次元での寿命。

### §4.3 Tangherlini スケーリング法則の総合

$D = 5$（4+1次元）で：
- $r_h \propto M^{1/2}$
- $T_H \propto r_h^{-1} \propto M^{-1/2}$
- $A_3 = 2\pi^2 r_h^3 \propto M^{3/2}$
- $S_{BH} \propto A_3 \propto M^{3/2}$
- $dM/dt \propto -M^{-1/2}$
- $\tau \propto M^{3/2}$

これら全てが $R = r_h$ で中心投影構成により回復される。

---

## §5. 蒸発レート指数：当初の混同の解明と真の未解決問題

### §5.1 当初要旨の「不一致」の正体

プログラム要旨で言及された「指数 $n=3/2$ vs $n=-1$ の不一致」は、3+1 と 4+1 次元の混同であった。

- **3+1次元 Hawking**：$dM/dt \propto -A_2 T_H^4 \propto -M^2 \cdot M^{-4} = -M^{-2}$、寿命 $\tau \propto M^3$。
- **4+1次元 Tangherlini**：$dM/dt \propto -A_3 T_H^4 \propto -M^{3/2} \cdot M^{-2} = -M^{-1/2}$、寿命 $\tau \propto M^{3/2}$。

中心投影構成（4+1次元）は4+1次元 Tangherlini 結果と一致し、3+1次元 Hawking 結果とは一致しない（自然な比較対象が4+1だから）。

### §5.2 真の未解決問題：4+1 と 3+1 の関係

中心投影構成が物理的 3+1 ブラックホールの記述であろうとするなら、4+1 蒸発が射影された理論で 3+1 蒸発として現れる機構が必要。可能性：

1. *次元削減*。4空間次元の1つが $R$ より遥か小さいスケールでコンパクト（Kaluza–Klein 様）なら、有効 3+1 理論は異なる蒸発レートを持つ。
2. *ブレーン局在化*。物質が 3+1 次元ブレーンに局在し重力が 4+1 で伝播すれば、コンパクト化スケール逆数より高エネルギーの Hawking 量子はバルクへ脱出し4次元蒸発レートを修正。
3. *ホログラフィー的射影*。4次元理論を 3+1 理論にホログラフィー様の射影で関係づけるなら、蒸発レートは射影で修正される。

現段階でいずれの解決も特権化しない。4+1（$\tau \propto M^{3/2}$）と 3+1（$\tau \propto M^3$）の寿命の差異は、中心投影構成が観測物理学への接触を主張するなら扱うべき真の物理的問題。

### §5.3 プログラムが主張しないこと

プログラムは 4+1 寿命 $\tau \propto M^{3/2}$ が物理的ブラックホールの観測寿命であるとは主張しない。4+1 次元と中心投影構成で寿命が Tangherlini 結果と整合的にスケールすると主張するのみ。3+1 物理学との関係は追加の問いで、未解決と明示。

---

## §6. 結論

中心投影構成の曲率半径 $R$ が物質内容から $E_{\mathrm{total}} \propto R^2$ により自己整合的に決定され、Tangherlini スケーリング法則
- $R \propto M^{1/2}$
- $T_H \propto 1/\sqrt{M}$
- $\tau \propto M^{3/2}$
- $dM/dt \propto -M^{-1/2}$

を導くことを示した。4つすべてが4+1次元 Schwarzschild–Tangherlini ブラックホールと定量的に一致。物理的 3+1 次元ブラックホール蒸発との関係は未解決問題として残す。

本論文の主要寄与は $M$ から $R$ を固定する自己整合論証で、論文4 で用いた同定 $R = r_h$ を正当化する。構成は内的に整合し 4+1 次元で Tangherlini 結果を再現；3+1 次元への更なる削減は将来研究の対象。

---

## 参考文献（選定、**後で検討**）

- F. R. Tangherlini, *Nuovo Cimento* 27, 636 (1963).
- R. C. Myers and M. J. Perry, *Annals Phys.* 172, 304 (1986).
- L. Randall and R. Sundrum, *Phys. Rev. Lett.* 83, 3370 (1999).
- N. Arkani-Hamed, S. Dimopoulos, and G. Dvali, *Phys. Lett. B* 429, 263 (1998).
