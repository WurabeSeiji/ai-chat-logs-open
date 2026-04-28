# 論文4：物理的接続 — $\Delta(R)$ と Bekenstein–Hawking エントロピーの対応 構成草案

**仮タイトル**: A Quantitative Correspondence between the Lattice Volume Deficit $\Delta(R)$ and the Tangherlini Bekenstein–Hawking Entropy

**分量目安**: 10–15 ページ

---

## 目次

| § | タイトル | 分量 | 状態 |
|---|---|---|---|
| §1 | Introduction | 1.5 ページ | 未着手 |
| §2 | Volume Deficit as Discrete-Drift Entropy | 2 ページ | 草稿要点 |
| §3 | Tangherlini Bekenstein–Hawking in 4+1 dimensions | 2 ページ | 草稿要点 |
| §4 | Quantitative Correspondence | 3 ページ | 草稿要点 |
| §5 | Interpretation of the Coefficient Ratio | 2 ページ | 草稿要点 |
| §6 | Discussion and Open Problems | 2 ページ | 未着手 |

---

## §2 $\Delta(R)$ の離散ドリフトエントロピーとしての解釈

論文3 で確立された $\Delta(R) = V_4(R) - N(k) \sim (16\pi/3) R^3$ を、Planck 単位の自由度として解釈する：

$$\frac{\Delta(R)}{\ell_P^4} = (\text{4次元体積に対応する Planck 単位数のうち離散単位で詰められない部分})$$

これを、4次元球（上位時空構造の有限領域）から「ドリフトする」自由度の数として解釈する。

### 仮定

- 中心投影由来の4次元主観空間（論文2 由来）における離散最小単位は4次元 Planck 体積 $\ell_P^4$ を持つ
- $\Delta(R)$ は単位立方体充填では捉えられない「境界自由度」を表す

### 本論文の主張

$\Delta(R)/\ell_P^4$ が、対応する4+1次元 Tangherlini ブラックホールの Bekenstein–Hawking エントロピーと、$R^3$ スケーリングおよび先行係数の両方で対応する。

---

## §3 4+1 次元 Tangherlini Bekenstein–Hawking エントロピー

### §3.1 Tangherlini 解（再掲・簡略化）

4+1次元 Schwarzschild–Tangherlini 解：
$$ds^2 = -f_T(r) dt^2 + f_T(r)^{-1} dr^2 + r^2 d\Omega_3^2, \qquad f_T(r) = 1 - \frac{\mu}{r^2}$$
where $\mu = 8GM/(3\pi)$ (5次元 Newton 定数 $G_5$ を使用、$d=4$ 次元空間想定で標準正規化).

地平面：$r_h = \sqrt{\mu} = \sqrt{8GM/(3\pi)}$

### §3.2 Bekenstein–Hawking エントロピー

地平面の3次元体積 $A_3 = 2\pi^2 r_h^3$（$S^3$ の体積、半径 $r_h$）。

Bekenstein–Hawking エントロピー：
$$S_{BH} = \frac{A_3}{4 G_5 \hbar} = \frac{2\pi^2 r_h^3}{4 G_5} = \frac{\pi^2 r_h^3}{2 G_5}.$$

（Planck 単位 $\ell_P^3 = G_5 \hbar/c^3$ 等の正規化は **後で検討**）

スケーリング：$S_{BH} \propto r_h^3 \propto R^3$（$R$ を $r_h$ と同一視するとき）。

---

## §4 定量的対応

### §4.1 $R^3$ スケーリングの一致

論文3 の漸近結果：$\Delta(R) \sim (16\pi/3) R^3$
4+1 次元 BH エントロピー：$S_{BH} \propto r_h^3$

両者の $R^3$ 依存性は完全一致。

### §4.2 先行係数の比較

$R = r_h$ と同定すると（中心投影による主観空間の曲率半径 $R$ を 5次元 BH の地平面半径 $r_h$ に対応させる仮定）：

- $\Delta(R) / R^3 \to 16\pi/3 \approx 16.755$
- $S_{BH} / r_h^3 = \pi^2/2 \approx 4.935$（Planck 単位の選択依存）

比：
$$\frac{\Delta(R)}{S_{BH}} \to \frac{16\pi/3}{\pi^2/2} = \frac{32}{3\pi} \approx 3.397.$$

### §4.3 数値（有限 $k$）

$k=10, R=21$: $\Delta = 141580$, $S_{BH} = (\pi^2/2) R^3 = 45712$. 比 = $3.097$.
$k=60, R=121$: $\Delta = 29303165$, $S_{BH} = (\pi^2/2)(121)^3 = 8740869$. 比 = $3.353$.

漸近比 $32/(3\pi) \approx 3.40$ への収束が確認される。

---

## §5 係数比 $32/(3\pi)$ の解釈

### §5.1 数学的構造

$$\frac{\Delta}{S_{BH}} = \frac{16\pi/3}{\pi^2/2} = \frac{32}{3\pi}$$

分子：論文3 の包除原理で導出された境界補正 $16\pi/3$。
分母：4次元球の表面（3次元体積）/ 体積比から派生する $\pi^2/2$。

### §5.2 物理的解釈（仮説的）

比 $32/(3\pi)$ は**離散単位の境界格子分率**に対応する可能性がある。すなわち：

- $S_{BH}$ は地平面（3次元表面）の自由度
- $\Delta(R)$ は境界スラブ厚 $1/2$ の中の体積
- 比 $32/(3\pi) \approx 3.40$ は「境界スラブ／表面薄殻」の幾何学的比

詳細な解釈は **後で検討**。

### §5.3 残課題

- $G_5$, $\ell_P$ の正規化と $R$ の同定の精密化
- 比 $32/(3\pi)$ が物理的に意味を持つか、単なる係数の偶然か
- 非漸近的補正（副主要項 $-6\pi R^2$, etc.）の物理的意味

---

## §6 議論

- 論文3 から論文4 への「数学から物理へ」の移行
- 「$\Delta(R)$ は離散ドリフトエントロピー」という解釈の限界
- 完全な定量一致を主張せず、$R^3$ スケーリング一致を主要結果とする

---

**自己引用なし** / **論文3 への引用は本論文公表後に許容（執筆段階では参照しない）**

**未確定事項（後で検討）**:
- §3 の Tangherlini 解の正規化の選択（$G_5$, $\ell_P$, $r_s$ の関係）
- §4 の $R = r_h$ 同定の正当化（論文5 で動的に導出）
- §5 の解釈の選択肢
