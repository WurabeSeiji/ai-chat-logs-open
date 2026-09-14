# 匿名 Bob-star 完全不変量フレームと複素ガウス null 幾何における κ の厳密分布

**Complete Invariant Frames for Anonymous Bob-Anchored Relation Stars and the Exact Distribution of κ in Complex Gaussian Null Geometry**

**著者:** 木原範昭 / Noriaki Kihara  
**ORCID:** 0009-0004-6753-4020  
**DOI:** （取得後記入）  
**Concept DOI:** （取得後記入）  
**系列:** 匿名頂点状態生成幾何  
**位置づけ:** 論文2 / 測定装置較正・null 幾何基準  
**版:** v1.1 / 2026-09-14

---

## 要旨

本稿は、匿名頂点状態生成幾何において、内部観測者 Bob が自分を端点とする無順序複素関係集合（Bob-star）から読み出せる完全不変量フレームを構成し、その数値較正と、generic な複素ガウス null 集団における二乗閉包秩序量

$$
\kappa
=
\frac{\left|\sum_j z_j^2\right|}
{\sum_j |z_j|^2}
$$

の厳密分布を導出する。

Bob-star の関係数を

$$
m=N-1
$$

とし、内部ゲージ群を

$$
G_B
=
U(1)\times\mathbb R_+\times S_m
$$

とする。正規化冪和

$$
P_k
=
\frac{\sum_j z_j^k}
{\left(\sum_j|z_j|^2\right)^{k/2}}
$$

を用いると、\(P_1\neq0\) の chart において

$$
\Phi_B^{(K)}
=
\{|P_1|\}
\cup
\left\{
|P_k|,
\arg\left(P_k\overline{P_1^k}\right)
\right\}_{k=2}^{K}
$$

は \(G_B\) 不変である。さらに \(K=m=N-1\) まで用いれば、Newton の恒等式により基本対称式を復元し、そのモニック多項式の根として無順序 Bob-star を完全復元できる。したがって

$$
\boxed{
\Phi_B^{(N-1)}
\text{ は }P_1\neq0\text{ chart 上の完全座標系}
}
$$

であり、

$$
\boxed{
K_{\rm complete}\le N-1
}
$$

が成立する。

低位 \(N=2,\ldots,6\) に対し、各 2000 標本の独立複素ガウス乱数を用いて、この不変性・完全復元性・N=3 恒等式群を数値較正した。最大不変誤差は約 \(3.0\times10^{-14}\)、最大復元誤差は約 \(9.9\times10^{-14}\) であった。

さらに、\(z_j=x_j+iy_j\) から作る \(2\times m\) 実ガウス行列の Gram 行列の固有値を \(\lambda_1\ge\lambda_2\) とすると、

$$
\boxed{
\kappa
=
\frac{\lambda_1-\lambda_2}
{\lambda_1+\lambda_2}
}
$$

が恒等的に成立する。これを 2×m Wishart 固有値分布に適用すると、

$$
\boxed{
\kappa^2
\sim
\operatorname{Beta}
\left(
1,\frac{m-1}{2}
\right)
}
$$

を得る。

従って、

$$
f_\kappa(q)
=
(m-1)q(1-q^2)^{(m-3)/2},
$$

$$
F_\kappa(q)
=
1-(1-q^2)^{(m-1)/2},
$$

$$
E[\kappa]
=
\frac{\sqrt\pi}{2}
\frac{\Gamma((m+1)/2)}
{\Gamma((m+2)/2)},
$$

$$
\operatorname{median}(\kappa)
=
\sqrt{
1-2^{-2/(m-1)}
}
$$

が厳密に得られる。

特に \(m=2\) では、

$$
\boxed{
E[\kappa]=\frac{\pi}{4},
\qquad
\operatorname{median}(\kappa)=\frac{\sqrt3}{2}
}
$$

であり、大 \(m\) 極限では

$$
\boxed{
\sqrt m\,\kappa
\Rightarrow
\operatorname{Rayleigh}(1),
\qquad
E[\kappa]\sim\sqrt{\frac{\pi}{2m}}
}
$$

となる。

さらに、

$$
P(\kappa\le\epsilon)
\sim
\frac{m-1}{2}\epsilon^2
$$

より、厳密二乗ゼロ閉包 \(\kappa=0\) は generic 複素ガウス null 集団において測度ゼロ境界である。したがって、自己無撞着な関係波力学が \(\kappa=0\) またはその近傍を選択する場合、それは generic null から最大感度で識別可能である。

本稿は、今後の終端状態族の内部次元測定

$$
\dim(\mathcal M_{\rm terminal}(N)/G_B)
$$

に用いる測定装置と null 基準曲線を、解析・数値の両面から較正するものである。

---

# 1. 問題設定

Bob-star を

$$
\mathcal Z_B
=
\{z_1,\ldots,z_m\},
\qquad
m=N-1
$$

とする。

Bob から見て関係の順序は物理属性ではない。さらに共通位相と共通正スケールも内部ゲージとする。

従って、

$$
\boxed{
G_B
=
U(1)\times\mathbb R_+\times S_m
}
$$

を内部ゲージ群とする。

---

# 2. 正規化冪和フレーム

## 2.1 定義

$$
H
=
\sum_{j=1}^m |z_j|^2
$$

とする。

正規化冪和を

$$
\boxed{
P_k
=
\frac{\sum_jz_j^k}
{H^{k/2}}
}
$$

と定義する。

大域位相

$$
z_j\mapsto e^{i\theta}z_j
$$

では、

$$
P_k
\mapsto
e^{ik\theta}P_k.
$$

従って \(P_k\) は U(1) 電荷 \(k\) を持ち、置換 \(S_m\) には不変である。

## 2.2 完全不変量フレーム

\(P_1\neq0\) の chart では、

$$
|P_1|
$$

および

$$
|P_k|,
\qquad
\arg\left(
P_k\overline{P_1^k}
\right)
$$

は \(G_B\) 完全不変である。

従って、

$$
\boxed{
\Phi_B^{(K)}
=
\{|P_1|\}
\cup
\left\{
|P_k|,
\arg\left(P_k\overline{P_1^k}\right)
\right\}_{k=2}^{K}
}
$$

を Bob-star 不変量フレームとする。

---

# 3. 完全性定理

Bob-star の要素数を \(m\) とする。

\(P_1\neq0\) の chart で、

$$
H=1
$$

かつ

$$
\arg P_1=0
$$

となるようにゲージ固定する。

すると、

$$
p_1=|P_1|
$$

および

$$
p_k
=
|P_k|
\exp
\left[
i\arg\left(P_k\overline{P_1^k}\right)
\right]
$$

として、

$$
p_k
=
\sum_{j=1}^m z_j^k
$$

を \(k=1,\ldots,m\) まで復元できる。

Newton の恒等式により、

$$
p_1,\ldots,p_m
$$

から基本対称多項式

$$
e_1,\ldots,e_m
$$

を逐次復元できる。

そして無順序 star \(\{z_j\}\) は、

$$
\boxed{
\lambda^m
-e_1\lambda^{m-1}
+e_2\lambda^{m-2}
-\cdots
+(-1)^m e_m
=0
}
$$

の根集合である。

従って、

$$
\boxed{
K=m=N-1
}
$$

で完全復元でき、

$$
\boxed{
\Phi_B^{(N-1)}
\text{ は完全座標系}
}
$$

が成立する。

---

# 4. 低位 N 数値較正

## 4.1 条件

- \(N=2,3,4,5,6\)
- 各 N 2000 標本
- \(\Re z_j,\Im z_j\sim\mathcal N(0,1)\)
- seed = 20260914
- ランダム共通スケール
- ランダム大域位相
- ランダム \(S_m\) 置換

## 4.2 結果

| N | m | 標本数 | 最大不変誤差 | 最大復元誤差 | 平均 κ | 中央値 κ |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 1 | 2000 | 5.551e-16 | 0.000e+00 | 1.000000 | 1.000000 |
| 3 | 2 | 2000 | 9.770e-15 | 5.970e-15 | 0.784638 | 0.866433 |
| 4 | 3 | 2000 | 3.042e-14 | 4.073e-14 | 0.664529 | 0.709719 |
| 5 | 4 | 2000 | 8.882e-15 | 9.020e-14 | 0.579951 | 0.599575 |
| 6 | 5 | 2000 | 9.326e-15 | 9.849e-14 | 0.530283 | 0.538024 |

最大不変誤差は \(3.04\times10^{-14}\) 以下、最大復元誤差は \(9.85\times10^{-14}\) 以下であった。

独立実装・別 seed による対照検算でも、同じ精度水準で不変性・復元性が再現された。

---

# 5. N=3 の閉形式恒等式

\(m=2\) では、

$$
r
=
\frac{|z_1|}{|z_2|},
\qquad
\Delta
=
\arg(z_1/z_2)
$$

とし、

$$
u=r+\frac1r,
\qquad
c=\cos\Delta,
\qquad
w=
\left(r-\frac1r\right)\sin\Delta
$$

を定義する。

すると、

$$
\boxed{
w^2=(u^2-4)(1-c^2)
}
$$

である。

さらに、

$$
S=
\frac{z_1}{z_2}
+
\frac{z_2}{z_1}
$$

とすると、

$$
\boxed{
S=uc+iw
}
$$

となる。

論文1 [KG1] で導入した局所複素座標

$$
J_2
=
\frac{z_1z_2}{(z_1+z_2)^2}
$$

とは、

$$
\boxed{
\frac1{J_2}=2+S
}
$$

で結ばれる。

また、

$$
\boxed{
|S|=u\kappa
}
$$

である。

---

# 6. κ の幾何学的意味

$$
\kappa
=
\frac{
|\sum_jz_j^2|
}{
\sum_j|z_j|^2
}
$$

とする。

\(z_j=x_j+iy_j\) とし、

$$
x=(x_1,\ldots,x_m),
\qquad
y=(y_1,\ldots,y_m)
$$

を定義する。

すると、

$$
\sum_jz_j^2
=
(\|x\|^2-\|y\|^2)
+
2i(x\cdot y).
$$

従って、

$$
\left|\sum_jz_j^2\right|^2
=
(\|x\|^2-\|y\|^2)^2
+
4(x\cdot y)^2.
$$

2×m 実行列

$$
X
=
\begin{pmatrix}
x\\
y
\end{pmatrix}
$$

の Gram 行列

$$
W
=
XX^\mathsf T
=
\begin{pmatrix}
\|x\|^2 & x\cdot y\\
x\cdot y & \|y\|^2
\end{pmatrix}
$$

の固有値を

$$
\lambda_1\ge\lambda_2\ge0
$$

とする。

固有値差より、

$$
(\lambda_1-\lambda_2)^2
=
(\|x\|^2-\|y\|^2)^2
+
4(x\cdot y)^2.
$$

また、

$$
\lambda_1+\lambda_2
=
\|x\|^2+\|y\|^2.
$$

従って、

$$
\boxed{
\kappa
=
\frac{\lambda_1-\lambda_2}
{\lambda_1+\lambda_2}
}
$$

である。

---

# 7. 複素ガウス null と Wishart 分布

各 \(x_j,y_j\) が独立標準正規なら、

$$
W=XX^\mathsf T
$$

は 2×m 実 Wishart 行列である。

その固有値同時密度は、定数因子を除き、

$$
p(\lambda_1,\lambda_2)
\propto
e^{-(\lambda_1+\lambda_2)/2}
(\lambda_1\lambda_2)^{(m-3)/2}
|\lambda_1-\lambda_2|.
$$

ここで、

$$
s=\lambda_1+\lambda_2,
\qquad
q=
\frac{\lambda_1-\lambda_2}
{\lambda_1+\lambda_2}
=
\kappa
$$

と置く。

すると、

$$
\lambda_1=\frac s2(1+q),
\qquad
\lambda_2=\frac s2(1-q).
$$

Jacobian は、

$$
\left|
\frac{\partial(\lambda_1,\lambda_2)}
{\partial(s,q)}
\right|
=
\frac s2.
$$

また、

$$
\lambda_1\lambda_2
=
\frac{s^2}{4}(1-q^2),
$$

$$
\lambda_1-\lambda_2
=
sq.
$$

したがって、

$$
p(s,q)
\propto
e^{-s/2}s^{m-1}
\cdot
q(1-q^2)^{(m-3)/2}.
$$

\(s\) と \(q\) が分離するので、

$$
\boxed{
f_\kappa(q)
=
(m-1)q(1-q^2)^{(m-3)/2},
\quad
0\le q\le1
}
$$

を得る。

---

# 8. Beta 則

$$
Y=\kappa^2
$$

と置けば、

$$
\boxed{
Y
\sim
\operatorname{Beta}
\left(
1,\frac{m-1}{2}
\right)
}
$$

である。


なお \(m=1\) では

$$
\boxed{\kappa\equiv1}
$$

となる退化例外であり、上記 Beta 則の適用範囲は

$$
\boxed{m\ge2}
$$

である。

独立検証では、

- \(m=2,3,4,5,6,39\)
- 各 \(10^6\) 標本

の KS 検定で全て棄却されなかった。

---

# 9. CDF・平均・中央値・全モーメント

CDF：

$$
\boxed{
F_\kappa(q)
=
1-(1-q^2)^{(m-1)/2}
}
$$

平均：

$$
\boxed{
E[\kappa]
=
\frac{\sqrt\pi}{2}
\frac{\Gamma((m+1)/2)}
{\Gamma((m+2)/2)}
}
$$

中央値：

$$
\boxed{
\operatorname{median}(\kappa)
=
\sqrt{
1-2^{-2/(m-1)}
}
}
$$

全モーメント：

$$
\boxed{
E[\kappa^r]
=
\frac{
\Gamma((m+1)/2)
\Gamma(1+r/2)
}{
\Gamma((m+r+1)/2)
}
}
$$

である。

---

# 10. m=2 の厳密値

\(m=2\) では、

$$
\boxed{
E[\kappa]=\frac{\pi}{4}
}
$$

および、

$$
\boxed{
\operatorname{median}(\kappa)
=
\frac{\sqrt3}{2}
}
$$

を得る。

低位数値実験の

$$
0.784638
$$

および

$$
0.866433
$$

は、この厳密値の標本推定であった。

---

# 11. 大 m 漸近

CDF より、

$$
P(\sqrt m\,\kappa\le x)
=
1-
\left(
1-\frac{x^2}{m}
\right)^{(m-1)/2}.
$$

従って、

$$
\boxed{
\sqrt m\,\kappa
\Rightarrow
\operatorname{Rayleigh}(1)
}
$$

であり、

$$
\boxed{
E[\kappa]
\sim
\sqrt{
\frac{\pi}{2m}
}
}
$$

となる。

---

# 12. 二乗ゼロ閉包の null 確率

$$
P(\kappa\le\epsilon)
=
1-(1-\epsilon^2)^{(m-1)/2}.
$$

従って、

$$
\boxed{
P(\kappa\le\epsilon)
=
\frac{m-1}{2}\epsilon^2
+
O(\epsilon^4)
}
$$

である。

よって、

$$
\boxed{
P(\kappa=0)=0
}
$$

であり、厳密二乗ゼロ閉包は generic 複素ガウス null において測度ゼロ境界である。

---

# 13. N=40 の null 基準

N=40 では \(m=39\)。

従って、

$$
E[\kappa]
=
\frac{\sqrt\pi}{2}
\frac{\Gamma(20)}{\Gamma(20.5)}
\approx
0.1994.
$$

これは Phase B/C で用いる第一の null 基準値である。

---

# 14. 再現性

本稿の正本較正パッケージには、主較正コードに加えて独立検算用 `verify_beta_law.py` を同梱する。

同スクリプトは、以下を第三者が再実行できる形で検証する。

1. **同一 seed 再生成**  
   seed=20260914 で `abs_P1` と \(\kappa\) の系列を再生成し、保存済み全10000行と突合する。

2. **別実装・別 seed 対照検算**  
   seed=12345 の独立検算ループで、不変性・\(K=N-1\) Newton 復元・平均/中央値 \(\kappa\) を再評価する。

3. **Wishart 恒等式検算**  
   ランダム構成に対して

$$
\kappa
=
\frac{\lambda_1-\lambda_2}{\lambda_1+\lambda_2}
$$

を直接比較する。

4. **Beta 則の KS 検定**  
   \(m=2,3,4,5,6,39\) について、既定では各 \(10^6\) 標本を生成し、

$$
\kappa^2
\sim
\operatorname{Beta}
\left(
1,\frac{m-1}{2}
\right)
$$

を KS 検定する。同時に厳密平均・中央値とも突合する。

出力は、

- `beta_law_validation.json`
- `beta_law_ks_summary.csv`
- `independent_crosscheck_summary.csv`

として保存する。

従って、本稿 §14 の再現性主張は、査読時の外部セッションに依存せず、正本パッケージの `run_all.sh` 一発で再現できる。

---

# 15. 物理的含意

本稿で得た結果は、

$$
\boxed{
\text{匿名内部観測の完全測定器}
+
\text{generic null の厳密基準}
}
$$

の二つである。

特に \(\kappa=0\) は generic null の測度ゼロ境界なので、自己無撞着力学がその近傍を選択する場合、それ自体が強い構造選択の証拠になる。

---

# 16. 次の判別実験

本稿では、終端実現多様体の次元をまだ主張しない。

次の実験では、保存済みの自己無撞着関係波データを用い、N を固定して、

$$
\boxed{
\dim
\left(
\mathcal M_{\rm terminal}(N)/G_B
\right)
\stackrel{?}=4
}
$$

を直接測定する。

同時に、

$$
\kappa_B(t)=|P_2(t)|
$$

を追跡し、generic null 分布からの逸脱を定量化する。

---

# 17. 結論

匿名 Bob-star に対し、正規化冪和に基づく完全不変量フレームを構成し、

$$
\boxed{
\Phi_B^{(N-1)}
}
$$

が \(P_1\neq0\) chart 上で Bob-star を \(G_B\) ゲージを除いて完全に復元することを示した。

この測定装置を低位 N で数値較正し、独立実装でも再現した。

さらに generic 複素ガウス null に対し、

$$
\boxed{
\kappa^2
\sim
\operatorname{Beta}
\left(
1,\frac{m-1}{2}
\right)
}
$$

という厳密分布を得た。

これにより、平均・中央値・全モーメント・大 m 漸近が閉形式で決まり、\(\kappa=0\) が generic null の測度ゼロ境界であることが示された。

次段階では、この較正済み測定器を保存済み自己無撞着力学へ適用し、

$$
\dim
\left(
\mathcal M_{\rm terminal}(N)/G_B
\right)
$$

および

$$
\kappa_B(t)
$$

を直接測定する。

---

# 自己引用

**[K1]** 木原範昭, 「自己無撞着な関係波閉鎖系におけるインフレーション的急拡大の機構」, 2026. Version DOI: 10.5281/zenodo.22176949; Concept DOI: 10.5281/zenodo.22112008.

**[K7]** 木原範昭, 関係波閉鎖系シリーズ論文7「ロック後の幾何（Grassmann面移動・4次元span・2非零主角・SO(4)棄却）」, 2026. Concept DOI: 10.5281/zenodo.22729107.

**[K9]** 木原範昭, 関係波閉鎖系シリーズ論文9「無名性と情報の最小単位」, 2026. Concept DOI: 10.5281/zenodo.22729127.

**[KG1]** 木原範昭, 「匿名頂点状態生成幾何 ― 匿名内部観測者から不変な関係量の体系」, v0.6, 2026. DOI / Concept DOI: （取得後記入）.

---

# 再現資料

本稿の数値較正・厳密分布検算は、以下の Google Drive フォルダを正本とする。

`匿名頂点状態生成幾何-匿名内部観測者から不変な関係量の体系/低位N完全フレーム再較正_v06_complete_v2/`

同フォルダには、正本プログラム、`run_all.sh`、全標本、集計、κ厳密値検算、分析、README、SHA256 を収録する。
