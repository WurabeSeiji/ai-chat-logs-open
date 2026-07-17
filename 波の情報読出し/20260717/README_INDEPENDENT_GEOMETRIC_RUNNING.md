# 独立数値実験：4次元 `1 + 8 + 128` 幾何からの有効結合曲線

## 目的

`BH熱力学プログラム/papers/paper7_alpha/BH_Paper7_Alpha_Identity_ja.md`
の完全包含条件から137セルを再生成し、

```text
137 = 中心1 + 内層8 + 外殻128
```

を既存の Phase 4--6 交換散乱コードとは独立に検査する。
外殻が中心反転 `u -> -u` のもとで閉じ、64組の反対点対、
すなわち向き付きには `128 = 64 + 64` と分解されることも検査する。

プログラム：

```text
independent_geometric_running_137_to_128_v1.py
```

## モデル

半径3の4次元球に完全包含される整数中心の単位4次元セルを、Paper 7 の条件

$$
\sum_{j=1}^{4}\left(|c_j|+\frac12\right)^2\leq 3^2
$$

だけから列挙する。得られた集合を内核9セルと外殻128セルに分ける。

4次元の方向平均核

$$
\left\langle e^{iq\hat n\cdot r}\right\rangle_{\hat n\in S^3}
=\frac{2J_1(qr)}{qr}
$$

を使い、内核と外殻の正規化相互コヒーレンスの二乗を
$W_{\mathrm{core}}(q)$ とする。連続補正には4次元単位球の形状因子

$$
F_{\mathrm{ball}}(q)=\frac{8J_2(q)}{q^2}
$$

の二乗 $W_{\mathrm{ball}}(q)$ を用いる。

各 $q$ で

$$
\frac1{\alpha(q)}
=128+9W_{\mathrm{core}}(q)
+\frac{\pi^2}{2}\alpha(q)W_{\mathrm{ball}}(q)
$$

の正根を計算する。

## 分類

- 137セルと `1 + 8 + 128` の列挙：Paper 7 の幾何からの厳密な帰結。
- $W_{\mathrm{core}}$ の選択：内核が外部チャネルとして見える度合いに関するモデル仮説。
- $W_{\mathrm{ball}}$ の選択：連続補正の解像度依存性に関するモデル仮説。
- 自己整合式：Paper 7 と高エネルギー極限128を接続する作業仮説。
- 128.946、129.394：曲線作成後にのみ使う診断値。

## 実行

```bash
python3 independent_geometric_running_137_to_128_v1.py
```

既定の出力先：

```text
results_geometric_running_137_to_128_v1/
```

出力：

- `cell_catalog_v1.csv`：137セルの全座標と層分類。
- `distance_histograms_v1.csv`：方向平均に使う厳密な二点間距離分布。
- `geometric_running_curve_v1.csv`：$q$ ごとの形状因子、可視度、$\alpha^{-1}$。
- `geometric_running_summary_v1.json`：極限、極値、診断値との交点。
- `geometric_running_curve_v1.png`：曲線と可視度の図。

## 解釈上の注意

有限セル集合の形状因子は回折振動を持つため、得られる曲線は一般に単調ではない。
128.946との交点が存在することだけでは物理スケールの予言にならない。
無次元波数 $q$ と物理エネルギー $Q$ の対応

$$
q=\frac{Q\ell}{\hbar c}
$$

に現れるセル長 $\ell$ は、この実験では決めていない。
