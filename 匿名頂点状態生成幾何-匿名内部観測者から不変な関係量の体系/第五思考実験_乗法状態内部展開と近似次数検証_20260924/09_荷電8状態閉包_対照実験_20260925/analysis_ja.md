# 荷電8状態閉包・対照実験 v1 — 結果分析

作成日: 2026-09-25

## 1. 目的

重力-only 対照実験で第6状態として内部化できた対称質量比

$$
\mathcal N=\eta
$$

に加え、leading-order 荷電準円軌道で独立に必要となる

$$
\mathcal C=Z=1-\lambda_A\lambda_B,
$$

$$
\mathcal D=(\lambda_A-\lambda_B)^2
$$

を第7・第8状態として内部化し、外部物理引数なしで状態生成を閉じられるかを検査した。

持続状態は

$$
Z_8=(U,P,E,H,Q,\mathcal N,\mathcal C,\mathcal D).
$$

本実験は leading-order adiabatic quasi-circular sector に限定する。完全 Einstein-Maxwell 二体系や R4 の完全乗法閉包を証明するものではない。

## 2. 状態だけから用いる生成式

準円軌道では $P=r$ とし、参照式

$$
\dot r=-\frac{4}{3}\frac{\eta(\Delta\lambda)^2Z}{r^2}
        -\frac{64}{5}\frac{\eta Z^2}{r^3},
$$

$$
\dot\phi=\sqrt{\frac{Z}{r^3}}
$$

から

$$
\boxed{
\frac{dP}{d\phi}
=-\frac{4}{3}\mathcal N\mathcal D\frac{\sqrt{\mathcal C}}{\sqrt P}
-\frac{64}{5}\mathcal N\frac{\mathcal C^{3/2}}{P^{3/2}}
}
$$

および

$$
\boxed{
\frac{dt}{d\phi}=\frac{P^{3/2}}{\sqrt{\mathcal C}}
}
$$

を用いた。

したがって生成関数は外部の $\eta,Z,\Delta\lambda$ を受け取らず、旧状態中の $\mathcal N,\mathcal C,\mathcal D$ だけを読む。

## 3. 保存3状態の恒等写像

11相 one-hot 位相状態 $q_j$ に対して

$$
\mathcal N' = \sum_j q_j(1\cdot\mathcal N),
$$

$$
\mathcal C' = \sum_j q_j(1\cdot\mathcal C),
$$

$$
\mathcal D' = \sum_j q_j(1\cdot\mathcal D)
$$

とした。

全11相で自己係数は正確に1であり、外部コピーは使っていない。

## 4. 6状態・7状態が閉じない反例

初期の

$$
(U,P,E,H,Q,\mathcal N)
$$

を完全に同一にし、$P=50,\mathcal N=0.25$ としたまま荷電関係量だけを変えた。

### A: 基準

$$
\mathcal C=1.09,\quad \mathcal D=0.36
$$

$$
\frac{dP}{d\phi}=-0.0280177282391.
$$

### B: C は同じ、D だけ変更

$$
\mathcal C=1.09,\quad \mathcal D=0.4225
$$

$$
\frac{dP}{d\phi}=-0.0310937330433.
$$

A と B は最初の6状態が同じなのに次状態が異なる。したがって6状態系、および第7状態として $\mathcal C$ だけを追加する系は閉じない。

### C: D は同じ、C だけ変更

$$
\mathcal C=1.08,\quad \mathcal D=0.36
$$

$$
\frac{dP}{d\phi}=-0.0277948500093.
$$

A と C は $\mathcal D$ まで同じだが $\mathcal C$ が異なり、次状態が異なる。第7状態として $\mathcal D$ だけを追加する系も閉じない。

したがって一般のこの LO 荷電準円セクターでは、少なくとも

$$
\boxed{\mathcal N,\mathcal C,\mathcal D}
$$

の3関係量が必要であり、旧5状態に加えて8状態が必要である。

## 5. 外部引数系との厳密対照

4ケースについて4000 macro step を比較した。

- 基準: $N=0.25,C=1.09,D=0.36$
- 同じ $C$・異なる $D$
- 同符号電荷: $C=0.91,D=0$
- 非等質量: $N=0.16,C=1.09,D=0.36$

結果は全ケースで

```text
max |Delta(U,P,E,H,Q)| = 0
bitwise macro mismatch = 0
N/C/D identity failure = 0
```

であった。

さらに seed=20260925 の500ランダム状態、合計5500 microstep の監査でも

```text
max core5 difference    = 0
bitwise micro mismatch  = 0
identity micro failure  = 0
phase failure           = 0
```

であった。

## 6. 明示11-microstep と折り畳みmacroの一致

長時間基準計算を高速化する前に、明示11-microstep 状態機械と、その代数的に折り畳んだ macro 更新を12000 macro step 比較した。

8状態すべてについて

```text
max difference = [0,0,0,0,0,0,0,0]
bitwise mismatch = 0
```

であった。

したがって長時間の物理解比較では折り畳み macro を用いた。ただしこれは高速化用の対照経路であり、状態閉包そのものの合否は明示 microstep 監査で判定している。

## 7. 保存済み解析基準との比較

基準条件は

$$
\mathcal N=0.25,\quad
\mathcal C=1.09,\quad
\mathcal D=0.36,\quad
r:50\to20.
$$

8状態生成器は $r=20$ を跨ぐまで

```text
488345 macro steps
= 5,371,795 explicit-microstep equivalents
```

進んだ。

$r=20$ への交点補間では

$$
t_{8}=169032.31760704893,
$$

$$
t_{ref}=169032.31760706357,
$$

したがって

$$
\boxed{|\Delta t|=1.46\times10^{-8}}
$$

であった。

位相は

$$
\phi_{8}=767.088999401507,
$$

$$
\phi_{ref}=767.0889994015913,
$$

より

$$
\boxed{|\Delta\phi|=8.44\times10^{-11}\ \mathrm{rad}}
$$

であった。

サンプリング時系列全体では

```text
max |t - t_closed(r)|        = 1.153e-8
max |H - exp(i phi_closed)|  = 4.248e-11
max |r - r_saved(t)|         = 9.069e-9
max xy residual              = 3.754e-6
```

であった。

## 8. 恒等状態の長時間保持

全基準軌道で

```text
max |N-N0| = 0
max |C-C0| = 0
max |D-D0| = 0
```

であった。

したがって荷電相互作用を導入しても、今回の LO モデルでは

$$
\boxed{
\mathcal N'=\mathcal N,
\quad
\mathcal C'=\mathcal C,
\quad
\mathcal D'=\mathcal D
}
$$

でよい。ただし、それぞれの列は $P,Q$ 等の更新へ作用する。

## 9. 重要な判定

今回確認できたのは

$$
\boxed{
\text{LO 荷電準円セクターは、}
(U,P,E,H,Q,N,C,D)
\text{で状態閉包できる}
}
$$

ということと、

$$
\boxed{
6状態および単一の荷電状態だけを足した7状態では一般に閉じない
}
$$

ということである。

これは「荷電により $N$ が変化する」という結果ではない。必要なのは $N$ を変化させることではなく、独立な荷電関係量 $C,D$ を追加状態として明示することである。

## 10. 残る課題

本実験は次をまだ解決しない。

1. eccentric charged orbit で同じ8状態だけで閉じるか。
2. charged 1PN/2PN conservative dynamics で追加状態が必要か。
3. full Einstein-Maxwell 二体系での状態数。
4. 既存B方式内部の `P+hK` 型加法演算を含む R4 完全乗法閉包。
5. $C,D$ をさらに根源的な匿名状態から導出できるか。

したがって次段階は、8状態を固定したまま eccentricity を解放し、状態数が増えるかを監査するのが自然である。
