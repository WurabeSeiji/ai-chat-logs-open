# 追加実験総括：seedless onset、Floquet 選択、保存則、時間再パラメータ化、N=5 モジュライ

**作成日:** 2026-08-26  
**対象:** N=5 raw-K corrected dynamics（既存物理則を変更せず、解析・初期固定点精度のみ変更）

## 1. Cayley 更新の厳密保存則とポンプ枯渇

実反対称生成子 $K_t^T=-K_t$ に対する Cayley 写像

$$
C_t=(I-\gamma K_t)^{-1}(I+\gamma K_t)
$$

は各 step で厳密に実直交である。

$$
C_t^T C_t=I.
$$

したがって exact arithmetic では

$$
Z_{t+1}^\dagger Z_{t+1}=Z_t^\dagger Z_t,
\qquad
Z_{t+1}^T Z_{t+1}=Z_t^T Z_t
$$

が厳密保存される。

既存 raw-K N=5 データで実装検証すると、

- $H_{\rm total}$: 1.0000000000000000 ～ 1.0000000000000058
- $\max|H_\parallel+H_\perp-H_{\rm total}|=4.44\times10^{-16}$
- $\max|Z^T Z|=6.61\times10^{-15}$

であった。

最終 step 5000 では

$$
H_\parallel=0.3220054054,
\qquad
H_\perp=0.6779945946,
$$

であり、$H_\perp$ の増大分は外部注入ではなく親平面内成分の減少から来ている。

この結果から

$$
0\le H_\perp\le H_{\rm total}=\text{const.}
$$

が運動学的に強制され、指数増大が無限に続かないこと自体は定理となる。

## 2. $K/\sigma_{\max}$ 正規化は「有限 step で完全な時計換え」ではなかった

同一状態で

$$
C(K/\sigma,\gamma)=C(K,\gamma/\sigma)
$$

は恒等的に成立し、連続極限では正規化あり・なしのベクトル場は状態依存のスカラー倍として時間再パラメータ化関係にある。

一方、既存 N=5 の finite-step データを、提案された累積支配モード Cayley 位相

$$
\Phi=\sum_t2\arctan(\gamma\sigma_t)
$$

で重ねた結果、成長域の $\log_{10}H_\perp$ RMSE は **0.944 decade** であり、曲線は完全には重ならなかった。

したがって、現行刻み $\gamma=\tan(\pi/144)$ では

- 連続極限での時間再パラメータ化解釈は成立する
- しかし finite-step の二軌道を「完全に同一」とする主張は支持されない

と区別する必要がある。

従来表の成長率比と到達 step 比の一致は強い近似的証拠ではあるが、有限刻みの状態依存写像では離散化差が残る。

## 3. seedless onset の tol 掃引：潜伏期は固定点残差の対数で決まる

`make_parent` の tolerance のみを変更し、明示 perturbation を一切加えず $Z_0=v$ から N=5 raw-K を走らせた。

| requested tol | parent residual | $f\ge10^{-8}$ onset step | $\log f$ growth rate/step |
|---:|---:|---:|---:|
| $10^{-6}$ | $3.87\times10^{-7}$ | 72 | 0.17484 |
| $10^{-8}$ | $1.82\times10^{-9}$ | 134 | 0.17253 |
| $10^{-10}$ | $5.08\times10^{-11}$ | 176 | 0.17251 |
| $10^{-12}$ | $2.38\times10^{-13}$ | 238 | 0.17251 |

回帰

$$
t_{\rm onset}=a[-\ln(\varepsilon_{\rm parent})]+b
$$

に対し

$$
a=11.6162,\qquad R^2=0.999992
$$

を得た。一方、指数成長率は tolerance に依存せず

$$
\lambda_H\simeq0.17251/\text{step}
$$

で一定である。

これは「seedなし」の初期床が `make_parent` の有限固定点残差／浮動小数点残差であり、力学的成長率そのものは残差に依存しないことを直接示す。

## 4. rotating-frame Jacobian/Floquet 解析：rank-4 開始機構を特定

$10^{-12}$ parent を用い、親相対平衡の一 step 写像から全体位相回転を除いた rotating-frame map を構成し、20実次元（N=5, M=10 complex）の Jacobian を中心差分で評価した。

固定点 defect は

$$
\|R(v)-v\|=1.04\times10^{-14}
$$

である。

有限差分幅 $3\times10^{-6}$ ～ $10^{-7}$ の全てで、上位 multiplier は安定して

$$
\mu_1=\mu_2=1.090086569,
$$

$$
\mu_3=\mu_4=1.052603212
$$

となった。

最大不安定固有空間は **2実次元** である。したがって parent の rank-2 平面に最速不安定2方向が加わることで

$$
\boxed{\text{rank }2+2=4}
$$

が選択される。

さらに最大 multiplier の振幅指数は

$$
\lambda_A=\ln\mu_1=0.0862571143,
$$

なので二乗振幅指数の予測は

$$
2\lambda_A=0.172514229.
$$

これは tol 掃引で直接測定した

$$
\lambda_H\simeq0.1725136
$$

と一致する。

さらに onset の理論式

$$
H_\perp(t)\propto \varepsilon_{\rm parent}^2e^{2\lambda_A t}
$$

から

$$
t_{\rm onset}\sim\frac{1}{\lambda_A}[-\ln\varepsilon_{\rm parent}]+\text{const.}
$$

を予測し、

$$
\frac1{\lambda_A}=11.5932
$$

は実測傾き 11.6162 と 0.2% 程度で一致した。

したがって N=5 について、

1. parent は相対平衡
2. residual が内在 perturbation floor
3. 最大不安定2実方向が rank-4 を開く
4. Floquet/Jacobian 指数が観測された $H_\perp$ 指数を定量的に与える

という開始機構が一本につながった。

第2不安定対 $\mu=1.052603212$ は、二乗振幅指数

$$
2\ln\mu\simeq0.10253
$$

のより遅い副成長を予測する。

## 5. N=5 の相対位相モジュライ seed 掃引

8個の parent random seed で 5000 step を実行した。全 run で最終4群の群サイズは

$$
\boxed{3+3+2+2}
$$

となり、二つの距離族の modulus はともに

$$
|v_A|\simeq|v_B|\simeq0.1=1/M
$$

へ収束した。

一方、符号を商に取った相対位相

$$
\arg(v_B/v_A)\bmod\pi
$$

は run 間でおよそ

$$
-0.0985\ \text{rad} \;\text{から}\; 0.2659\ \text{rad}
$$

に分布した。

したがって現時点の8 seed では、3+3+2+2 と等 modulus は強く固定される一方、二距離族間の相対位相は一意にロックされない。これは平坦方向（modulus）の存在を支持する。ただし一部 run は step 5000 で完全収束前であり、長時間追跡での再確認が必要である。

## 6. 等分配の定量化：スペクトルエントロピー

$$
p_m=|z_m|^2/\sum_k|z_k|^2,
\qquad
S=-\sum_mp_m\ln p_m
$$

を既存 N=5 全 step データから計算した。

初期は

$$
S/\ln M=0.97472,
$$

step 375 で一度 0.97322 まで低下した後、最終 step 5000 では

$$
\boxed{S/\ln M=1.000000}
$$

となった。

したがって第2段階は振幅分布の完全等分配へ到達する。ただし $S(t)$ は厳密単調ではなく、単純な H 定理ではない。

## 7. 数理構造として論文へ格上げすべき定理

### 定理A：Cayley 保存則

実反対称 $K_t$ に対する Cayley 更新は $Z^\dagger Z$ と $Z^T Z$ を厳密保存する。

### 定理B：局所頂点閉包とヌル simplex の同値

重心ゼロの複素埋込み $x_i$ と

$$
z_{ij}^2=(x_i-x_j)\cdot(x_i-x_j)
$$

について、全頂点 star 閉包

$$
\sum_{j\ne i}z_{ij}^2=0\quad\forall i
$$

は

$$
x_i\cdot x_i=0\quad\forall i
$$

と同値である。

等 modulus $|z_{ij}^2|=1/M$ と合わせると最終状態は

$$
\boxed{\text{equimodular null complex simplex}}
$$

として記述できる。

### 定理C：N=5 の 13 と12

4群が $A_+,A_-,B_+,B_-$、サイズが $3,3,2,2$、値が符号対なら、2辺ゼロ閉包数は

$$
3\times3+2\times2=\boxed{13},
$$

exact cover 数は

$$
3!\,2!=\boxed{12}.
$$

### 定理D：N=4 の120°

3つの等 modulus 複素値が和ゼロなら、位相差は厳密に $120^\circ$ である。

## 8. $U^n=I$ についての論理修正

自己無撞着性から今回厳密に得られるのは

$$
\text{self-consistency}
\Rightarrow
\text{complex rotating pair}
\Rightarrow
\sum z^2=0
\Rightarrow
S^1\text{ compact phase orbit}
$$

までである。

$$
U^n=I
$$

にはさらに

$$
\Delta\theta/2\pi\in\mathbb Q
$$

を保証する rational phase locking / discrete phase selection が必要である。したがって $U^n=I$ は本論文ではまだ完全導出定理とはしない。

## 結論

今回の追加実験で、N=5 の急拡大開始機構についてはかなり強く閉じた。

$$
\boxed{
\text{self-consistent rank-2 relative equilibrium}
\rightarrow
\text{2D dominant unstable eigenspace}
\rightarrow
\text{rank-4}
\rightarrow
H_\perp\text{ exponential transfer}
\rightarrow
\text{bounded saturation}
\rightarrow
\text{equimodular simplex ordering}
}
$$

特に、Floquet/Jacobian から得た $2\ln\mu_1=0.172514$ と、直接時間発展の $H_\perp$ 成長率 $0.172513$ の一致、および onset-vs-residual 傾きの一致は、急拡大開始を「相対平衡の線形不安定性」として定量的に特定する結果である。
