# AB 回転不変量から theta を生成する派生トイモデル v1

## 目的

既存の系統 A が外部入力していた反射率 $R$ を廃止し、保存中の二つの複素波
$A(\tau),B(\tau)$ だけから散乱角 $\theta$ を生成する。

初期倍音関係は状態の初期条件として残すが、コマンドラインにもコード定数にも
散乱用の $R$ または $\theta$ は置かない。

## 更新式

各衝突時点で現在の AB 波から $\theta$ を再計算し、同じ時点の波へ

$$
\begin{pmatrix}
A'\\
B'
\end{pmatrix}
=
\begin{pmatrix}
\cos\theta & -\sin\theta\\
\sin\theta & \cos\theta
\end{pmatrix}
\begin{pmatrix}
A\\
B
\end{pmatrix}
$$

を作用させる。

この実直交回転は、各 FFT ビン $k$ の AB 合成パワー

$$
P_k^{AB}=|\widehat A_k|^2+|\widehat B_k|^2
$$

を保存する。したがって、固定したビン集合上のパワー和から作る $\theta$ も、
相互作用の反復中に保存される。

## v1 の theta 読出し

現行状態生成器には $q_A=+1,q_B=-1,p_0=1$ の搬送波が含まれる。
このため内在的な奇数倍音は、生の FFT 上では偶数ビンへ一つずれる。

基本波が占める低域 $|k|\leq2$ を除き、偶数の高域ビンをフェルミオン関係候補
$\mathcal F$ とする。

$$
\mathcal F=\{k:\ |k|\geq4,\ |k|\ {\rm is\ even}\}
$$

$$
P_f=\sum_{k\in\mathcal F}P_k^{AB},
\qquad
P_b=\sum_{k\notin\mathcal F}P_k^{AB}
$$

$$
\theta=\operatorname{atan2}\!\left(\sqrt{P_f},\sqrt{P_b}\right)
$$

$$
R=\sin^2\theta=\frac{P_f}{P_f+P_b},
\qquad
T=\cos^2\theta
$$

ここで $\mathcal F$ の選び方は v1 のトイ仮説であり、既存実験から確定した
粒子統計の定理ではない。ただし、選択後の $P_f,P_b,\theta$ が AB 回転で
保存されること自体は数値モデル内で直接検証できる。

## 初期条件

比較可能な三条件を組み込んでいる。

1. `fundamental_control`: $A=(1),B=(1)$
2. `even_boson_control_B62`: $A=(1),B=(1,2,4,\ldots,62)$
3. `odd_fermion_candidate_B63`: $A=(1),B=(1,3,5,\ldots,63)$

`--high-n` を変更すれば、初期倍音数だけを変更できる。これは散乱率の外部入力ではなく、
残す必要がある初期状態の指定である。

## 保存量と正規化

既存系統 A の各チャネル個別正規化は使用しない。実直交回転そのものが

$$
\lVert A'\rVert^2+\lVert B'\rVert^2
=
\lVert A\rVert^2+\lVert B\rVert^2
$$

および非共役二乗和

$$
\sum_i (A'_i)^2+\sum_i (B'_i)^2
=
\sum_i A_i^2+\sum_i B_i^2
$$

を保存するためである。後者が初期値でゼロなら、ゼロ閉塞もそのまま保存される。

## 依存範囲

元の `波の情報読出し/` 配下は変更しない。この派生ランナーが読むのは、
隔離環境内の次の無修正コピーだけである。

```text
ab_invariant_theta_toy_v1/run_ab_invariant_theta_toy_v1.py
  -> ../20260715/run_system_A_localization_exchange_R_sweep_preliminary_v1.py
  -> ../20260713/run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1.py
  -> ../20260711/...result_v2.json
```

したがって `対照実験_波束収縮_実行環境_v1/` 全体の中で依存関係が閉じる。
ただし、この子フォルダだけを単独で移動すると兄弟フォルダへの参照が失われる。

## 実行

次の既定条件で三条件を走らせる。

```bash
python3 ab_invariant_theta_toy_v1/run_ab_invariant_theta_toy_v1.py \
  --high-n 63 \
  --max-collision 32
```

生成物は既定で `ab_invariant_theta_toy_v1/result_v1/` に保存される。

## 初回実行結果（2026-07-29）

`--high-n 63 --max-collision 32` で実行した。

| 条件 | 生成された $\theta$ | 生成された $R$ | 最大 $\theta$ ドリフト | 判定 |
|---|---:|---:|---:|---|
| `fundamental_control` | 0 | 0 | 0 | PASS |
| `even_boson_control_B62` | 0 | 0 | 0 | PASS |
| `odd_fermion_candidate_B63` | 0.761952071831 | 0.4765625 | $2.22\times10^{-16}$ | PASS |

奇数倍音B63条件では、2回目の衝突で

$$
P(B_0\rightarrow A_2)=0.997802734375
$$

$$
P(A_0\rightarrow B_2)=0.997802734375
$$

となり、A/B起源成分がほぼ完全に交換された。

同条件の32衝突中の最大保存誤差は、ABノルムで
$5.78\times10^{-15}$、非共役二乗和で
$7.70\times10^{-17}$、AB合成スペクトルの各ビンで
$6.94\times10^{-18}$ だった。

初期の非共役二乗和自体も絶対値
$3.50\times10^{-16}$ で数値的なゼロ閉塞にあり、その後も保存された。

この内生値 $R=0.4765625$ は v1 の偶奇セクター定義から直接出た値であり、
$R_{137}$ へ合わせる調整は行っていない。

## 図化（2026-07-29）

`make_ab_invariant_theta_figures_v1.py` でPNGとSVGを生成した。

1. `result_v1/figures_v1/fig1_theta_R_transfer_invariants_v1`
   - 三条件の $\theta$ と $R$
   - 奇数倍音B63条件の反復移乗
   - $\theta$、ABノルム、二乗閉塞、AB合成スペクトルの数値誤差
2. `result_v1/figures_v1/fig2_odd_B63_waveform_harmonic_exchange_v1`
   - 衝突0・1・2回後のA/B波形
   - 同時点のA/B倍音分布
   - 2衝突後に局在波形と高次倍音構造がBからAへ、基本波構造がAからBへ
     ほぼ交換される様子

図は次のコマンドで再生成できる。

```bash
python3 ab_invariant_theta_toy_v1/make_ab_invariant_theta_figures_v1.py
```

## 広域ボゾンA＋局在ボゾン／フェルミオンBの長時間比較

衝突0–256について、次の2条件を同じ軸で比較した。

$$
(b:w,\ b:n)
$$

$$
(b:w,\ f:n)
$$

数値モデル上の対応は次のとおり。

1. 広域ボゾンA＋局在ボゾンB:
   $A=(1),B=(1,2,4,\ldots,62)$
2. 広域ボゾンA＋局在フェルミオン候補B:
   $A=(1),B=(1,3,5,\ldots,63)$

長時間データは `result_longrun_v1/`、比較図は
`result_longrun_v1/comparison_figures_v1/` に保存した。

1. `fig3_boson_boson_vs_boson_fermion_dynamics_v1`
   - B起源成分のAへの移乗と走行平均
   - A/Bの局在度
   - A/Bの有効倍音数
2. `fig4_boson_boson_vs_boson_fermion_waveforms_v1`
   - 衝突前、1衝突後、2衝突後、尾部時間平均の波形
   - 上段がボゾン同士、下段がボゾン・フェルミオン

結果は単純な「両方が定常状態へ収束」ではなかった。

- ボゾン同士は $\theta=0$ のため、衝突0から波形・局在度・倍音分布が不変。
- ボゾン・フェルミオンは固定 $\theta$ の損失なし回転なので、256衝突後も
  A/B交換が持続し、瞬時状態は定常化しない。
- 衝突193–256でも、B起源成分のAへの瞬時移乗率は
  $0.000145059939494$ から $0.999463216261$ まで振動した。
- 一方、衝突1–256の走行平均は $0.499604830752$ となり、
  時間平均だけが約 $1/2$ へ定常化した。

したがって図の右端は「収束後の瞬時波形」ではなく、衝突193–256の
時間平均波形として明示している。瞬時定常化を得るには、現在の直交回転にはない
散逸・位相平均・観測読出しなどの追加機構が必要になる。

再生成コマンド:

```bash
python3 ab_invariant_theta_toy_v1/run_ab_invariant_theta_toy_v1.py \
  --high-n 63 \
  --max-collision 256 \
  --output-dir ab_invariant_theta_toy_v1/result_longrun_v1

python3 ab_invariant_theta_toy_v1/make_boson_boson_vs_boson_fermion_comparison_v1.py
```

## 元論文 R=0.70 図の再現基準

新トイモデルとの比較基準として、元20260713エンジンの
$R=0.70$ 波形発展図を隔離ミラー環境から再生成した。

条件は元論文と同一である。

- $A$: 基本波 $N=1$
- $B$: 等振幅奇数倍音 $1,3,\ldots,63$
- 衝突回数: $0,1,2,3,5,10,20,42$
- 各A/B波形を各時点の自身の最大値で正規化
- 4行×2列、160 dpi

生成物は `reference_R070_reproduction_v1/` に保存した。
元論文PNGとの比較結果は次のとおり。

- バイト一致: `True`
- 画素配列一致: `True`
- 最大画素絶対差: `0.0`

したがって、隔離ミラー環境は元論文の $R=0.70$ 波形発展図を完全再現できる。

## 初期状態探索による比較図

$R=0.50$ と次の $R_\alpha$ については、外部 $R$ を散乱へ与える図を
比較対象として採用しない。

$$
R_\alpha
=1-\sqrt{\frac{4\pi}{137.035999084}}
=0.697177879128\ldots
$$

`inverse_initial_conditions_v1/` の独立探索器が、Aを単位基本波、
Bを等振幅奇数倍音 $1,3,\ldots,63$ に固定したまま、Bの初期振幅だけを
ブラックボックス二分探索する。候補ごとに、このファイルの
`theta_from_ab(A,B)` を無変更で呼び、その返り値と目標値の差を評価する。

探索後の前進処理が受け取るのは完成した初期A/B配列だけである。
目標 $R$ や目標 $\theta$ は受け取らず、各衝突時に現在のAB波から
$R$ と $\theta$ を読み直す。

| 比較目標 | 探索されたB初期振幅 | 初期AB波から読んだR | 誤差 |
|---:|---:|---:|---:|
| $0.50$ | $1.05045146287739$ | $0.499999999999824$ | $1.76\times10^{-13}$ |
| $R_\alpha$ | $1.65042835030908$ | $0.697177879128463$ | $2.16\times10^{-13}$ |

散乱本体 `run_ab_invariant_theta_toy_v1.py` には探索用コードを追加していない。
前進42衝突までの最大Rドリフトは、それぞれ
$2.78\times10^{-16}$、$3.33\times10^{-16}$ で、保存量判定はいずれもPASSだった。
同じ4行×2列形式のPNG/SVGは `inverse_initial_conditions_v1/` に保存した。

## 厳密有限位数根の1500衝突対照

以前の論文で $U^n=I$ から得られた二つの根

$$
R_{124,23}
=\cos^2\left(\frac{23\pi}{124}\right)
=0.697177927556659290\ldots
$$

$$
R_{620,117}
=\cos^2\left(\frac{117\pi}{620}\right)
=0.687825191114145187\ldots
$$

について、`inverse_initial_conditions_v1/exact_finite_order_roots_longrun_v1/`
で初期状態を逆算し、1500衝突まで実行した。探索されたB63波の初期振幅は
それぞれ $1.6504285637713139$ と $1.6101656970067282$ であり、初期AB波から
読み出したRは倍精度で各目標値と一致した。

散乱本体は変更しておらず、前進処理には完成した初期A/B配列だけを渡す。
現在の実直交回転では同じRに対する閉鎖周期が248衝突と1240衝突になり、
最初の周期での初期AB波への規格化回帰残差は
$5.19\times10^{-15}$ と $9.11\times10^{-14}$ だった。

図は指定点
$0,1,2,3,5,10,20,42,100,200,300,500,1000,1500$ に加え、
厳密回帰点248と1240を含む。全1501時点の診断図とCSVも同じフォルダに保存した。

## 判定境界

- 確認対象: 外部散乱パラメータがないこと、AB 合成スペクトル・$\theta$・二乗和の保存、
  チャネル間の波形移乗、$L$ と $N_{\rm eff}$ の変化。
- まだ確認しないこと: この偶奇分割が実在のボゾン・フェルミオンを一意に分類すること。
- $\theta$ の保存は反復中に同じ $U$ を与えるが、有限整数 $n$ に対する
  $U^n=I$ は $\theta=2\pi m/n$ をさらに要するため、この v1 では仮定しない。
