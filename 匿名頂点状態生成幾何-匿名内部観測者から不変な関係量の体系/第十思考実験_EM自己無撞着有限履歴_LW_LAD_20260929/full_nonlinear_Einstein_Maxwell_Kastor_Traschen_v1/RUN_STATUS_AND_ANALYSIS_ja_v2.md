# full nonlinear Einstein-Maxwell 二体系：実行状況と結果

## 1. 今回実際に実行したもの

今回この計算環境で実際に実行したのは、二中心 Kastor-Traschen (KT) 厳密解である。

\[
 ds^2=-\Omega^{-2}dt^2+a(t)^2\Omega^2(dx^2+dy^2+dz^2),
\]

\[
 \Omega=1+\sum_A \frac{M_A}{a(t)r_A},\qquad a(t)=e^{Ht},\qquad \Lambda=3H^2.
\]

各中心は \(|Q_A|=M_A\) を持つ。今回は

\[
M_1=M_2=1,\quad Q_1=Q_2=1,\quad d=8,\quad H=-0.05
\]

とした。\(H<0\) は収縮枝であり、二中心は背景 de Sitter スケールに対して接近する。

KT は Einstein-Maxwell-\(\Lambda\) 方程式の full nonlinear exact solution であり、線形摂動ではない。

## 2. 3+1 変数

\[
W=a\Omega=a+\sum_A\frac{M_A}{r_A}
\]

と置くと

\[
\gamma_{ij}=W^2\delta_{ij},\qquad \alpha=\frac{a}{W},\qquad \beta^i=0,
\]

\[
K_{ij}=-H\gamma_{ij},\qquad K=-3H.
\]

電場は符号規約を除いて

\[
E_i=-\partial_i\ln W,
\qquad
E^i=-\frac{\partial_i W}{W^3}.
\]

## 3. 制約式の数値検証

格子上で Hamiltonian constraint と Maxwell Gauss constraint を独立に有限差分評価した。

Hamiltonian constraint:

\[
{}^{(3)}R+K^2-K_{ij}K^{ij}-16\pi\rho_{EM}-2\Lambda=0.
\]

Gauss constraint:

\[
D_iE^i=0
\]

（点電荷位置を除く）。

格子サイズは N=33,49,65。中心から半径 1.2 未満と外側境界1セルを評価から除外した。

### t=0 の RMS residual

| N | h | Hamiltonian RMS | Gauss RMS |
|---:|---:|---:|---:|
| 33 | 0.75 | 5.5740e-3 | 1.3935e-3 |
| 49 | 0.50 | 1.9527e-3 | 4.8819e-4 |
| 65 | 0.375 | 1.1803e-3 | 2.9507e-4 |

推定収束次数は successive grids で約 2.59, 1.75。点源近傍を固定物理半径で除外した単純2次差分として、概ね2次収束と整合する。

Momentum constraint は

\[
K_{ij}=-H\gamma_{ij}
\]

かつ Poynting flux がゼロなので解析的に 0 である。

## 4. 動的二中心

背景物理距離 \(a(t)d\) は

| t | a(t) | a(t)d |
|---:|---:|---:|
| 0 | 1.000000 | 8.000000 |
| 4 | 0.818731 | 6.549846 |
| 8 | 0.670320 | 5.362560 |

へ減少した。

したがって今回の実行は、静的 double-RN ではなく full nonlinear かつ動的な二中心 Einstein-Maxwell 解を実際に評価している。

## 5. aa, ab, bb の非線形同時存在

中点で

\[
W^2=(a+u_a+u_b)^2
\]

を展開すると

\[
W^2=a^2+2au_a+2au_b+u_a^2+2u_au_b+u_b^2.
\]

したがって

\[
aa=u_a^2,\qquad ab=2u_au_b,\qquad bb=u_b^2
\]

は一つの full nonlinear metric factor 内に同時に存在する。

直接展開値と \(W^2\) の差は全時刻で最大

\[
2.22\times10^{-16}
\]

であった。

## 6. 重要な境界：今回まだ実行していないもの

今回 **実行していない** のは、asymptotically-flat な一般 charged binary merger を 3D numerical relativity で時間発展し、\(\Psi_4\) と Maxwell Newman-Penrose scalar を null infinity / extraction sphere で抽出する計算である。

理由は物理ではなく現在の実行環境である。このセッションには Einstein Toolkit / GRChombo / SpECTRE と MPI runtime がインストールされておらず、利用可能資源は 5 CPU core、約 5.8 GiB RAM である。production NR の charged binary evolution をここで実行済みと主張することはできない。

Einstein Toolkit の現行 Hypatia release は 2026-07-10 公開で、BSSN/Z4c spacetime evolution と matter support を含む。公式 gallery の比較的大きな merger 例では 256 cores / 140 GB 程度が記載されており、一般3D merger はこの環境の規模を大きく超える。

## 7. 本研究に対する現在の到達点

実行済み：

1. static exact Einstein-Maxwell: RN / MP / double-RN
2. dynamic exact full nonlinear Einstein-Maxwell-Lambda: two-center Kastor-Traschen（今回）
3. spin-2 radiation on exact Schwarzschild background: RW/Zerilli
4. coupled spin-2/spin-1 radiation on exact RN background: RN coupled perturbation
5. retarded Maxwell + Lorentz-Dirac finite-history evolution

未実行：

6. asymptotically-flat full nonlinear charged two-body merger in 3D numerical relativity

したがって、full nonlinear Einstein-Maxwell 自体への対応は今回実行した。しかし「一般 charged binary の放射まで含む full nonlinear 3D evolution」は依然として別計算資源を要する。

## 8. 参考

- D. Kastor and J. Traschen, Cosmological Multi-Black Hole Solutions, Phys. Rev. D 47, 5370 (1993), arXiv:hep-th/9212035.
- Einstein Toolkit, current release Hypatia (2026-07-10), https://einsteintoolkit.org/

## 9. 2026-09-29 後続実験：strict rules による状態生成系の全面再監査

KT 実験の後、本研究固有の `STRICT_EXPERIMENT_RULES_ja` に照らして、EM finite-history、RWZ、RN coupled、KT の各実装を再監査した。

その結果、標準理論参照計算としては有効である一方、以下の旧実装は universal-interaction state generator の正本としては扱わないことにした。

- EM finite-history: 内部 RK4 を使用。
- RWZ: previous/current の時間配列を persistent state として使用。
- RN coupled: 同じく leapfrog 型 previous/current state を使用。
- KT: exact solution の時刻サンプリングを次状態生成に使用。

これらは削除せず、標準理論側の reference / benchmark として保存する。

strict 再構築版では、永続状態を有限履歴を表す複素倍音係数だけに限定し、全次状態を固定 sparse quadratic interaction array による一回積和

\[
Z_{n+1}=\mathcal T[Z_n,Z_n]
\]

だけで同期生成する構成へ変更した。

実行ループ中の transition は本質的に

```python
mon = z[pi] * z[pj]
z_next = M @ mon
```

のみである。読出し量は transition に帰還しない。

## 10. strict harmonic RN coupled：匿名二チャネル放射実験

同一 interaction array、同一 transition のまま、初期状態の Q/M だけを変更して RN coupled benchmark を再実行した。

条件：

- Q/M = 0, 0.3, 0.6, 0.9, 0.99
- anonymous channel 0 initial / channel 1 initial の両方
- nx = 61
- harmonic history K = 4
- dt = 0.12
- 700 steps

Q=0 では他チャネル flux は 0 となり、Q>0 では Q/M の増加に伴って他チャネルへの変換が連続的に現れた。

channel 0 initial の他チャネル境界 flux 比率：

| Q/M | other-channel fraction |
|---:|---:|
| 0.30 | 0.0016139064 |
| 0.60 | 0.0062901234 |
| 0.90 | 0.0137090482 |
| 0.99 | 0.0164396584 |

channel 1 initial：

| Q/M | other-channel fraction |
|---:|---:|
| 0.30 | 0.0014597175 |
| 0.60 | 0.0057071048 |
| 0.90 | 0.0125027843 |
| 0.99 | 0.0150217465 |

生成側では二チャネルは匿名であり、重力的 spin-2 / 電磁的 spin-1 という意味付けは Q→0 極限に基づく読出しでのみ行う。

## 11. strict Q=0 放射・吸収バランス

同じ strict map のまま、初期波束位置だけを変更して inner / outer flux のバランス点を探索した。

\[
x_0\approx2.97890625
\]

で

\[
\frac{F_{out}-F_{in}}{F_{out}+F_{in}}
\approx-6.44\times10^{-5}
\]

を得た。

また strict harmonic map と外部 validation 用 direct sample recurrence を 120 step 比較し、

\[
\max|\Delta|=3.0654\times10^{-17}
\]

で一致した。

## 12. strict Kastor–Traschen full nonlinear exact sector

旧 KT 実装では exact state を時刻ごとに外部注入していたため、strict 再構築では scale state 自身を状態に含め、

\[
a' = a h
\]

を同じ quadratic map の内部更新として実装した。

結果：

- 32 step まで外部 exact benchmark との差：0 〜 1.1e-16
- Hamiltonian / Gauss residual：N=33,49,65 の格子細分化に伴って減少
- aa / ab / bb 再構成誤差：0

したがって KT sector でも、exact solution は状態生成器ではなく外部 benchmark としてのみ使用する構成へ移行した。

## 13. strict spin hierarchy

同じ一回 quadratic map で half-angle state と ordinary-angle state を更新し、位相 weight の階層を検証した。

- 2pi sign flip error = 1.16e-14
- 4pi return error = 2.40e-14
- weight-2 covariance error = 4.48e-16

これは spin-1/2 型二重被覆と weight-2 読出しを、同一の乗法的状態更新の中で再現できることを示す数値監査である。

## 14. 追加実験：GW / EM 二読出しチャネルの放射・吸収同時バランス

さらに strict harmonic RN coupled 系について、transition と interaction array を一切変更せず、二つの初期波束位置だけを独立に変更して、二読出しチャネルの inner / outer flux が同時に釣り合う初期状態を探索した。

各チャネルについて

\[
\delta_i=
\frac{F_{i,out}-F_{i,in}}
{F_{i,out}+F_{i,in}}
\]

を定義し、\(\delta_0\simeq0\)、\(\delta_1\simeq0\) を同時に要求した。

探索変数は

\[
x_{0,\mathrm{ch0}},\qquad x_{0,\mathrm{ch1}}
\]

のみであり、初期振幅比は \(\eta=-1\) に固定した。

共通条件：

- nx = 61
- K = 4
- dt = 0.12
- steps = 700
- sigma = 1.5
- channel 0 amplitude = 1e-4
- channel 1 amplitude = -1e-4

最終結果：

| Q/M | x0(ch0) | x0(ch1) | delta(ch0) | delta(ch1) | simultaneous residual norm |
|---:|---:|---:|---:|---:|---:|
| 0.30 | 2.9802290757 | 2.8080751966 | 3.0518e-9 | -9.1561e-11 | 3.0532e-9 |
| 0.60 | 2.9931696323 | 2.9372689819 | 7.1925e-10 | -5.5462e-11 | 7.2138e-10 |
| 0.90 | 3.0291188627 | 3.0626034735 | 3.2937e-13 | -4.4875e-14 | 3.3241e-13 |
| 0.99 | 3.0407376807 | 3.0948830400 | 1.4451e-12 | -9.8395e-13 | 1.7483e-12 |

全4ケースで inner / outer flux は非ゼロであり、両チャネルで同時にほぼ一致した。

Q→0 極限に基づく読出しを用いる限り、これは

\[
F_{GW,out}\simeq F_{GW,in}\neq0
\]

と同時に

\[
F_{EM,out}\simeq F_{EM,in}\neq0
\]

に対応する初期状態が、この離散化・境界条件・有限時間窓の中に存在することを示す。

なお最初の探索では両チャネルを同じ初期位置に拘束していたため、Q/M=0.9 で residual norm が約 5.62e-3 残った。transition を変更せず、初期波束位置を独立にしただけで

\[
5.62\times10^{-3}\rightarrow3.32\times10^{-13}
\]

まで低下した。したがって最初の失敗は高電荷側で解が存在しないことを示すものではなく、初期状態族の制約が強すぎたことによる。

## 15. 再現性監査

保存済み最適初期状態を用いて `verify_dual_balance_repro_v1.py` から独立に直接再実行したところ、4ケースすべてで保存値との差は

\[
0.0
\]

であり、`all_exact_reproduction = true` を確認した。

さらに `RUN_DUAL_BALANCE_ALL_v1.sh` を用いて探索そのものを最初から再実行し、4個の `dual_q*_v2.json` が再実行前の保存版と SHA-256 まで完全一致した。

したがって、最終残差だけでなく探索履歴を含む JSON 全体がバイト単位で再現された。

## 16. 2026-09-29 時点の実験系の位置づけ

現在は、以下の二層を明確に分離して保存する。

### A. 標準理論 reference / benchmark

1. retarded Maxwell + Lorentz-Dirac finite-history
2. Fourier/harmonic finite-history representation
3. RW/Zerilli spin-2 scattering
4. RN coupled spin-2/spin-1 perturbation
5. full nonlinear exact two-center Kastor–Traschen

これらは標準理論との対応確認には使用するが、本研究の strict universal-interaction state generator の正本ではない。

### B. strict universal-interaction implementation

1. persistent state = 複素倍音係数のみ
2. fixed sparse quadratic interaction array
3. 一回積和だけで全次状態を同期生成
4. Q/M を含むパラメータは状態側へ持たせる
5. 読出しは transition へ帰還しない
6. RN coupled anonymous two-channel conversion
7. Q=0 放射・吸収バランス
8. KT exact-sector benchmark 再現
9. spin hierarchy
10. GW / EM 二読出しチャネル同時放射・吸収バランス
11. 完全再実行による byte-exact reproduction

この strict 系を、2026-09-29 時点の本研究における再現実験の正本とする。

## 17. 未実行範囲

依然として未実行なのは、一般の asymptotically-flat full nonlinear charged binary merger を 3D numerical relativity で時間発展し、重力波・電磁波を同時抽出する production-level 計算である。

今回の strict 同時バランス実験はこの full nonlinear 3D merger を代替するものではなく、現在の離散 universal-interaction benchmark 内で二チャネル同時バランス状態が存在することを確認した段階である。
