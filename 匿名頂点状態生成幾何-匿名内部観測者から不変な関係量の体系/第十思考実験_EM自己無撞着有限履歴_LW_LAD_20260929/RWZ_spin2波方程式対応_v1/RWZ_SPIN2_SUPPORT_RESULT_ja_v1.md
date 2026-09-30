# RWZ spin-2 波方程式対応 v1

## 目的
既存の重力波近似放射式を使わず、Schwarzschild 背景の実際の GR 線形摂動方程式である Regge–Wheeler / Zerilli–Moncrief 系を時間領域で解き、同じ master field から無限遠への放射と事象地平面への吸収を同時に計算する。

## 実装方程式
Schwarzschild tortoise 座標 r_* で

(-d_t^2 + d_{r_*}^2 - V_l) Psi_lm = 0

を直接時間発展する。

Odd parity (Regge–Wheeler):

V_RW = f [l(l+1)/r^2 - 6M/r^3],  f=1-2M/r.

Even parity (Zerilli–Moncrief):

n=(l-1)(l+2)/2,

V_Z = 2 f [n^2(n+1)r^3 + 3n^2Mr^2 + 9nM^2r + 9M^3]
      / [r^3 (nr+3M)^2].

r(r_*) は Lambert W を用いて Schwarzschild tortoise relation を逆変換した。

## flux
Moncrief master-field 規格化に対応する共通係数

C_l = (l+2)! / [(l-2)! 64 pi]

を用い、時間発展中の master-field energy current から左右境界へ出る flux を積分した。左境界は horizon 側、右境界は future-null-infinity 側の数値近似である。

## l=2 実行結果
初期 Gaussian displacement: x0=5, sigma=7, amplitude=1e-3.

### Regge–Wheeler odd
- horizon: 0.3446364009
- infinity: 0.6553430101
- domain remaining: 5.43e-6
- total closure: 0.9999848402

### Zerilli even
- horizon: 0.3167411406
- infinity: 0.6832385452
- domain remaining: 5.42e-6
- total closure: 0.9999851086

したがって、同じ spin-2 master field の時間発展から吸収と放射が同時に得られている。

## 放射=吸収の探索
初期 Gaussian の中心 x0 のみを変更して、その他の方程式・境界条件は固定した。

### odd parity
x0 ~= 2.39014 で
- horizon = 0.49997667
- infinity = 0.49998529

### even parity
x0 ~= 1.89697 で
- horizon = 0.50012628
- infinity = 0.49983570

したがって、少なくとも RWZ の spin-2 scattering problem には、同じ master equation の中で horizon absorption と infinity radiation がほぼ等分される状態が存在する。

## 重要な境界
これは quadrupole formula ではない。一方で full nonlinear two-body Einstein equation の厳密解でもない。Schwarzschild exact background 上の一次 GR perturbation equation を直接解いている。

次段階では、この RWZ master field を現在の有限履歴/倍音主状態へ接続し、aa, ab, bb の二次状態から source / boundary data を生成する。その後 Kerr 背景へ進む場合は s=-2 Teukolsky equation を実装する。
