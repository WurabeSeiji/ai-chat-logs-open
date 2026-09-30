# EM有限履歴：倍音係数を主状態とする実行結果 v3

## 目的

配列として保持していた有限履歴を persistent state から除き、有限時間窓の全履歴を Fourier/倍音係数そのものとして保持する。必要な過去状態は retarded Lienard-Wiechert 場を評価するときだけ一時的に読み出し、次状態の書込みは履歴配列を更新せず倍音係数空間で直接行う。

## 物理系

- 相互場：点電荷の retarded Lienard-Wiechert field
- 自己力：共変 Lorentz-Dirac (LAD)
- 単位系：c=1, 4 pi epsilon_0=1
- dipole radiation-power 近似は使用しない
- 重力は今回含めない

## 倍音主状態

N点履歴 x_j (j=0,...,N-1) の係数を

C_m = (1/N) sum_j x_j exp(-i 2 pi m j/N)

とする。persistent state は C_m のみとし、履歴配列は保持しない。

1 step 後に x_0 を捨て、新状態 y を末尾に書くとき、係数は直接

C'_m = exp(+i 2 pi m/N) [ C_m + (y-x_0)/N ]

で更新できる。ここで

x_0 = sum_m C_m

なので、履歴配列は不要である。

## 実行条件

対照実験 v2 と完全に同じ条件を使用した。

### radial control
- q1=+0.1, q2=-0.1
- m1=m2=1
- separation=1
- v1=v2=0
- history=12000
- dt=1e-4
- steps=50

### transverse radiation case
- q1=+0.1, q2=-0.1
- m1=m2=1
- separation=1
- v1=+0.05, v2=-0.05
- history=12000
- dt=1e-4
- steps=200

## 結果

### radial control
最終状態の配列版との差：
- r: 7.55e-15
- v: 6.48e-18
- u: 1.31e-14
- a4: 2.62e-15

最大履歴書込み/読出し誤差：1.61e-13。

### transverse radiation case
最終状態の配列版との差：
- r: 8.33e-16
- v: 1.10e-15
- u: 2.44e-15
- a4: 7.56e-16

最大履歴書込み/読出し誤差：5.24e-14。

放射場：

max |E_rad| = 5.0923943725305686e-05

共変拘束：

max |u.u + 1| = 1.05e-13

max |u.a| = 7.70e-18

## 解釈

この結果は、有限履歴を persistent な時間配列として保持しなくても、全モードの倍音係数を主状態とし、必要時に履歴を読み出すことで、同じ retarded Maxwell + LAD 力学を機械精度で再現できることを示す。

ただし今回は全 N モードを保持しているため、情報圧縮を示した実験ではない。示したのは「履歴の物理的表現を時間セル配列から倍音状態へ置換できる」という表現同値性である。

次の課題は、(1) 何モードまで削減しても retarded field と反作用が保たれるか、(2) 全モード保持を必要とせず状態生成則そのものが有限倍音系を選ぶか、を検証することである。
