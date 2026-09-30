# Reissner-Nordström 背景上の coupled spin-2 / spin-1 放射実験 v1

## 目的

既存の Schwarzschild RWZ 実験系を、電荷を持つ Reissner-Nordström (RN) 背景上の coupled gravitational-electromagnetic perturbation 系へ拡張し、同じコード・同じ初期波束で Q/M を変えるだけで、spin-2 と spin-1 の放射・吸収・相互変換が連続的に現れるかを検証する。

本実験は full nonlinear radiating Einstein-Maxwell solution ではない。背景 RN metric は exact Einstein-Maxwell solution であり、その上で用いる perturbation equations は linearized Einstein-Maxwell 系の標準的な Moncrief/Zerilli/Chandrasekhar master equations である。

## 実装した方程式

odd-parity sector の二成分 master field を

\[
\mathbf\Psi=(G,E)^T
\]

とし、

\[
(-\partial_t^2+\partial_{r_*}^2)\mathbf\Psi
-\mathbf V(r)\mathbf\Psi=0
\]

を時間発展した。

背景は

\[
f(r)=1-\frac{2M}{r}+\frac{Q^2}{r^2},
\qquad
\frac{dr_*}{dr}=f^{-1}.
\]

\(L=\ell(\ell+1)\), \(\mu^2=(\ell-1)(\ell+2)\) として、物理的 gravitational / electromagnetic basis の対称 potential matrix を

\[
V_{GG}=f\left(\frac{L}{r^2}-\frac{6M}{r^3}+\frac{4Q^2}{r^4}\right),
\]

\[
V_{EE}=f\left(\frac{L}{r^2}+\frac{4Q^2}{r^4}\right),
\]

\[
V_{GE}=V_{EG}
=f\frac{2Q\sqrt{\mu^2}}{r^3}
\]

とした。

この 2x2 matrix の固有値は、文献で知られている decoupled RN odd-parity potentials

\[
V_j=\frac{\Delta}{r^5}
\left[Lr-q_k+\frac{4Q^2}{r}\right]
\]

と一致する。ここで

\[
q_1=3M+\sqrt{9M^2+4Q^2(\ell-1)(\ell+2)},
\]

\[
q_2=3M-\sqrt{9M^2+4Q^2(\ell-1)(\ell+2)}.
\]

実装上の固有値最大誤差は全ケースで約 1e-16 以下であった。

## Q -> 0 極限

Q=0 では結合項 V_GE が厳密に 0 となる。

重力初期波束では electromagnetic flux は 0、電磁初期波束では gravitational flux は 0 となった。

したがって同一コードが Q=0 で自動的に

- Regge-Wheeler spin-2 sector
- Maxwell spin-1 sector

へ分離することを確認した。

## 実験条件

- M = 1
- ell = 2
- Q/M = 0, 0.3, 0.6, 0.9, 0.99
- initial channel = grav / em の双方
- Gaussian displacement pulse
- x0 = 5
- sigma = 7
- amplitude = 1e-3
- tortoise-coordinate domain approximately [-40,220]
- tmax = 360

各 run で horizon 側と infinity 側の gravitational / electromagnetic canonical master-field flux を独立に積分した。

## 結果

### gravitational initial pulse

| Q/M | GW->EM total conversion | same-channel GW | closure |
|---:|---:|---:|---:|
| 0.00 | 0.000000 | 0.999983 | 0.999988 |
| 0.30 | 0.026658 | 0.973325 | 0.999988 |
| 0.60 | 0.102052 | 0.897929 | 0.999987 |
| 0.90 | 0.206083 | 0.793894 | 0.999984 |
| 0.99 | 0.220886 | 0.779088 | 0.999982 |

ここで conversion は horizon + infinity の EM channel flux の合計である。

### electromagnetic initial pulse

| Q/M | EM->GW total conversion | same-channel EM | closure |
|---:|---:|---:|---:|
| 0.00 | 0.000000 | 0.999976 | 0.999983 |
| 0.30 | 0.020556 | 0.979419 | 0.999982 |
| 0.60 | 0.078206 | 0.921767 | 0.999981 |
| 0.90 | 0.155013 | 0.844953 | 0.999976 |
| 0.99 | 0.163836 | 0.836125 | 0.999972 |

## 主要結果

### 1. spin-2 <-> spin-1 conversion は Q=0 で消える

これは coupling が背景電荷 Q によってのみ生じることと一致する。

### 2. Q/M の増加とともに conversion が連続的に増える

同じ方程式・同じコードで、背景電荷を変更するだけで gravitational / electromagnetic mixing が自然に現れた。

### 3. 放射先は horizon と infinity の双方

例えば Q/M=0.9 の gravitational initial pulse では、初期 canonical energy に対して

- horizon GW: 0.087631
- horizon EM: 0.069997
- infinity GW: 0.706262
- infinity EM: 0.136086

となった。

したがって一つの coupled master field から、吸収・放射・spin conversion が同時に得られる。

### 4. energy closure

全 run で

\[
E_{domain}+E_H^{GW}+E_H^{EM}+E_\infty^{GW}+E_\infty^{EM}
\]

を initial canonical energy と比較し、closure は 0.99997 以上であった。

## 今回の意味

今回初めて、静的 exact RN / double-RN と、電荷なし RWZ を別々に扱う段階から進み、

\[
\boxed{\text{charged exact RN background}
+\text{coupled spin-2/spin-1 radiation}}
\]

を同じ実験系で時間発展できるようになった。

これにより

\[
Q=0
\rightarrow
\text{pure spin-2 / spin-1 separation}
\]

から

\[
0<Q<M
\rightarrow
\text{gravito-electromagnetic conversion}
\]

までが一つの実装で連続につながった。

## 注意点

1. これは full nonlinear radiating Einstein-Maxwell two-body exact solution ではない。
2. RN background 自体は exact solution である。
3. 放射 equations は linearized Einstein-Maxwell equations を master variables に還元した標準的 perturbation equations である。
4. 本 v1 は odd parity sector を実装した。even parity は isospectral だが、物理的 field reconstruction を含めて別途追加する価値がある。
5. 本表の flux は symmetric master system の canonical flux であり、遠方での実際の SI 単位の gravitational/electromagnetic luminosity への変換は別の normalization/reconstruction が必要である。
6. Q/M >= 1 の over-extremal RN はブラックホール horizon scattering ではなくなるため、本 v1 の horizon-flux 実験は Q/M < 1 に限定した。

## 本論への接続

この結果は、万能相互作用候補が要求する

\[
\boxed{\text{質量/電荷比を変えるだけで spin sector の寄与が連続的に変わる}}
\]

という構造を標準 Einstein-Maxwell perturbation theory 側で実際に持つことを示す。

次に必要なのは、この coupled RN master state を有限履歴・倍音主状態から直接生成し、さらに aa, ab, bb の relation sector と対応付けることである。
