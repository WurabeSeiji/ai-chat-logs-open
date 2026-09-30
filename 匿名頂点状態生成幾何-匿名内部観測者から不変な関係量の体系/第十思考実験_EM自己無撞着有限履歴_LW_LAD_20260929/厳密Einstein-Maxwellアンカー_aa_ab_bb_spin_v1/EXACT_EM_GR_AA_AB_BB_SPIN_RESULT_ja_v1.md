# 厳密 Einstein–Maxwell アンカーと aa / ab / bb スピン読出し検証 v1

## 目的
同一の状態生成系が、電荷なし・質量優勢・臨界・電荷優勢の領域を既知の厳密 Einstein–Maxwell 解と整合して扱えるかを確認し、さらに aa / ab / bb を一つの厳密解の内部成分として読み出せるか、各成分が spin-2 型位相変換を持ち得るかを検証した。

## 1. Reissner–Nordström 厳密族による比率スイープ
幾何単位 G=c=4πeps0=1 で

f(r)=1-2M/r+Q^2/r^2

を使い、M=1, r=20 で Q/M = 0, 0.1, 0.5, 1, 2, 10 を評価した。

- Q/M=0: Schwarzschild 極限
- 0<Q/M<1: subextremal RN
- Q/M=1: extremal RN
- Q/M>1: overextremal RN（数学的には厳密 Einstein–Maxwell 解だが、ブラックホールではなく naked singularity）

これは「電荷なし→質量優勢→臨界→電荷優勢」が同じ厳密族の中で連続に記述できることを確認するアンカーである。

## 2. 二中心 Majumdar–Papapetrou 厳密解と交差項
二中心では

U=1+u_a+u_b,

u_a=M_a/r_a,  u_b=M_b/r_b,

で、厳密 metric は

ds^2=-U^{-2}dt^2+U^2 d\mathbf{x}^2.

空間共形因子は厳密に

U^2 = 1 + 2u_a + 2u_b + u_a^2 + 2u_a u_b + u_b^2.

したがって

aa = u_a^2,

ab = 2u_a u_b,

bb = u_b^2

が一つの厳密解の中に同時に含まれる。

今回3点で数値評価し、aa+ab+bb を含む分解から U^2 を再構成した最大誤差は

4.44e-16

であった。

## 3. aa / ab / bb の spin-2 型読出し
複素状態

a=A_a exp(i phi_a),  b=A_b exp(i phi_b)

を考え、共通の横断面基底回転 psi に対して

a -> exp(i psi)a,

b -> exp(i psi)b

とする。このとき

aa=a^2 -> exp(2i psi) aa,

ab=2ab -> exp(2i psi) ab,

bb=b^2 -> exp(2i psi) bb.

さらに

(a+b)^2 = a^2 + 2ab + b^2

全体も同じ exp(2i psi) を持つ。

33個の回転角で数値検証した最大共変誤差は

1.35e-15

であった。

従って、a,b 自身が transverse-frame rotation に対して weight 1 の位相状態であるなら、aa, ab, bb は個別にも全体としても weight 2 を厳密に読み出せる。

ただし、これは「aa や bb が物理的 graviton である」ことの証明ではない。物理的 spin-2 の確認には、Einstein 方程式から得られる放射的 Weyl curvature / h+ , h× の変換則との一致が必要である。静的 Majumdar–Papapetrou 解そのものには重力波はない。

## 4. spin-1/2 型二重被覆
V=exp(i phi/2) とすると

V(phi+2pi)=-V(phi),
V(phi+4pi)=V(phi).

数値検証では

2pi 符号反転誤差 = 1.57e-16
4pi 復帰誤差 = 0

であった。

また V^4=exp(2i phi) なので、同じ乗法系列の中で half-angle carrier から weight 2 まで到達できる。

## 5. 現段階の結論

1. 電荷/質量比の全領域をテストする厳密アンカーとして RN 厳密族を使用できる。
2. aa / ab / bb を一つの厳密 Einstein–Maxwell 解の中に同時に含める要求は、少なくとも extremal 二中心 MP 解で成立する。
3. aa / ab / bb の spin-2 型位相読出しは数学的に可能で、3項は同じ weight-2 表現に属する。
4. ただし静的厳密解では「放射としての spin-2」は検証できない。次段階では、放射的な exact Einstein solution または full Einstein–Maxwell 数値解から Weyl scalar / h+,h× を読み、同じ weight-2 則と比較する必要がある。
5. 任意の二つの質量・電荷・距離を持つ静的二体系には Alekseev–Belinski の5パラメータ厳密 Einstein–Maxwell 解が存在し、今回の一般二体静的アンカーとして次に実装する候補である。
