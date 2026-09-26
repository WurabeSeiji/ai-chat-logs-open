# 導出式: 方法Aと方法B

## A: 局所 exp/log

加算

$$
a+b
$$

を

$$
a+b=\lambda^{-1}\log\left(e^{\lambda a}e^{\lambda b}\right)
$$

で置換する。数値的同値性の比較用であり、`log` が内部演算に入るので最終物理モデルには採用しない。

## B: 内部状態展開

RK4を11相のmicrostate machineへ分解する。

0. $k_1$ 評価
1. stage2形成
2. $k_2$ 評価
3. stage3形成
4. $k_3$ 評価
5. stage4形成
6. $k_4$ 評価
7. 重み付き合成 $d$
8. 角度評価中点形成
9. $\Delta\phi$ 評価
10. macro commit

相番号そのものも one-hot 状態 $q\in\{e_0,\ldots,e_{10}\}$ とし、

$$
q_{m+1}=C_{11}q_m,\qquad C_{11}^{11}=I
$$

で循環させる。
