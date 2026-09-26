# 再現条件

## 角度打切り
- $N=12,16,20,24$
- $p_0\in\{24.3,60,120\}$
- $e_0\in\{0.20,0.45,0.55\}$
- $q=1$, $\nu=0.25$
- 3 radial orbits
- 参照角度因子: $(1-2u)^{-1/2}$ をDOP853で同じ$p,e$発展上に積分。

## 積分次数
- 同じ5複素/正状態 $U,P,E,H,Q$。
- Taylor order $p=4,6,8$。
- 9ケース版CSVでは上記3x3全ケース。
- steps/orbit = 8,16,32,64。
- 参照: scipy DOP853, rtol=$2\times10^{-13}$, atol=$2\times10^{-15}$。
- `log` は最終$Q\to t$誤差読出しだけ。

## 図
`plot_order_convergence.py` はSVGだけを生成する。
