# 巡回二乗波の因数分解定理（検証パッケージ）

日付: 2026-09-13。読出しのみ・物理無変更・新規走行なし。論文4 §3–5 の一本化根拠。

## 目的

巡回Fourier床 $z_{ij}=r_d e^{i2\pi(i+j)/N}$ の零閉塞 $\sum z^2=0$ を、二乗波 $w=z^2$ の巡回構造から一本化する：

$$q=e^{i8\pi/N},\quad L=\operatorname{ord}(q)=N/\gcd(N,4),\quad p=\operatorname{spf}(L)\ (\text{最小ブロックサイズ}),\quad p=2\iff 8\mid N.$$

$N=8k$ の明示 $\pm i$ 構成（$m=N/8$、$z_{i+m}=\pm i z_i$、グラフ探索なしで完全 matching）と、全空間で「厳密2波±i対 $\iff 8\mid N$」（異距離クラス間の位相格子＋パリティ＋振幅反射 $r_d=r_{N/2-d}$）を確認する。

## 実行

```
bash run_all.sh
```

`verify_cyclic_square_wave_factorization_20260913.py`：距離クラス振幅 $r_d$ を committed の
`../高対称理論床_偶数系光子対精密解析_20260908/distance_class_amplitudes.csv` から読み、全ペア ±i 監査を
`even_N_summary.csv`（既往）と独立に照合。numpy 依存。終了コード 0 = 全 OK。

## 出力

- `results/cyclic_factorization_summary.csv`：全 $N$ の $L,p$・$\operatorname{ord}(q)$・p角形和・±i対数（本検証＝既往監査）・反射誤差・$N=8k$ 構成。
- `results/実行ログ_20260913.log`：数値表＋最小ブロック表。

## 結果（要約）

$N=3$–$40$：$\operatorname{ord}(q)=L=N/\gcd(N,4)$ 全一致、$p$ 角形和 $\le4.9\times10^{-15}$、$p=2\iff8\mid N$、
反射 $r_d=r_{N/2-d}$（誤差 $\le1.2\times10^{-16}$）、全ペア ±i 対数が既往監査と完全一致（$N=8k$：$84,424,1020,1872,2980$；他 0）、
$N=8k$ 明示構成が厳密 ±i かつ完全 matching（$14$–$390$ 対）。詳細は `分析_巡回二乗波因数分解定理_20260913.md`。
