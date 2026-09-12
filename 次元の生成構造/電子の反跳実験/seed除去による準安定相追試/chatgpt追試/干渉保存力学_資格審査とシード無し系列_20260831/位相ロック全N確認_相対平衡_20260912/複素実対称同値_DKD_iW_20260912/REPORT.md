# 90°骨格：複素自己無撞着問題 ≡ 実対称0-1行列 W の固有値問題（$D^{-1}KD=iW$）と Perron–Frobenius

日付: 2026-09-12。読出しのみ・物理無変更・新規走行なし（既存床 npz、$N=3$–$40$）。論文1 §5 の核心の追加結果。

## 主張（導出）

$D=\operatorname{diag}(e^{i\varphi_e})$ とする。90°骨格では位相差 $\psi_{ef}=\varphi_f-\varphi_e\in\frac{\pi}{2}\mathbb Z$ なので $\sin2\psi_{ef}=0$、
したがって恒等式 $\sin\psi\,e^{i\psi}=\sin\psi\cos\psi+i\sin^2\psi=\tfrac12\sin2\psi+i\sin^2\psi=i\sin^2\psi$ が成り立つ。ゆえに
$$(D^{-1}KD)_{ef}=e^{-i\varphi_e}K_{ef}e^{i\varphi_f}=A_{ef}\sin\psi_{ef}\,e^{i\psi_{ef}}=i\,A_{ef}\sin^2\psi_{ef}=i\,W_{ef},\qquad W_{ef}=A_{ef}\sin^2\psi_{ef}.$$
すなわち $\boxed{D^{-1}KD=iW}$（$W$ は実対称・成分 $0/1$・非負）。よって床 $z=Dr$（$r$ 実正）について
$$Kz=i\omega z\iff D^{-1}KD\,r=i\omega r\iff iWr=i\omega r\iff \boxed{Wr=\omega r}.$$
**90°複素自己無撞着固有問題は実対称非負整数行列 $W$ の固有値問題と厳密同値**である。$W$ が既約（対応グラフが連結）なら
Perron–Frobenius により $\omega_{\max}=\rho(W)$ の主固有ベクトルは $r_e>0$・スケールを除き一意で、$z=Dr$ が90°骨格上の
自己無撞着単一回転 $Kz=i\rho(W)z$ をそのまま与える。

## 検証結果（`results/DKD_iW_summary.csv`、整数K床・make_parent 床、$N=3$–$40$ 各38床）

| 床 | max$\|D^{-1}KD-iW\|$ | $Wr=\omega r$ 一様性 | $\omega=\rho(W)=\sigma$(iK) | $r>0$ | Perron 重なり | $W$ 連結 |
|---|---|---|---|---|---|---|
| 整数K床 | $6.4\times10^{-14}$ | std $\le10^{-13}$ | 全一致 | 全床 | $=1.000000$ | 全床（成分1） |
| make_parent 床 | $1.9\times10^{-12}$ | std $\le10^{-12}$ | 全一致 | 全床 | $=1.000000$ | 全床（成分1） |

- $D^{-1}KD=iW$ が機械精度で成立（整数K床は厳密90°で $\sim10^{-14}$、make_parent 床は90°残差 $\sim10^{-12}$ 限界）。
- $Wr=\omega r$：$r=|z|$ で $\omega$ は全辺一様（std $\sim10^{-13}$）、$\omega=\rho(W)$（$W$ 最大固有値）$=\sigma$（iK スペクトル半径）が厳密一致。
- 床の振幅 $r$ は $W$ の**主固有ベクトル**（Perron 重なり $=1.000000$）で $r>0$。$W$ は全 $N=3$–$40$ で連結（成分1）＝Perron–Frobenius が全 $N$ で適用可。

## 判定（探究型物理学者ロール）

- 分類: 導出済み帰結（解析的同値 $D^{-1}KD=iW$＋Perron–Frobenius）＋全 $N$ 数値確認。判定: 保持（断定）。
- **「90°位相骨格の上に自己無撞着な正振幅波を作れるか」への数学的回答＝作れる**：$W$ が既約なら正の主固有ベクトルが一意に存在し、$z=Dr$ が自己無撞着床。本系では $N=3$–$40$ で $W$ は全て既約ゆえ、全 $N$ で構成可能。
- 論文1 はこの一般原理（複素自己無撞着 ≡ 実対称 $W$ の Perron 問題）を与え、論文3 は実際の高対称骨格で $W$ の形・$\rho(W)$ の閉形式（偶数 $N(N-2)$／奇数 $N^2-2N-1$）・振幅の2/3値縮退・整数 $K$ の $\sigma_{\max}$ との一致を解く。

## 未決 / 留保

- $W$ の既約性（連結性）は骨格のラベル付けに依存する一般的性質で、本系の全 $N=3$–$40$ で成立を確認したが、任意の90°ラベル配置での既約性の閉形式条件は本追試の範囲外。

## 成果物

- `verify_DKD_iW_equivalence_20260912.py`（正本 adjacency を SHA 照合の上 import・読出しのみ）
- `results/DKD_iW_summary.csv`, `results/実行ログ_20260912.log`
- `README.md`, `SHA256SUMS.txt`, `run_all.sh`
- 入力: 既存 整数K床（`../../90度理論床_整数K解析解_全N_20260909/parents_90deg_floor_analytic/`）と make_parent 床（`../../make_parent型初期値_自己無撞着構造_20260907/full_N3_N40_sweep/states/`）。依存: numpy, scipy。
