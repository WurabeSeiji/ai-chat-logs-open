# Phase A 実装一式 v1 — Phase 0（力学カーネル監査）完了 2026-09-15

正本: 実装仕様書 v1.2（`../有限完全関係系_匿名因果写像_自己無撞着閉路_実装仕様書_v1.2_20260915.md`、ChatGPT 最終承認済み・変更禁止）。
実行範囲: 承認された実装順序 1〜8（Phase 0 まで）。**Phase 1/2/3 は未実行**（ChatGPT への Phase 0 結果再提示後に着手）。

## ファイル

- `kernel_canonical_dense.py` — 正本相互作用カーネルの忠実コピー（密参照）。出典: `run_N3_N40_stage123_v1_CANONICAL_REFERENCE.py` の edges/adjacency/H_of/one_step（14–28 行目）。物理式・数値手順は無変更。
- `kernel_sparse_matrixfree.py` — 疎/matrix-free 実装（仕様 §11）。$Ay=B(B^Ty)-2y$、$Kz=[\bar u\odot A(u\odot z)-u\odot A(\bar u\odot z)]/2i$、$\exp(\Delta\tau K)z$ は Hermitian Lanczos（tol=1e-14、完全再直交化）。
- `generator_phaseA.py` — 初期条件生成器（仕様 §4: 等振幅・半位相刻み）。力学へは匿名状態 $\{z_e\}$ のみを渡す。`analysis_basis` は解析専用（provenance）。
- `run_phase0_audit.py` — Phase 0 監査（検査 10 項）。
- `results/phase0_results.json` — 全数値。

## Phase 0 結果（L=6,12 × den∈{L（主）, 40（副）}、全合格）

| 検査 | 最大残差 |
|---|---|
| (1) A=0 ⇒ F=Id | 8.2e-17 |
| (2) 全位相同一 ⇒ F=Id | 1.6e-15 |
| (3) θ差∈{0,π} ⇒ F=Id | 1.4e-15 |
| (4) $S_P$ 置換同変（20置換×2状態×4系） | 8.6e-15 |
| (5) U(1) 同変（20位相） | 8.4e-15 |
| (6) 符号反転同変 | 5.9e-15 |
| (7) H 保存（500 step） | 4.3e-14 |
| (7) Q2 保存（500 step） | 7.6e-14 |
| (8) 生成器の事前登録値照合 | 差 ≤6.5e-7 = 仕様 §5 登録値の小数6桁丸め粒度内で一致 |
| (9) 密（正本）vs 疎（Lanczos） | step1 ≤4.6e-15 / step50 ≤1.9e-13 |

注記: 実行時 macOS Accelerate BLAS が複素 matmul で RuntimeWarning（divide by zero / overflow / invalid）を間欠的に出すが、全出力に NaN/Inf は皆無で最小再現でも消える spurious warning であることを確認済み（`results/phase0_results.json` の全値検査）。

## 確定した閾値（仕様 §10 の規則による。以後、本実験の結果を見て変更禁止）

実測 $\epsilon_{\rm floor}$ 構成: $\sqrt{M}\epsilon_{\rm mach}=1.8\times10^{-15}$、$\epsilon_{\rm expm}=10^{-14}$、$\epsilon_H=4.3\times10^{-14}$、$\epsilon_{Q_2}=7.6\times10^{-14}$、$\epsilon_{\rm equiv}=8.6\times10^{-15}$。

$$
\epsilon_{\rm floor}=7.609\times10^{-14}
$$

| 定数 | 採用値 | 結果 |
|---|---|---|
| $C_{\rm close}=10^2$ | $\varepsilon_{\rm close}=7.609\times10^{-12}$ | |
| $C_{\rm move}=10^6$（上限 $10^{-3}$） | $\varepsilon_{\rm move}=7.609\times10^{-8}$ | |
| $C_a=10^3$ | $a_{\min}/\|Z\|=7.609\times10^{-11}$ | |

**宣言: 上記 3 閾値は本コミットで確定し、以後 Phase 1〜4 の結果を見て変更しない。**

## Δτ 副系列の固定値（事前決定）

$$
\Delta\tau_0=\frac{2\pi}{40}\quad(\text{den}=40)
$$

根拠: 既存 N=40 正本走行（hm_N40_den_40）と同一刻み。旧規則 $2\pi/(\text{全頂点数})$ の歴史的値をそのまま継承した $L$ 非依存定数であり、閉路が出るように選んだ値ではない。

## Phase 1/2 結果（2026-09-15 追記。L∈{6,12}×den∈{L,40}×init∈{A,B,AB}=12 走行、T=500）

実装: `run_phase12.py`。出力: `results/phase12/`（時系列 CSV 12 本＋Δ_NL 4 本）、`results/phase12_summary.json`。
閾値は Phase 0 確定値を使用（変更なし）。dense_sparse_drift は指示どおり別監査量（閾値不算入）。

主要所見（数値的事実。解釈は ChatGPT 再提示後）:

1. **L=6 単一モード B（m=3）は力学的に厳密固定点**: 生成子ノルム ≤1.1e-13、D_orbit ≤2.8e-11。解析的理由: $U^{3\Delta}=(-1)^\Delta$ より状態は実×大域位相 → θ差∈{0,π} → K=0。$q=5$ で $\epsilon^{U(1)}=0$ が出るが $D_{\rm orbit}<\varepsilon_{\rm move}$ により非自明条件が正しく除外（自明解除外設計の実働例）。運動学予測（ord=2 の周期軌道）と力学（固定点）の差の最初の観測。
2. **L=12 単一モード B（m=3）は準回帰**: $\epsilon^{U(1)}_{\min}=6.3\times10^{-4}$@q=191（den=12）/ $5.4\times10^{-4}$@q=161（den=40）。閉路閾値には遠い。ただしこの 2 走行のみ dense_sparse_drift が 500 step で最大 0.373 に達し（軌道敏感性の累積）、回帰値の精密解釈は限定される。
3. 単一モード A（m=2）は全系で活発に動く（gen≈5〜11、D_orbit≈2）。単一モードでも力学は自明でない。
4. **AB 系の coverage は主刻み den=L で運動学 CRT 予測と数値一致**: L=6: 1.000（予測 1.000・CRT 独立）、L=12: 0.500（予測 0.500・CRT 非独立）。副刻み den=40 では 0.833 / 0.417 と予測から乖離。**ただし解釈上の重要な限定**: 全 AB 走行で $\eta_\perp$ 最大 0.75〜0.95（二波分解の大幅な崩れ、警告閾値 0.1 を大きく超過）、off-lattice residual も大きい（$r_B$ 平均 0.39〜0.77 rad、最大は格子間隔の半分＝セル境界レベル）。したがって coverage の数値一致を「離散 CRT 格子軌道の力学的保存」と読むことは現データでは支持されない（§8.3 の二段診断が設計どおり機能）。
5. **Δ_NL は即座に O(1)**: 最大 1.87〜2.19。二モード同時入力は単独モード和から強く逸脱し、非線形相互作用が支配的。
6. $\varepsilon_{\rm close}$ 級の非自明閉路（absolute / U(1) とも）は全 12 走行で未検出（$\epsilon^{U(1)}_{\min}$: AB 系で 0.43〜0.73）。
7. 保存量は全走行 $\Delta H/H\le5.9\times10^{-14}$、$\Delta Q_2\le2.7\times10^{-14}$（floor 内）。AB 系の dense_sparse_drift は ≤3.3e-8 と小さく、単一モード系（最大 0.373）との顕著な非対称がある。

状態: **Phase 3 主系列は未実行**（ChatGPT への Phase 1/2 再提示待ち）。

## 再実行

```sh
./run_all.sh          # Phase 0
python3 run_phase12.py  # Phase 1/2
```
