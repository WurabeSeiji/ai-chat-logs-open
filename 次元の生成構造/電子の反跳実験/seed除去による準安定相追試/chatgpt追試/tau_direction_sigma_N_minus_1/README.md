# τ方向閉包の検証と σ = N−1 定理（再現パッケージ）

作成日 2026-08-27。対象メモ `../tau_direction_closure_and_UnI_consideration.md` の中心仮説（τ 方向閉包から U^n=I が創発するか）を、論文 v6（DOI 10.5281/zenodo.22112009）と同一のエンジン・同一の seed 規約で検査した記録。**新しい物理は追加していない**（エンジン無変更、γ = tan(π/144)、raw-K Cayley、seedless Z₀ = 親固有モード）。

## 結論（詳細は `検討結果_τ方向閉包とσ=N−1定理_20260827.md`）

1. N=5 の準安定相は剛体回転（相対位相一定）。step 単位の有限回帰は存在しない（メモの仮説を反証）。
2. 後期の集団回転数は **σ_eff = N−1**（N=4〜16 で σ_eff² = (N−1)² が 6 桁一致、seed 非依存）。N=3 は頂点閉包が不可能で σ² = 9/4 の別アトラクタ。
3. 定理：等モジュラー＋頂点 star 閉包 ⇒ 自分の生成子の固有モードで σ = N−1。数値検証済（Σ_{e'~e} e^{2iΔ} = −2、固有モード残差 ~1e-14）。
4. U^n=I の正体は「位相の有理ロック」ではなく「スペクトルの整数性」（自然時間で 2π 周期）。
5. Floquet：不安定モードは親にロック（角 0）、中立モード角は Δφ の整数倍——全 N で普遍。

## ファイル

| ファイル | 内容 |
|---|---|
| `検討結果_τ方向閉包とσ=N−1定理_20260827.md` | 分析結果・定理の証明・メモ各節への回答・次の作業 |
| `run_n_scaling_lowrank_v1_no_sigma_norm.py` | エンジン原本（`../complex_simplex_decompactification_N5_N16_20260826.zip` 収録と同一、SHA は `SHA256SUMS.txt`） |
| `analyze_n5_recurrence.py` → `results/n5_recurrence.json` | N=5 生データの剛体性・回帰・有理性検査、Floquet 角/Δφ |
| `sweep_sigma2.py` → `results/sweep_sigma2.json` | N=3〜16：親 σ²、後期 σ_eff²、剛体性、最終 K スペクトル、等分配 |
| `sweep_floquet.py` → `results/sweep_floquet.json` | N=3〜16：親相対平衡まわりの rotating-frame Jacobian（不安定・中立・減衰、角/Δφ、Σln|μ|） |
| `verify_theorem.py` → `results/verify_theorem.json` | seed 依存性、σ_eff²(t) 軌道、N=3 長時間、定理の中間量の検証 |
| `data/N5_phase_by_edge_5000steps.csv`, `data/N5_phase_increments_5000steps.csv` | `../N5_complex_simplex_complete_analysis_20260826.zip` からのコピー |
| `data/floquet_spectrum.csv` | `../N5_dynamics_followup_theorems_and_stability_20260826.zip` からのコピー |
| `results/*.log` | 各スクリプトの標準出力 |
| `run_all.sh` | 全再現（約 2 分） |
| `SHA256SUMS.txt` | 全ファイルのハッシュ |

## 再現

```bash
cd tau_direction_sigma_N_minus_1
bash run_all.sh          # python3 + numpy のみ（pandas 不要）
```

`sweep_sigma2.py` は引数で step 数と N リストを指定できる：`python3 sweep_sigma2.py 5000 3,4,5`。乱数 seed は decompactification パッケージと同じ `40260722 + 1000·N + seed`。

## 環境

macOS / Python 3.9 / numpy 2.0.2。倍精度。実行時間は全体で 2 分未満（N=16 の Floquet は 240×240 有限差分 Jacobian、数秒）。
