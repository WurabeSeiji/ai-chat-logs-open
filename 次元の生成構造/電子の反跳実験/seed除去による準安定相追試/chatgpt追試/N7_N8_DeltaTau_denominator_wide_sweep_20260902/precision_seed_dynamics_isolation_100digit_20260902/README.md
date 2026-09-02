# 初期シード精度と時間発展精度の切り分け（64bit vs 100桁）

実行日: 2026-09-02
指図: `CLAUDE_CODE_EXPERIMENT_INSTRUCTION_precision_seed_dynamics_isolation_100digit_20260902.md`
（Drive file ID 1Y6WSnL29JHWXjpQteUsiX9XdOqyqHKws。実行時ローカル未同期のため
Google Drive API から全文取得して忠実実行した）

## 3条件

- **A: IC64 + Dynamics64** — 既存親 parent_v.npz、float64/complex128、numpy.linalg.eigh。
  既存 sweep（N8 D8 / N7 D7）との再現確認付き。
- **B: IC64 + Dynamics100** — A と数学的に同一の初期値（as_integer_ratio による
  exact binary64 lift、全成分検証済み）、力学のみ mpmath dps=100（eighe、π も100桁）。
- **C: IC100 + Dynamics100** — B と同一の100桁力学（実装は run_mp100_same_ic64.py から
  import、力学の同一性をコード共有で担保）、初期値のみ100桁解析生成。
  位相は正本生成器の設計規則（N=8: 1-factor、N=7: 距離クラス）を100桁で再現。
  振幅は legacy norm（解析式が存在しない凍結 binary64 定数）を exact lift して
  v·NORM/‖v‖ の規約を100桁で再現（指図 §4.2 の規則 1-3 に従う）。

走行順序: N=8 の A→B→C 完了・中間報告ののち N=7 replication（§16 遵守）。
全 run 2000 step で onset 到達、延長は不要だった。

## 入力

- `../data/N7/parent_v.npz` SHA256 `39f2a241…95a840b`
- `../data/N8/parent_v.npz` SHA256 `9412ec6d…8b691696`
（親 sweep の検証済みコピー。正本と同一 SHA256、読み取りのみ）

## 主結果（詳細は ANALYSIS.md）

- γ_τ が A/B/C で一致（N=8 で8桁、N=7 で6桁）— 指数不安定性は力学に内在。
- A≈B（onset +6 step のみ）— 発火源は IC64 に含まれる種、発展中の float64 注入は二次的。
- C は closure ~1e-102 の100桁初期値から同一 γ で約179桁を単一指数則で成長し発火
  （N=8: step 967、N=7: step 871）。onset 遅延は seed 桁数から定量的に予測可能。

## フォルダ構成

- `program/` — run_float64_control.py（A）/ run_mp100_same_ic64.py（B・共有力学）/
  build_ic100.py（IC100生成＋資格審査）/ run_mp100_ic100.py（C）/
  analyze_precision.py / plot_precision.py
- `data/N{8,7}_D{8,7}/{A,B,C}_*/` — timeseries.csv（100桁 run は110桁文字列保存）、
  run_info.json、checkpoint.json（mpf タプルによる正確な状態、再開可能）、
  B: lift_verification.json・precision_selftest.json、C: ic100_qualification.json ほか
- `results/` — precision_summary.csv（§7 の seed floor 記録を含む）/ growth_fits.csv /
  qualification_ic100.csv
- `figures/figP1..P7`

## 数値監査メモ

- 100桁 self-test: mpmath 1.3.0、1+1e-80≠1 / 1+1e-110==1、eighe 残差 ~2e-100。
- ノルム drift: float64 run ~1e-13、100桁 run ≤ 8.5e-98。補正なし。
- A の bitwise 非再現（相対差 ~1e-12、onset は完全一致）は Accelerate BLAS の
  プロセス間 ulp 非決定性として記録（ANALYSIS 問1・問10）。
- 100桁 CSV の極小値は固定小数表記になる場合がある（Decimal/float で解析可能）。
