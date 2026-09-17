# 二波99走行 正本Dense再実験（2026-09-17）

準拠指示書（本フォルダ内）:

- `ClaudeCode_99run_theory_correct_dense_rerun_instruction_20260917_v3.md`（実験本体・v3 が確定版）
- `ClaudeCode_plotting_instruction_20260917.md`（図化）

旧99走行系列（二波サーベイ_生軌道データベース_20260915、Sparse/Lanczos 実装）は破棄。
比較基準でも正本でもない。本実験は理論式 K_ef = A_ef·Im(z̄_e z_f)（振幅込み・係数1）を
直接実装した新正本カーネルで、full dense 全固有分解（np.linalg.eigh）による
非近似時間発展をゼロから実行する。

## 状態

- Phase 1（コード作成・監査）完了。
- **Phase 2 実行中（2026-09-17）**: 実行ゲート撤去後、構造診断（記録のみ）→ 最小系列
  第1走行 L8_ma1_mb2_den8 → 3図生成。残り2本は木原確認後。
- 最小系列の全体像: L=8, den=8, 4096 steps の3条件
  `L8_ma1_mb2_den8`（奇偶）→ 確認 → `L8_ma1_mb3_den8`（奇奇）・`L8_ma2_mb4_den8`（偶偶）。
  99走行へは進まない。手順は `run_all_phase2.sh`。

## ファイル構成と出自

| ファイル | 内容 | 出自・変更 |
|---|---|---|
| `interaction_kernel_theory_v1.py` | 新正本カーネル（edges/adjacency/K_of/one_step） | v3 指示書 §2 コードの**文字単位一致コピー**（証跡 `provenance/kernel_vs_instruction_v3.diff` = 空） |
| `generator_twowave_v1.py` | 二波初期条件生成器 | `PhaseA_実装一式_v1_20260915/generator_phaseA.py` の物理コピー。変更は import 元1行（外部 import 禁止のため）＋docstring への出自追記のみ（証跡 `provenance/generator_vs_phaseA.diff`）。数式・辺順序（triu 列挙）は無変更 |
| `gen_manifest.py` | 99走行 manifest 定義（唯一の定義源） | 旧サーベイ正本と**バイト同一**（証跡 `provenance/gen_manifest_vs_survey.diff` = 空） |
| `structural_audit_v1.py` | §4 構造診断（新規） | 8項目を生値で記録するだけ（判定・閾値・実験開始ゲートなし。2026-09-17 指示でゲート撤去） |
| `run_dense_rerun.py` | ランナー（新規） | run_id 明示指定制。実行ゲートなし（2026-09-17 指示）。スレッド1固定・環境 fingerprint 保存。力学は kernel の one_step のみ |
| `run_all_phase2.sh` | Phase 2 手順書 | 未実行 |
| `plot_scripts_original/` | 旧図化正本4本（初期複素・終了複素・インフレーション・残差時系列） | 旧 figures_v1 と**バイト同一** |
| `plot_scripts_run/` | 図化実行用コピー | 出力パス文字列（各2行）＋**存在 run のみ図化ガード**（木原承認 2026-09-17、証跡 `plot_diffs/*.diff`） |
| `provenance/` | 出自証跡 diff 一式 | — |
| `SHA256SUMS_code.txt` | コード SHA256 固定（Phase 1 時点） | — |

## 力学（v3 §12 の二式のみ）

- K_ef = A_ef·Im(z̄_e z_f)（旧カーネルとの差: `np.angle` による振幅破棄を**しない**。旧は sin(θ_f−θ_e)）
- Z_{n+1} = exp(Δτ·K(Z_n))·Z_n を H=iK の全固有分解によるスペクトル指数写像で厳密評価
- Δτ = 2π/den。正規化・射影・補正・近似・低ランク化・疎化・打切りは挿入しない
- K^T=−K はコードで作らず、診断・invariants で観測する（実測 ~1e-16: numpy SIMD/FMA の丸め非対称による ulp 水準の残差。停止条件ではなく観測値）

## 保存形式（旧図化コードに合わせた上位互換）

- `runs/<run_id>/states.npz` — キー Z（4097×M complex128）, L, ma, mb, den, T（旧と同一）
- `runs/<run_id>/invariants.csv` — 列 `t,H,Q2_re,Q2_im,generator_norm,antisym_residual`
  （旧5列と同名・同順＋末尾に反対称残差1列追加。generator_norm は新Kでの ‖K(Z)Z‖/‖Z‖）
- `runs/<run_id>/metadata.json` / `SHA256SUMS.txt`
- run_id 形式 `L{L}_ma{ma}_mb{mb}_den{den}`（旧と同一）
- QA は機械的確認のみ（4097状態・NaN/Inf・metadata・SHA256）。H/Q2 drift・反対称残差は
  生値記録し閾値分類しない。FAIL でも states を削除しない。

## 図化の段階実行対応（木原承認 2026-09-17）

`plot_scripts_run/` の各本へ「**存在する states.npz の run だけ読込・図化し、
欠損パネルは空欄（axis off）にする**」最小変更を承認の上で追加した。

- 変更内容は各スクリプト共通: `have()` 存在チェック関数の追加、データ読込ループの
  存在ガード、パネルループの `if key in data / else axis('off')`、空データ時の即終了、
  `存在する run: n/99` の表示のみ。**パネル描画関数・軸・スケール式・色・凡例は無変更**。
- 全99走行が揃った時点では、出力はパス変更のみのコピーと同一になる（ガードが全て真）。
- 注意: 共通スケール S・重心 max|Σz|・Σz² の成分別 max は「読み込めた run の集合」上で
  計算される（式は無変更）。部分集合段階の図と99完備後の図でスケール値が変わり得る。
- 証跡: `plot_diffs/*.diff`（出力パス2行＋上記ガード変更。それ以外の差分なし）。
