# 干渉保存力学：資格審査とシード無し第 1 系列（2026-08-31）

干渉を捨てない新フレーム（read: H_ef=A_ef·conj(z_e)z_f、step: exp(−iΔH) ユニタリ）の資格審査と、
シード無し走行 40 本（hm/ne N=3〜16、rb N=5〜16、40000 step、Δ=2π/124、‖v‖²=M/15）。

- 設計文書: ../干渉保存力学_設計固定_20260831.md（パッケージ直下）
- program/interference_dynamics.py … 統一相互作用（unified_interference_step）と統一読出し。力学はこの関数のみ
- program/common.py, state_provider.py … v2補完実験_4生成法_N3toN16_統一プロトコル_20260830/program/ のコピー（シリーズ内完結規約）
- program/pass0_qualification.py … 資格審査 V0.1〜V0.6（ユニタリ性・頂点形恒等式・M1/M2 凍結対照・現行フレーム対照・アンカー定理・閉塞不変）
- program/pass1_parents.py … 親 40 個の生成・受け入れ・走行前予測の固定（results/parents_predictions.csv）
- program/pass2_run.py <tag> … シード無し走行（timeseries.csv、snapshots.npz、final_state.npz、summary.json）
- program/pass3_analysis.py … 機械判定（results/matrix_N_by_method.{csv,md}）
- program/pass4_figures.py … figures/fig1〜fig4

再現: `bash run_all.sh`（約 15 分）。順序は pass0→pass1（予測固定）→pass2×40→pass3→pass4 で固定。

注：data/<tag>/timeseries.csv（全 step 時系列、計 307MB）は git 管理外。`bash run_all.sh` で再生成可。親・スナップショット・最終状態・summary・results は全て管理内。
