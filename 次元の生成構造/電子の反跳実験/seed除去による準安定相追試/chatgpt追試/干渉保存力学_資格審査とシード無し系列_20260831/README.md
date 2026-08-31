# 干渉保存力学：資格審査とシード無し第 1 系列（2026-08-31、やり直し版）

干渉を捨てない新フレーム（read: H_ef=A_ef·conj(z_e)z_f、step: exp(−iΔH) ユニタリ）の資格審査と、
**v2補完実験（../v2補完実験_4生成法_N3toN16_統一プロトコル_20260830/）と同一プロトコル**のシード無し走行。

## v2 との対応（変更は力学 1 点のみ）

- 生成法・rng・NORM 規約（各 N で mp 親のノルムに揃える）・受け入れ検査・判定規則（ρ−1>1e-3）・
  λ_f=2lnρ・t50・40000 step・L=124・種なし・列名・ファイル名・図の様式（fig1〜fig7）：**全て v2 と同一**
- 変更点（各ファイル冒頭に明記）：exp((2π/L)K) → exp(−i(2π/L)H)。sigma/ノルム/対称性検査列は H の同量。
  共回転角は実測位相進み −arg⟨v,Φ(v)⟩（旧版の解析式 (2π/L)μ と同値の手続き）
- 追加（置換ではない）：fig8（閉塞率グリッド、fig3 と同一様式。旧フレームで保存量だった |ΣZ²|/H が本フレームでは力学量のため）、
  pass7（星型アンカー検証）、新フレーム自己無撞着列（mu_new, residual_new_over_r2, is_equilibrium）

## 初版（2026-08-31 コミット e2225c01）の訂正記録

初版は許可なく (1) mp 系列を落とし 3 生成法 40 走行に縮小 (2) スケール規約を ‖v‖²=M/15 に変更
(3) 図を独自様式（閉塞率・重なり欠損・時計・最終状態の 4 枚）に差し替えており、v2 と比較不能だった。
本版で v2 同一プロトコル・同一図様式に全面やり直し。初版の資格審査（pass0）は規約非依存のため有効。

## ファイル

- 設計文書: 干渉保存力学_設計固定_20260831.md／結果: 実験結果_シード無し第1系列_20260831.md（やり直し後に更新）
- program/interference_dynamics.py … 統一相互作用と統一読出し（力学はこの関数のみ）
- program/pass0_qualification.py … 資格審査 V0.1〜V0.6
- program/pass1_parents.py, pass2_run.py, pass5_analysis.py, pass6_figures.py … v2 の同名/対応ファイルのコピー＋力学置換（冒頭に差分明記）
- program/pass2_embed_random.py, common.py, state_provider.py, original_engine.py … v2 からのコピー（無変更）
- program/pass7_final_structure.py … 本フレーム追加分析（星型アンカー検証）

再現: `bash run_all.sh`（pass0→pass1→埋め込み→走行54→pass5→pass6→pass7、約 40 分）。

注：data/<tag>/ の全 step 状態 states_treatment.npz と時系列 csv は容量のため git 管理外（run_all.sh で再生成可）。
親・summary・key_steps・results・図は管理内。
