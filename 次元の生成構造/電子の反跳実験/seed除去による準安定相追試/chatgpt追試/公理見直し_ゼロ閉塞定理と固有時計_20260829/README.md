# 公理見直し_ゼロ閉塞定理と固有時計_20260829

シード無し論文 v1 の版上げに向けた公理構成の見直し。新走行なし（read-only 再集計）。

- `公理見直し_ゼロ閉塞は定理_τは時間でない_20260829.md` — 本文（定理・証明・実測・公理構成案）
- `analyze_axioms.py` — 既存パッケージの step 0 と時系列から全数値を再生成
- `run_all.sh` — 実行と SHA 更新
- `results/closure_step0.csv` — 4 系統 × N=3〜16 の Σa², Σb², Σab, |Σz²|/H
- `results/closure_conservation_and_phase_advance.csv` — directHperp N=5,8,10,16,20 の |ZᵀZ|, σ₁, 位相進み（要所 step）
- `results/phase_advance_summary.csv` — 同・要約
- `results/sigma_spectrum_equimodular_parent.csv` — 等モジュラー親の σ/r² と有理近似
- `results/analyze_axioms.log` — 実行ログ

入力（相対パス）：`../論文v1_全再現テスト_20260828/original/`、`../論文v1_全プログラム修正版_20260828/{fixed_baseline,fixed,fixed_equimodular}/`、`../N{5,8,10,16,20}_linear124_equimodular_selfconsistent_directHperp_treatment_only_20260828/data/`
