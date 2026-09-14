# sens_v1 — SENS-2〜6 感度検査

5seed 主推定の定性的結論の頑健性検査（設計書 v1.2 §14、代表 N={6,10,16,20,30,40}）。

- run_sens.py — SENS-3/4/5/6（A/B × S0/S1/S2 × ε3値、K曲線、内蔵対照=主推定一致 assert）
- sens_generate_window_states.py — SENS-2 用窓状態（stage1 忠実コピー＋窓保存、bit一致検証内蔵）
- run_sens2_window.py — SENS-2（窓混入効果とスケール切替の判定）
- analysis_sens.md — 判定：rank-4 非支持は全感度軸に頑健。窓混入は上振れでなく測定対象の変更。
  副産物=軌道沿い局所次元 ≈4・N非依存（新仮説、要専用実験）

再現: `sh run_all.sh`（窓生成 約10分）。読出し専用・力学正本無変更。
