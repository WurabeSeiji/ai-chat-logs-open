# 旧実装の位置づけ変更通知

2026-09-29 再監査により、以下の旧実装は標準理論参照計算としては有効だが、`STRICT_EXPERIMENT_RULES_ja` に従う universal-interaction state generator としては不適合であることを明示する。

- `倍音主状態_配列除去_v3/em_retarded_harmonic_lad_v3.py`
- `RWZ_spin2波方程式対応_v1/.../rwz_spin2_flux_v1.py`
- `RN_coupled_spin2_spin1_放射対応_v1/rn_coupled_gem_v1.py`
- `full_nonlinear_Einstein_Maxwell_Kastor_Traschen_v1/kastor_traschen_two_center_v1.py`

理由：RK4、leapfrog previous/current persistent arrays、または external exact-state sampling が strict rules R2/R4/R6 の少なくとも一部に適合しないため。

新しい strict 実装は本フォルダの `strict_harmonic_*_tensor_v1.py` を正本とする。
