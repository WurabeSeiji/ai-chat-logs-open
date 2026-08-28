# N5_linear124_all3fix_control_rerun_40000_20260828

`../N5_linear124_all3fix_seedless_parentnorm_removed_40000_20260828/`（ChatGPT 追試 zip、40000 step）の対照実験。参照プログラムを無変更で再実行して突合し、旧プログラムの 4 つのミス（親の振幅正規化／初期 seed と正規化／Cayley 回転／位相のみ生成子）が修正版で正しく直っていることを静的・動的に検証した。

## 結論（詳細は `対照実験結果_linear124_40000step_20260828.md`）

- **4 修正すべて実装済み**（`results/verify_fixes.json`: all_four_fixes_ok = true）。線形回転は ‖log U − (2π/L)K‖ = 4×10⁻¹⁵、振幅込み K は K(2z) = 4K(z)。
- treatment は step 6834 まで 10⁻⁸ で一致、以後は丸め誤差の指数増幅で個別軌道が逸脱（統計的性質は一致：H⊥/H_total の 0.92〜0.94 への再上昇、PR/M ≈ 0.30、最小振幅 0.004）。baseline は step 256 から逸脱（不変量は一致）。
- 親は計算機ごとに対称性で同値な別解が選ばれる（振幅多重集合は 10⁻¹⁶ 一致）。

## 再現

```bash
cd N5_linear124_all3fix_control_rerun_40000_20260828 && bash run_all.sh
```
参照パッケージ `../N5_linear124_all3fix_seedless_parentnorm_removed_40000_20260828/data/` が必要。python3 + numpy + matplotlib、約 1 分。
