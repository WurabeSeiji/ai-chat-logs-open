# 手作り親_設計点走行_N4toN8_20260829

手作り自己無撞着親（偶数=1-因子分解 A／奇数=ハミルトン閉路 C）N=4〜8 の種なし 40000 step 走行。変更点は親のみ（力学は N8 テンプレートと同一）。3 パス構成。

- `実験結果_手作り親N4toN8_偶奇の縞_20260829.md` — 結果（予測との照合、偶奇の縞）
- `program/pass1_make_parents.py` — 親生成＋受け入れ検査＋保存（決定論・乱数ゼロ）
- `program/pass2_run.py` — 走行（テンプレート同一力学、図なし）。`python3 pass2_run.py <N>`
- `program/pass3_figures.py` — 図化・局面 τ 判定・平面占有（読出し専用、再実行可）
- `program/original_engine.py` — テンプレートの byte コピー（participation_ratio 等の参照用）
- `data/N{4..8}/` — parent_v.{npz,csv}, parent_checks.json, states_treatment.npz(全 step), timeseries, key_steps, summary.json, plane_occupation.csv
- `figures/N{4..8}/` — 従来図＋平面占有＋5 局面×(複素平面図・幾何構造図)
- `results/` — stage_taus.csv, analysis_summary.csv, 各パス log
- `run_all.sh` — 全パス実行＋SHA 更新
