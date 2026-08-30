# 複素シンプレックス_重心閉塞_非等モジュラー族_20260830

等モジュラー条件を外した複素シンプレックス（自己無撞着状態）の組み直し。基礎文書 `../複素シンプレックス基礎_N別全展開_20260830/` の位相配置・K・B・Takagi を一切変えず、振幅だけを自由にした。

**主結果**：(1) 自己無撞着は各辺で S_i+S_j=(2μ+W_i+W_j)z_ij²/|z_ij|² と同値（厳密恒等式）。局所閉塞＋頂点重み均等 ⇒ 自己無撞着（μ=−W）。(2) 等モジュラー点を通る自己無撞着解集合は N²−4N+1 次元の多様体（N≥5、ヤコビアン核次元＝{S=0, W=const} の拘束数 3N−1 から出る次元と一致、Newton で 100/100 到達）。クラス重み付き族（Σa_cω^c=0、次元 q−3）はその小さな部分族。(3) 力学：族メンバー（N=6,8）は全て等モジュラー親と同じくインフレーション、λ は共回転モノドロミー予測と 4 桁一致。対称性なしの乱数均衡親 20 個（N=5〜8）も全てインフレーション（予測比 0.98〜1.02）。**床（中立）は対称親（N=4,5,7、丸い）だけの非一般的性質**。(4) N=9 は対称親・族メンバーとも床（事前予測の閾値誤り→較正で修正、§7.3）。

- `複素シンプレックス_重心閉塞_非等モジュラー族_20260830.md` / `.pdf` — 本文（導出全展開・仮説→反証→修正の記録）
- `program/common.py` — 辺・隣接・位相クラス・K・自己無撞着検査・B/Takagi（基礎文書と同一規約）
- `program/pass1_family.py` — クラス重み付き族の実測 N=3〜16（`results/family_table.csv`, `family_members.csv`, `data/family_path_N*.csv`）
- `program/pass2_jacobian.py` — 自己無撞着写像のヤコビアン核次元（`results/jacobian_nullity.csv`）
- `program/pass2b_nullspace_probe.py` — 核方向摂動→Newton 精密化（`results/nullspace_probe.csv`, `nullspace_probe_pca_N*.csv`）
- `program/pass3_parents.py` — 走行用親 9 個＋走行前予測（`data/N*_eps*_k*/`, `results/parents_predictions.csv`）
- `program/pass3b_monodromy_calibration.py` — ρ の較正（走行後追加、`results/monodromy_calibration.csv`）
- `program/pass7_balanced_random.py` — 恒等式検証＋乱数均衡親 20 個＋予測（`results/identity_check.csv`, `balanced_random_parents.csv`, `data/random_N*_s*/`）
- `program/run_dynamics.py` — 走行（手作り親パッケージ pass2_run.py のコピー、変更 3 行）、`program/original_engine.py` — テンプレート byte コピー
- `program/pass5_analysis.py` — 集計（`results/dynamics_summary.csv`）、`program/pass6_figures.py` — 図 11 枚
- `run_all.sh` — 全再実行＋PDF＋SHA、`build_pdf.sh` — PDF のみ

**データの所在**：`data/<親>/` に parent_v.{npz,csv}・parent_checks.json・summary.json・key_steps.csv・走行 CSV（gzip、各 4 MB）・全 step 状態 `states_treatment.npz`（各 8〜23 MB、計 323 MB）。**states_treatment.npz は git に入れていない**（`run_all.sh` で 12 分で再生成、Drive 上の正本には残す）。図 7 だけがこれを読む。
