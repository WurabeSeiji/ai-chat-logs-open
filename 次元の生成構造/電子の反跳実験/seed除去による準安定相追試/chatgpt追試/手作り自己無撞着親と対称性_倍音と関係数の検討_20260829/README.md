# 手作り自己無撞着親と対称性_倍音と関係数の検討_20260829

検討途中の保存（2026-08-29）。既存データの読み取り＋手作り構成＋線形化＋短い直接走行のみ（長時間走行なし）。

- `検討記録_手作り親3分類_対称性とインフレーション_倍音と関係数_20260829.md` — 本文（定理・構成・観察・未決）
- `construct_and_analyze.py` — 和則検算（fixed_equimodular step 0, N=3〜16）、設計 A/B/星型の構成と残差、解析ヤコビアンによる共回転線形化（a/r²）、直接走行
- `run_all.sh` — 実行と SHA 更新
- `results/sum_rules_and_neighbor_angles.csv` — 和則誤差と隣接位相差の分類
- `results/stability_by_parent_type.csv` — 親の型ごとの残差・μ/r²・a/r²・不安定モード周波数
- `results/direct_runs.csv` — 直接走行の H⊥/H（τ 刻み）
- `results/run.log` — 実行ログ

入力：`../論文v1_全プログラム修正版_20260828/fixed_equimodular/`、`../N{5,6,7,8,10,16,20}_linear124_equimodular_selfconsistent_directHperp_treatment_only_2026082{8,9}/data/states_treatment.npz`
関連：`../公理見直し_ゼロ閉塞定理と固有時計_20260829/`、`../飽和ステップ数とNの関係_固定点ヤコビアン解析_20260829/`
