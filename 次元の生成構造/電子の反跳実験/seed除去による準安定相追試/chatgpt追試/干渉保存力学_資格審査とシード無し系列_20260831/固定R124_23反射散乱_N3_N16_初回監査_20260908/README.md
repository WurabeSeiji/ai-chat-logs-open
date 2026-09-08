# 固定R124,23反射散乱 N=3..16 初回監査

- `R反射系_固定R124_23_段123合成_実験指示書_v1_20260908.md`: 実験仕様
- `run_R12423_scatter_N3_N16_v1.py`: 実行コード
- `results/summary_N3_N16.csv`: 系列別要約
- `results/timeseries_N3_N16.csv`: 全step指標
- `results/global_scatter_generator_spectrum.csv`: G_sc固有値と多重度
- `results/N*_control_states_500.npz`: control状態
- `results/N*_scatter_states_500.npz`: scatter候補状態
- `results/RUN_METADATA.json`: 実行条件
- `固定R124_23反射散乱_N3_N16_初回監査_分析_20260908.md`: 分析

注意: 今回のglobal scatter候補は square closure を一般に保存せず、物理実装としては不合格。負の結果を含む再現可能な監査として保存する。
