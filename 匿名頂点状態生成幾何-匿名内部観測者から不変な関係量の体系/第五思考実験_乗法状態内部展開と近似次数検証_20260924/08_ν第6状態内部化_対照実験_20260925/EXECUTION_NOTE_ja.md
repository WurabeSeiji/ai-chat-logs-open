# 実行記録

ChatGPT 実行環境の1回あたりツール時間制限のため、検証は同一コード関数を使って次の2群に分割して実行した。

1. 12ケース x 4000 macro steps: `run_nu_internal_state_control_v1.py` の `run_case()` を用いて時系列 raw CSV を生成。
2. 代表3ケース x 12000 macro steps および保存済み Paper-4 C1 CSV 照合: 同じモジュールの `run_case()` と `compare_saved_c1_csv()` を別プロセスで実行。

結果は `nu_internal_state_control_summary.json` に統合した。力学式、初期条件、状態遷移関数、比較式を途中で変更していない。分割は実行時間上の都合だけであり、物理条件や数値条件の変更ではない。

ローカル再現環境で時間制限がなければ `python run_nu_internal_state_control_v1.py` 一回で同じ検査を連続実行できる。
