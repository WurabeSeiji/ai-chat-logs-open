# 再現条件

- `audit_two_paths_rk4.py`: A/Bを同一ファイルで定義した基底版。
- `audit_two_paths_rk4_fast.py`: one-hotで選択された相だけを実行する高速監査版。
- 代表3ケースを12000 macro steps。
- Paper 4 保存済みC1 CSVとの直接比較。
- $\tau_0=10^4$、$N=12$、$h=2\pi/4000$。

注意: 高速版の初期commit resetには同一microstepの新値参照が残っていた。これは次の `04_...厳密監査` で修正した。
