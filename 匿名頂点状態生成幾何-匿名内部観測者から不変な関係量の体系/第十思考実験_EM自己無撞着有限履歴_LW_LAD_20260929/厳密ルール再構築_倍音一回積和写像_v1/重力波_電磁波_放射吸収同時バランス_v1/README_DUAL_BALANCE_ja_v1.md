# 再現手順：重力波・電磁波の放射／吸収同時バランス

1. Python 3 と `numpy`, `scipy` を用意する。
2. このフォルダで `bash RUN_DUAL_BALANCE_ALL_v1.sh` を実行する。
3. `dual_q0p3_v2.json`, `dual_q0p6_v2.json`, `dual_q0p9_v2.json`, `dual_q0p99_v2.json` が生成される。
4. `dual_balance_repro_verification_v1.json` で独立再実行一致を確認する。

状態更新本体は `strict_harmonic_rn_tensor_v1.py` の `transition()` のみである。探索プログラムは初期状態を選ぶだけで transition を変更しない。

`search_dual_balance_strict_v1.py` と `dual_q0p9.json` は、同位置2チャネルという狭い初期状態族では Q/M=0.9 で十分な同時バランスに到達しなかったことを残す監査用パイロットである。最終結果は v2 を使用する。
