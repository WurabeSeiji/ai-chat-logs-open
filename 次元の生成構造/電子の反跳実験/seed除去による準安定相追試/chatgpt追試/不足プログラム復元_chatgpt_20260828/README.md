# 不足プログラム復元_chatgpt_20260828

ChatGPT が Drive に置いた復元パッケージ（B-1〜B-7）の取得・展開・検証。結果は `検証結果_ChatGPT復元プログラム_20260828.md`。

```bash
cp -R ../論文v1_全再現テスト_20260828/original verify_root; mkdir -p compat/bits; cp ../論文v1_全再現テスト_20260828/compat_bits_stdc++.h compat/bits/stdc++.h
export CPLUS_INCLUDE_PATH=$PWD/compat; python3 extracted/論文v1_全プログラム修正版_20260828/run_all_reconstructed.py --root $PWD/verify_root --decompact-results $PWD/verify_root/complex_simplex_decompactification_N5_N16_20260826/results --n5-raw-k-source $PWD/verify_root/K_sigma_normalization_artifact_test_N4_N5_20260826/N5_raw_K_raw_observables.csv; python3 compare_reconstructed.py
```
