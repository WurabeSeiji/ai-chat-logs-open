# 高対称理論床・偶数系の倍音関係監査 v2

2026-09-08。**v2を再現正本とする。**

v1実行後、該当0件だった3つのCSVが0 byteとなりGoogle Driveが保存を拒否したため、v2では0件でもヘッダ行を必ず生成する。物理計算・判定値はv1と同一。

## 実行
```bash
./run_all_v2.sh
```

## 再現対象
- `analyze_highsym_harmonic_spectrum_v2.py`: 高対称床を定義式から再構成し、H=iKの正固有値を監査。
- `verify_highsym_harmonic_results_v1.py`: 全19偶数系列の陰性結果とスペクトル構造をassertで検証。
- `even_N_harmonic_summary.csv`: N別集計。
- `positive_frequency_modes.csv`: 正周波数と縮退度。
- `integer_harmonic_relations.csv`: 2,3,4倍音。0件なのでヘッダのみ。
- `low_denominator_rational_relations.csv`: q<=16非自明有理比。0件なのでヘッダのみ。
- `three_wave_sum_resonances.csv`: 三波和共鳴。0件なのでヘッダのみ。
- `RUN_METADATA.json`: 許容誤差・定義・限定。
- `高対称理論床_偶数系_倍音関係監査_20260908.md`: 分析本文。
- `SHA256SUMS_v2.txt`: v2正本ハッシュ。

## 主結果
全偶数N=4..40で2/3/4倍音、q<=16低分母有理比、三波和共鳴は0。N>=6の正固有周波数は縮退を除き3種類。N=8kの±i二波完全被覆は倍音出現と相関しない。占有床は単一carrierなので閉塞ブロック間は1:1共巻き。

## 限定
凍結生成子H=iKの監査であり、状態依存stage1+2+3写像のJacobian/Floquet監査は次段階。
