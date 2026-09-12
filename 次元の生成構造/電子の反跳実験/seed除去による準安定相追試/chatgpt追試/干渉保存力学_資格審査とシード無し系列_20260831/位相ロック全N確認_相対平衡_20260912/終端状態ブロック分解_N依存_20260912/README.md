# 追試B：ロック終端状態の最小ゼロ閉塞ブロック分解の N 依存（論文9 の提案検証）

日付: 2026-09-12。読出しのみ・物理無変更・新規走行なし。論文9（検討）が「終状態のブロック分解＝
種の目録が S_N 不変量に立つか（既存データで検定可能・未実施）」とした点を、既存 10000step 終端状態で実施。

## 実行
```
bash run_all.sh
```
`analyze_terminal_block_decomposition_20260912.py` が正本 `../../N3_N40_stage123_sweep_20260905/run_N3_N40_stage123_v1.py`
（SHA256 `1abf2353…`）から adjacency を import し、既存 10000step 終端状態
（`../../N3_N40_long10000_20260905/results/`）と床（`../../make_parent型初期値_自己無撞着構造_20260907/full_N3_N40_sweep/states/`）
のペアワイズ二波自己閉塞残差を N=3..40 で監査する。

## 結果（詳細は REPORT.md）
終端状態は全 N で厳密等モジュラー・Σz²≈0 保存だが exact ±i 対ゼロ。make_parent 床も exact ±i 対ゼロ
（論文4・実験15 と整合）。→ 最小ブロック分解は巡回Fourier/整数K 高対称床（N=8k、論文4）固有で、
データ系列の終端では立たない。

## 依存・入力
numpy のみ。入力は読み取りのみ（正本 adjacency＋既存 states npz）。出力は `results/`。
