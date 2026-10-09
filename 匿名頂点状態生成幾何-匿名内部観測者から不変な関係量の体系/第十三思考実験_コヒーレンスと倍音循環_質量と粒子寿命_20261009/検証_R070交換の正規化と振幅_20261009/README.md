# 検証：R=0.70 の再帰交換における正規化と振幅（2026-10-09）

第十三思考実験の図 1（0713「交換干渉散乱行列フェルミオン的衝突における加速度基底と局在性交換予備実験総括 v1」の R=0.70 波形発展図）について、元プログラムの再帰更新

```
a_next = normalize(r * a + t * b)
b_next = normalize(t * a + r * b)
```

の各衝突での `normalize` が何をしているか、各波の振幅（ピークの高さ）がどう変わるかを確かめた。

## 手順

1. **忠実コピー**：`original_copy/run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1.py` は 0713 の原本と一字一句同一（`cmp` で確認）。原本が相対パスで読む 0711 の結果 JSON を、同じ相対位置 `20260711/…/…v2.json` にコピーし、コードを変えずに走らせた。
2. **対照**：コピーの出力を原本の出力（`波の情報読出し/20260713/exchange_scattering_matrix_fermionic_localization_transfer_preliminary_result_v1/`）と照合した。
3. **記録**：`measure_normalization_trace.py` は、コピーしたプログラムをモジュールとして読み込み、その関数だけを使って、原本の `recursive_snapshot_states` と同じ更新を衝突 0〜128 まで行い、記録用の行だけを足した。記録が原本の関数と同じ値を出すことを先に確かめた。
4. **比較**：同じ記録を、`normalize` を外した純粋な $U_R$ の更新でも取った。

## 結果

**対照（原本との一致）**

- CSV 6163 行、JSON の判定・全行が同一。
- このプログラムが出力する図 7 枚がピクセル単位で同一（R=0.70 波形発展図を含む）。
- 原本の結果フォルダにある `exchange_scattering_matrix_asymmetric_transfer_v1.png` は、このプログラムからは出力されない（別プログラムの出力）。
- 記録用ループの N_eff・L は、原本の `recursive_snapshot_states` の値と衝突 0,1,2,3,5,10,20,42 ですべて同一（`results/control_vs_original_snapshots.csv`）。

**正規化**

- A と B の重なり $|\langle a,b\rangle|$ の最大は $5.4\times10^{-15}$。A（q=+1、識別振動 m=1）と B（q=−1、m=2）は直交している。
- 正規化の前の各チャネルのノルム² $\|ra+tb\|^2$、$\|ta+rb\|^2$ は、全衝突で 1 からのずれが最大 $2.7\times10^{-15}$。
- したがって `normalize` は丸め誤差の水準でしか働いておらず、更新は純粋な $U_R$ と同じ。正規化あり・なしの N_eff の差は最大 $3.6\times10^{-14}$、L の差は $7.8\times10^{-18}$。純粋な $U_R$ でも各波のノルム² は 1 に保たれる（$0.9999999999999999$〜$1.0000000000000089$）。

**振幅（ピークの絶対高さ）**

$\sum_\eta|\psi|^2$ の最大値（総和で割らない、各波の全パワーは 1）：

| 衝突 | A のピーク | B のピーク |
|---:|---:|---:|
| 0 | 0.0039 | 0.1250 |
| 1 | 0.0402 | 0.0887 |
| 2 | 0.1056 | 0.0233 |
| 3 | 0.1216 | 0.0073 |
| 5 | 0.0109 | 0.1180 |
| 10 | 0.0304 | 0.0985 |
| 20 | 0.0867 | 0.0422 |
| 42 | 0.0648 | 0.0642 |

各波の全パワーは 1 のまま、局在が相手へ移ると、ピークの高さは約 0.004〜0.125 の間で約 30 倍変わる。原本の図は各線を最大値 1 に揃えて描いているため、この変化は図からは見えない。

**有効次数の閉形式**

衝突 0〜128 のすべてで

$$
N_{\rm eff}^A(k)=1+31\sin^2\!\frac{k\omega}{2},\qquad \omega=\pi+2\arcsin\sqrt R,
$$

との差は最大 $1.5\times10^{-12}$、$N_{\rm eff}^A+N_{\rm eff}^B=33$ からの差は最大 $7.1\times10^{-14}$。$N_{\rm eff}$ は倍音次数の重み付き平均で、A（次数 1）と B（奇数倍音 1〜63 の平均次数 32）の間で、交換された割合 $\sin^2(k\omega/2)$ だけ混ざる。和 33 は厳密に保存される。

## ファイル

- `original_copy/`：原本の忠実コピー、その出力（`exchange_scattering_matrix_fermionic_localization_transfer_preliminary_result_v1/`、検証メモ md、`control_stdout.txt`）
- `20260711/`：原本が読む 0711 の結果 JSON のコピー（相対パスを保つため）
- `measure_normalization_trace.py`：記録・対照・比較
- `results/`：`trace_renormalized_R070.csv`（原本の更新）、`trace_pure_unitary_R070.csv`（normalize なし）、`control_vs_original_snapshots.csv`、`summary.json`、`normalization_trace_R070.png`、`measure_stdout.txt`
- `run_all.sh`：一式の再実行（対照走行に約 5.5 分）
- `SHA256SUMS`

環境：Python 3.9.6、numpy 2.0.2、matplotlib 3.9.4。乱数は使わない。
