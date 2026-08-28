# N5_linear124_all3fix_control_rerun_20260828

`../N5_linear124_all3fix_seedless_parentnorm_removed_20260828/`（ChatGPT 追試 zip を展開したもの）の**対照実験**。参照プログラムを無変更で再実行して保存データと突合し、「make_parent の正規化除去」の意味を解析的に予言・検証し、参照の生成子が公開エンジンと同一であることを確認した。

## 結論（詳細は `対照実験結果_linear124_parentnorm_removed_20260828.md`）

- treatment（振幅込み K）は 4.8×10⁻⁸ で再現。baseline（位相のみ K、seedless）は成長率・等分配・σ₁=4 の不変量が再現し、onset と軌道は丸め誤差依存で再現しない（347 / 334 / 196）。
- 報告された parent residual 0.487 は σ·c·|1−c²|（c = |v| = 0.9272）の恒等式。正規化除去はスケール c と時間 c² の再尺度化に尽きる（H 系列が 1.5×10⁻⁸ で一致）。非正規化親は正規化親と頂点置換＋辺符号反転で厳密一致（重なり 1.000000000000）。
- 参照の K と公開エンジンの K は同一（2.7×10⁻¹⁶）。final σ₁ = 4 = N−1 は指数写像・L=124 でも成立。
- 振幅込み K は等分配せず σ₁=3.61（位相のみ K の等分配・σ=N−1 は旧プログラムの隠れた振幅正規化の性質。**追記参照：振幅込み K が理論に忠実な実装**）。

## 再現

```bash
cd N5_linear124_all3fix_control_rerun_20260828 && bash run_all.sh
```
参照パッケージ `../N5_linear124_all3fix_seedless_parentnorm_removed_20260828/data/` が同階層に必要。python3 + numpy + matplotlib、全体 30 秒程度。

## ファイル

- `program/`：参照プログラム 2 本（無変更）＋公開エンジン＋検証スクリプト 3 本
- `data/`, `figures/`：無変更再実行の出力
- `results/`：JSON 3 本＋ログ 5 本
- `run_all.sh`, `SHA256SUMS.txt`
