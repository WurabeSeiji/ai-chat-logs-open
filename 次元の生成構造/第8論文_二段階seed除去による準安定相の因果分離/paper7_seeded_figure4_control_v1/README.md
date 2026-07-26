# 論文7・現行二段階seed図4対照再実行

論文7の原本プログラムを変更せず、現行の初期seedと図4用横摂動seedの
両方を有効にして、N=5, 40, 300 の図4を隔離再生成する。

## 固定条件

| N | 初期PRNG seed | 初期振幅 | 横摂動PRNG seed | 横方向数 | ε | 実行経路 | 記録間隔 |
|---:|---:|---:|---:|---:|---|---|---:|
| 5 | 40265722 | `1e-15` | 70005 | 3 | `1e-8,1e-10,1e-12,1e-14` | 非キャッシュ版 | 50 |
| 40 | 40300722 | `1e-15` | 70040 | 3 | `1e-8,1e-10,1e-12,1e-14` | キャッシュ版 | 50 |
| 300 | 40560722 | `1e-15` | 70300 | 3 | `1e-8,1e-10,1e-12` | キャッシュ版 | 100 |

- 横摂動開始: `t0 = crossing + 3000`
- Benettin再射影・再投入間隔: 500 step
- 絶対step上限: 55000
- 原本の出力先変数だけを `outputs/` へ差し替える。
- 原本ファイルと既存の論文7成果物は読み取り専用とし、上書きしない。

## 実行

```bash
python3 run_seeded_figure4_control_v1.py verify
python3 run_seeded_figure4_control_v1.py run 5
python3 run_seeded_figure4_control_v1.py run 40
python3 run_seeded_figure4_control_v1.py run 300
python3 run_seeded_figure4_control_v1.py figures
python3 run_seeded_figure4_control_v1.py compare
```

比較結果は `comparison/control_comparison.json` に保存する。

## 対照再実行結果

- N=5, 40, 300 の全計算を完走した。
- 数値CSVとメタJSONは既存論文7成果物とバイト単位で一致した。
- PNGは既存図とバイト単位・画素単位で一致した。
- SVGは生成日時とMatplotlib内部IDの正規化後に全文一致した。
- 総合判定: `all_passed=true`

詳細は `control_test_report.md` を参照。
