# 論文7・現行seed図3対照再実行

論文7の原本プログラムを変更せず、現行の初期seed条件で
`N=5, 40, 300` の五成分長時間時系列と図3を隔離再生成する。

## 固定条件

| N | PRNG seed | 初期seed振幅 |
|---:|---:|---:|
| 5 | 40265722 | `1e-15` |
| 40 | 40300722 | `1e-15` |
| 300 | 40560722 | `1e-15` |

- 初期状態: `Z0 = (v + 1e-15 g) / ||v + 1e-15 g||`
- 絶対step: `0..55000`
- サンプリング: N=5,40は25 step、N=300は100 step
- 原本の `run_paper7_5color_timeseries.py` と
  `make_paper7_figures.py::fig23()` を実行時に読み込む。
- 原本の出力先変数だけを `outputs/` へ差し替える。
- 原本ファイルと既存の論文7成果物は読み取り専用とし、上書きしない。

## 実行

```bash
python3 run_seeded_figure3_control_v1.py verify
python3 run_seeded_figure3_control_v1.py run 5
python3 run_seeded_figure3_control_v1.py run 40
python3 run_seeded_figure3_control_v1.py run 300
python3 run_seeded_figure3_control_v1.py figures
python3 run_seeded_figure3_control_v1.py compare
```

比較結果は `comparison/control_comparison.json` に保存する。
