# 論文7・現行seed図3対照テスト結果

## 結論

論文7の現行seed条件を変更せず、`N=5,40,300` を原本コードで
絶対step 0から55000まで再実行した。

五成分CSV、メタJSON、図3の個別PNG、およびN比較PNGは、
既存の論文7成果物とすべてバイト単位で一致した。

```text
all_passed = true
```

## 実行条件と結果

| N | PRNG seed | δ | crossing | 実行時間 |
|---:|---:|---:|---:|---:|
| 5 | 40265722 | `1e-15` | 1167 | 4.073 s |
| 40 | 40300722 | `1e-15` | 2011 | 19.520 s |
| 300 | 40560722 | `1e-15` | 4844 | 923.204 s |

全Nで五成分射影閉鎖誤差の最大値は `2.2e-16` だった。

## 照合結果

| 対象 | N=5 | N=40 | N=300 | N比較図 |
|:--|:--:|:--:|:--:|:--:|
| 五成分CSV・バイト一致 | ✓ | ✓ | ✓ | — |
| 五成分CSV・最大数値差0 | ✓ | ✓ | ✓ | — |
| メタJSON・バイト一致 | ✓ | ✓ | ✓ | — |
| 図3 PNG・バイト一致 | ✓ | ✓ | ✓ | ✓ |
| 図3 PNG・画素差0 | ✓ | ✓ | ✓ | ✓ |
| 図3 SVG・描画内容一致 | ✓ | ✓ | ✓ | ✓ |

SVGの生バイトには、生成日時とMatplotlibが実行ごとに付与する内部IDの差がある。
これらを正規化するとテキストは一致し、PNGは正規化なしでバイト単位・画素単位とも一致した。

## 検証ファイル

- `comparison/control_comparison.json`
- `source_verification.json`
- `logs/run_N00005.json`
- `logs/run_N00040.json`
- `logs/run_N00300.json`
- `outputs/raw/N00005/paper7_long_timeseries.csv`
- `outputs/raw/N00040/paper7_long_timeseries.csv`
- `outputs/raw/N00300/paper7_long_timeseries.csv`
- `outputs/figures/figure3_N00005_5color.png`
- `outputs/figures/figure3_N00040_5color.png`
- `outputs/figures/figure3_N00300_5color.png`
- `outputs/figures/figure3_compare_N5_N40_N300.png`

原本6ファイルのSHA-256は固定値と一致した。原本および既存成果物は変更していない。
