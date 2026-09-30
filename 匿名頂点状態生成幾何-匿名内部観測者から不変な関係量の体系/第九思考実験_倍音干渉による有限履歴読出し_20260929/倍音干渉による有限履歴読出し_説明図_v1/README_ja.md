# 倍音干渉による有限履歴読出し：説明図 v1

生成スクリプト: `generate_harmonic_history_figures_v1.py`

## 生成図

1. `figure01_monochromatic_vs_harmonic_time_resolution`
   - 単色波と倍音波の一周期内部時間分解能の比較。
2. `figure02_wave_period_fifo_phase_delay_memory`
   - 周期波を循環FIFO / phase-delay memory として読む模式図。
3. `figure03_frequency_comb_to_time_domain_resolution`
   - 周波数コムの帯域増加と時間領域局所構造の対応。
4. `figure04_odd_even_harmonic_roles_and_readout_gauge`
   - 奇数半整数系列と偶数系列の役割分解。

## 再現

```bash
python generate_harmonic_history_figures_v1.py
```

Python依存: numpy, matplotlib。
すべての図は同じスクリプトから PNG / SVG を同時生成する。

## 注意

図1・図3は解析式から生成した可視化であり、実測データではない。
図2・図4は概念構造を説明する模式図である。
