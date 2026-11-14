# 補遺 A・B 数値実験の検証（2026-11-14）

- 元プログラム `../sim_symmetric_crossing.py`・`../sim_asymmetric_transmission.py` を無改変でコピーして再実行。
- 図：赤い曲線の画素は元の PNG と完全一致（対称版 10082 画素、非対称版 7341 画素、不一致 0）。差は文字のフォントのみ。
- `check_values.py` → `check_values_result.txt`：ラベルの値、w と −w の対称性、36° の倍数以外の位置での山の高さ。
- `check_compartment_readout.py` → `check_compartment_readout_result.txt`：「観測者の分解能はゲージである」（Concept 10.5281/zenodo.24097611）の (T5)・図 3 の読みで、ω を 3° 刻みで連続に動かし、最大の山の位置と区画の歯の高さ・比を読む。
