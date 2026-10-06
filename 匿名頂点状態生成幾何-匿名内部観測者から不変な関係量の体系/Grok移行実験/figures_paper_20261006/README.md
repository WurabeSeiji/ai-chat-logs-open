# figures_paper_20261006：論文用の図（4 系統）

4 系統の実験（水素原子 e–p、電子–電子、陽子–陽子、水素分子 H–H）の図を一か所で生成する。物理の計算は既存のプログラムをそのまま使い、本フォルダは描画と、描いた数値の保存（`data/*.json`）だけを行う。ラベルは英語（日英版の論文で共用）。

## 図の一覧

`figures_index.md`（生成時に書き出す）。15 枚：

| 系 | 図 | 元になる計算 |
|---|---|---|
| 水素原子 | fig01 準位図とカスケード、fig02 放出記録と線幅、fig03 時間発展 | `../exchange_engine_20261006/engine.py` ＋ `system_hydrogen.py`（固定一式 `../exchange_rel_20261005` に回帰済み：`regression_hydrogen.md`） |
| 水素原子・吸収体 | fig04 層 2（温度）、fig05 層 1（Doppler・反跳）、fig06 層 3（偏光・整列） | 同上。fig05 右の一原子の履歴は `../exchange_general_20261006/absorber_direction_recoil.py` の関数（乱数列 20261006、`absorber_direction_recoil.md` と同じ） |
| ee・pp | fig07 Mott 断面積、fig08 四重極制動放射 | `../exchange_general_20261006/scattering_same_sign.py` の `mott`・`quadrupole_spectrum`（import 時に同じ md が再生成される） |
| H₂ | fig09 準位図とカスケード、fig10 放出記録、fig11 検証、fig12 放射会合、fig13 層 1、fig14 層 3、fig15 重力の別勘定 | `../exchange_engine_20261006/{audit_h2_*.json, validation_h2.json, radiative_association_h2.json}`（`run_all.sh` の出力）。fig15 右は `system_h2.build()` の GW の辺 |

## 走り方

```
bash run_all.sh                 # engine の run_all（約 5 分、JSON を含む）→ make_figures.py（約 3 分）→ SHA256SUMS
SKIP_ENGINE=1 bash run_all.sh   # engine の出力が揃っているとき
```

注意：重い Python を同時に複数走らせない。Accelerate（BLAS）が全コアを使うので、並走させると 10 倍以上遅くなる（本日の作業で負荷平均 100 を記録）。

## 再現性

- 乱数列：fig05 左・中は `default_rng(7)`（回帰テストと同じ）、右は `default_rng(20261006)`（層 1 の md と同じ）。fig13 は `run_h2.py` の `default_rng(11)`。
- 図に描いた数値は `data/fig*.json` に保存。`SHA256SUMS` は make_figures.py・run_all.sh・README・figures_index.md・make_figures.log・figures/*.png・data/*.json。
- 依存：Python 3、numpy、scipy、matplotlib。フォントは DejaVu Sans（英語のみ）。
