# Stage A2d 保存済み完全無seed N=5方向系譜 報告書

## 実行状態

**SAMPLED_ROTATING_OR_MIXED_LINEAGE**

本解析はStage A2aで保存済みの完全無seed軌道だけを使用した。
新しい時間発展、状態更新、seed注入、横摂動、Benettin再投入は行っていない。

## 入力と標本

- 軌道: `N=5 float64 Z0=v; both explicit seeds OFF`
- A2a exec 1/2 bitwise一致: `True`
- 保存点数: `1024`
- step範囲: `0..5000`
- step間隔: `{'1': 11, '2': 12, '3': 13, '4': 9, '5': 978}`
- crossing: `1166`
- B0: `saved Bdom at step 0`

## 後期代表D34

- 後期窓: `1800..2500`
- 使用標本数: `141`
- 後期窓内overlap中央値: `0.944477747578`
- 後期窓内overlap最小値: `0.570122871215`

## 急拡大前の解像可能D34

- 解像基準: `min(q3,q4)/q1 >= 1e-06`
- crossing前の標本数: `40`
- step範囲: `[985, 1165]`
- 無seedD34対late overlap中央値: `0.142520216234`
- 無seedD34対late overlap範囲: `[0.12096442289609637, 0.14448924223390616]`
- 同じ保存stepでのseedありA2c overlap中央値:
  `0.149318721301`
- 無seed−seedあり overlap差の中央値:
  `-6.784815360611e-03`

急拡大前の解像可能部分空間は、急拡大後の後期代表部分空間と
単純に同一とは分類されない。完全無seed条件でも、既存seedありA2cと
同程度の低い早期対後期overlapが得られた。

## 無seed軌道とseedあり軌道の同一step直接比較

これは各軌道の後期代表方向へのoverlap比較ではなく、同一stepにおける
二つの `P34(t)` 自体のoverlapである。

| 区間 | 標本数 | overlap中央値 | overlap最小値 | 最大主角中央値(rad) | 最大主角最大値(rad) |
|---|---:|---:|---:|---:|---:|
| crossing前・解像可能 | 40 | 0.999999987161 | 0.999976415076 | 1.344618638131e-04 | 6.865869594358e-03 |
| crossing〜1799 | 128 | 0.999997466740 | 0.994305564668 | 2.127071312468e-03 | 1.069221907517e-01 |
| 1800〜2500 | 141 | 0.999852246337 | 0.990154654167 | 1.562753758376e-02 | 1.117166228234e-01 |
| 2500より後 | 500 | 0.265979686425 | 0.203572930166 | 1.329006688605e+00 | 1.561509505623e+00 |

- overlapが0.95を初めて下回る解像後step:
  `2690`
- overlapが0.5を初めて下回る解像後step:
  `3185`

急拡大前から三方向閉包の成立直後まで、無seedとseedありのD34部分空間は
ほぼ同一である。したがって、明示的初期seedは最初に生成される方向部分空間を
選択していない。step 2500以後の分岐は、生成方向の初期選択とは分離して
長期軌道差として扱う。

## 標本間回転

- crossing窓の標本対数: `106`
- 最大標本間主角: `1.190043454812` rad
- 最大標本間射影距離:
  `0.928387414621`
- `angle/delta_step` の最大記述値:
  `0.238008690962` rad/step

これらは不等間隔標本間の記述量であり、未保存stepを含む1-step最大値ではない。

## 数値健全性

- Bdom直交誤差最大: `9.734e-16`
- D34直交誤差最大: `9.942e-16`
- P34冪等誤差最大: `1.031e-15`
- P34対称誤差最大: `0.000e+00`
- 全配列finite: `True`

## 判定

分類: **SAMPLED_ROTATING_OR_MIXED_LINEAGE**

分類範囲: `sampled early-vs-late lineage only; no Tperp and no one-step maximum`

保存済み完全無seed軌道においても、急拡大前の解像可能D34は後期D34の
単純な微小振幅版ではない。観測された低overlapと標本間回転は、
方向部分空間が急拡大過程で回転・混合を伴って再編される記述を支持する。

## データだけでは言えないこと

- dominant planes are sampled mainly every 5 steps
- sample-to-sample maxima are not one-step maxima
- Tperp comparison at exact seedless t0 is not performed
- no single physical direction-establishment step is selected
- H1/H2/H0 is not judged
