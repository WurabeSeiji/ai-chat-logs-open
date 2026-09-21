#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Analytic/numeric summary for E3-B0 symmetric three-body rigid control."""
from pathlib import Path
import csv, json, math

HERE=Path(__file__).resolve().parent
SUMMARY=HERE/'E3B0_symmetric_three_body_summary_v1.csv'
RESULTS=HERE/'E3B0_symmetric_three_body_results_v1.json'
REPORT=HERE/'E3B1_symmetric_three_body_analysis_v1.md'

with open(SUMMARY,encoding='utf-8') as f:
    rows=list(csv.DictReader(f))
with open(RESULTS,encoding='utf-8') as f:
    results=json.load(f)

nz=[r for r in rows if float(r['eps'])>0]
max_mass=max(float(r['mass_ratio_max_dev']) for r in nz)
max_cent=max(float(r['centroid_rel_max']) for r in rows)
max_pd=max(max(float(r['pair_distance_equality_error']),float(r['pair_distance_drift'])) for r in rows)
max_area=max(float(r['triangle_area_rel_drift']) for r in rows)
max_readout=max(max(float(r['readout_pair_distance_equality_error']),float(r['readout_pair_distance_drift'])) for r in rows)
max_cfit=max(float(r['directed_coupling_rel_error']) for r in nz)
max_bodysplit=max(max(float(r['body_radius_split']),float(r['body_q_split'])) for r in rows)

text=rf'''# E3-B1 完全対称三体・剛体対照実験 解析 v1

日付: 2026-09-21

## 1. 問い

三つの匿名状態 $a,b,c$ を完全対称に置き、完全グラフ $K_3$ 上で同じ結合を与えたとき、
三体になったことだけから個体差（質量比的・電荷比的読み）が自発的に生じるかを調べる。

候補則は

$$
X_{{i,k+1}}+X_{{i,k-1}}=\tau X_{{i,k}}+\varepsilon\sum_{{j\ne i}}X_{{j,k}}.
$$

初期状態は原点を重心とする正三角形であり、三状態は完全に置換対称である。

## 2. 解析的縮約

完全対称初期条件では

$$
X_a+X_b+X_c=0.
$$

したがって各 $i$ について

$$
\sum_{{j\ne i}}X_j=-X_i,
$$

ゆえに

$$
X_{{i,k+1}}+X_{{i,k-1}}=(\tau-\varepsilon)X_{{i,k}}.
$$

つまり完全対称・重心ゼロ部分空間では、三体完全グラフ則は厳密に
**同一の二値二階則を三つの位相だけずらして走らせる剛体回転**へ縮約する。

今回、

$$
\tau-\varepsilon=2\cos(2\pi/31)
$$

と固定したので、全 $\varepsilon$ で相対モードは同じ31歩周期を持つ。

## 3. 数値結果

走査: $\varepsilon=0,10^{{-5}},10^{{-4}},10^{{-3}},5\times10^{{-3}},10^{{-2}}$、各3100 step。

- 重心閉包の最大相対誤差: `{max_cent:.3e}`
- 三辺長の等値・剛体性の最大誤差: `{max_pd:.3e}`
- 三角形面積の最大相対ドリフト: `{max_area:.3e}`
- 二乗読み出し後の三辺長の最大誤差: `{max_readout:.3e}`
- 6方向結合の回帰値と入力 $\varepsilon$ の最大相対誤差: `{max_cfit:.3e}`
- 質量比的読み $C_{{i\leftarrow j}}/C_{{j\leftarrow i}}$ の 1 からの最大偏差: `{max_mass:.3e}`
- 状態ごとの半径・一歩面積の最大分裂: `{max_bodysplit:.3e}`

事前宣言した E3B1--E3B7 はすべて PASS した。

## 4. 質量比的読み

全ての方向付き結合は

$$
C_{{a\leftarrow b}}=C_{{b\leftarrow a}}=
C_{{b\leftarrow c}}=C_{{c\leftarrow b}}=
C_{{c\leftarrow a}}=C_{{a\leftarrow c}}=\varepsilon
$$

に数値精度内で一致した。したがって

$$
\boxed{{m_a:m_b:m_c=1:1:1}}
$$

以外の質量的個体差は読めない。

これは「読み出しに失敗した」のではなく、完全対称性から非対称な個体属性を作らないという対照実験の期待結果である。

## 5. 電荷比的読み

質量的非対称成分を除いても、三つの対に残る対称結合は全て同じである。
したがってこの実験だけからは三者を区別する電荷比は作れない。
もし対称残差を $P_{{ij}}\propto \sigma_i\sigma_j$ と因数分解するなら、読めるのは

$$
|\sigma_a|=|\sigma_b|=|\sigma_c|
$$

という等値までで、全体符号や共通オフセットは決まらない。

## 6. 結論

$$
\boxed{{
\text{{完全対称三体系では、三体になっただけでは個体差は生成されない。}}
}}
$$

さらに強く、この候補では重心ゼロ条件のため三体則自身が同一二体則へ厳密縮約する。
したがって、質量比・電荷比の非自明な読み出しには、**無名性を破らずに生じる非対称性**が必要である。

次の最小実験は、外部から粒子名を与えるのではなく、初期関係または保存量に最小の一パラメータ非対称を入れ、
その非対称が方向付き応答比へどう写るかを調べることである。
'''
REPORT.write_text(text,encoding='utf-8')
print(REPORT)
