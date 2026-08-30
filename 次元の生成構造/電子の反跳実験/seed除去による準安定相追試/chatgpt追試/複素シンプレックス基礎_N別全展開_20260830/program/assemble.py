#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""断片 sections/N{5..16}.md を本文 md に合体し、N=5 補足・最終一覧表を付け、N=3/4 の図をインライン化する。"""
import os, re
import numpy as np
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
P=os.path.join(ROOT,'複素シンプレックス基礎_N別全展開_20260830.md')
s=open(P).read()
# 1) N=3 / N=4 の図をインライン化
s=s.replace('''- 図 1 `figures/N3_fig1_z_complex_plane.png`：3 本の波の複素平面（同一円上、0/60/120°）
- 図 2 `figures/N3_fig2_d2_closure_triangle.png`：d² ベクトルの首尾接続＝閉塞の三角形
- 図 3 `figures/N3_fig3_geometry_rank1.png`：埋め込み（1 複素軸内の正三角形）''',
'''図 1：3 本の波の複素平面（同一円上、0/60/120°）
![N=3 波の複素平面](figures/N3_fig1_z_complex_plane.png)
図 2：d² ベクトルの首尾接続＝閉塞の三角形
![N=3 閉塞の三角形](figures/N3_fig2_d2_closure_triangle.png)
図 3：埋め込み（1 複素軸内の正三角形）
![N=3 幾何](figures/N3_fig3_geometry_rank1.png)''')
s=s.replace('''- `figures/N4_fig1_z_complex_plane.png`：6 本の波（3 ペア × 2、同一円上 0/60/120°）
- `figures/N4_fig2_d2_closure.png`：6 本の d² ベクトルの首尾接続（各角度 2 本ずつで閉じる）
- `figures/N4_fig3_geometry_axes.png`：埋め込み（3 本の等長な複素軸に 4 頂点。丸い rank 3）''',
'''6 本の波（3 ペア × 2、同一円上 0/60/120°）
![N=4 波の複素平面](figures/N4_fig1_z_complex_plane.png)
6 本の d² ベクトルの首尾接続（各角度 2 本ずつで閉じる）
![N=4 閉塞](figures/N4_fig2_d2_closure.png)
埋め込み（3 本の等長な複素軸に 4 頂点。丸い rank 3）
![N=4 幾何](figures/N4_fig3_geometry_axes.png)''')
# 2) プレースホルダ以降を差し替え
marker='## 4. N=5（以下、同じ要領で N=16 まで順次追記）'
idx=s.index(marker)
s=s[:idx]
parts=[]
for n in range(5,17):
    frag=open(os.path.join(ROOT,'sections',f'N{n}.md')).read()
    if n==5:
        sup='''
### 4.6 補足：N=5 だけの特別な性質——実の擬ユークリッド空間に置ける

位相が 0° と 90° の 2 種類しかないため、d² = z² は +r²（辺の 5 組）と −r²（対角線の 5 組）の**実数**になる。距離の 2 乗が実数なら、複素空間を使わずに**実の空間で符号だけ混ぜた幾何**（擬ユークリッド空間）に置ける。実測（機械精度）：

- 埋め込み先は空間的方向 2 本＋時間的方向 2 本の 4 次元空間 R^{2,2}（B が実行列になり、固有値が +2 個・−2 個・0 が 1 個）。
- 5 頂点の座標は「2 枚の円」に載る：空間的 2 平面では角度 144°×j（対角線＝五角星の回り方）、時間的 2 平面では角度 72°×j（辺＝五角形の回り方）、**半径は同じ**。これはユークリッド 4 次元の正五胞体（4 次元の正単体）の標準座標そのもので、違いは片方の円が時間的であることだけ。
- 全 5 頂点で x·x = 0（光円錐上）。この光円錐は 4 次元空間の中の 3 次元の面で、その光的方向の全体は球面ではなく**トーラス**（2 つの円の積）。5 頂点は共通半径のトーラス上、時間円 1 周につき空間円 2 周の**巻き数 1:2 の 1 本の閉曲線**の上に等間隔で並ぶ。
- 辺の組（+r²）は空間的間隔、対角線の組（−r²）は時間的間隔——「対角線は長さでなく**間隔の種類**で辺と区別される」。

この読み替えができるのは（N=3〜16 の設計の中では）N=5 だけである【実測】。他の N は d² が本当に複素になり、複素シンプレックスが必須になる。
'''
        frag=frag+sup
    parts.append(frag)
# 3) 最終一覧表を再計算して生成
rows=[]
for N in range(3,17):
    E=[(i,j) for i in range(N) for j in range(i+1,N)]; M=len(E)
    if N==3:
        rows.append('| 3 | 3 | 3 | 3 | 0 | 1 | 3（0/60/120°） | 3 | **1**（例外） | —（rank1） |'); continue
    if N%2==0:
        n2=N-1; col={}
        for rr in range(n2):
            col[tuple(sorted((rr,N-1)))]=rr
            for k in range(1,N//2): col[tuple(sorted(((rr-k)%n2,(rr+k)%n2)))]=rr
        ncls=n2; step=180.0/(N-1)
    else:
        nn=(N-1)//2; col={}
        for d in range(1,nn+1):
            for i in range(N): col[tuple(sorted((i,(i+d)%N)))]=d-1
        ncls=nn; step=180.0/nn
    r=1/np.sqrt(15); th=np.array([np.radians(step*col[e]) for e in E]); z=r*np.exp(1j*th); d2=z*z
    D=np.zeros((N,N),complex)
    for val,(i,j) in zip(d2,E): D[i,j]=D[j,i]=val
    J=np.eye(N)-np.ones((N,N))/N; B=-0.5*J@D@J
    sv=np.linalg.svd(B,compute_uv=False); rank=int((sv>1e-12).sum()); ax=np.sqrt(sv[sv>1e-12])
    rnd='丸い' if (ax.max()-ax.min())<1e-8*ax.max() else '等長でない'
    rows.append(f'| {N} | {N} | {M} | {N} | {M-N} | {N*(N-1)*(N-2)//6} | {ncls}（k×{step:.6g}°） | {ncls} | {rank} = N−1 | {rnd} |')
summary='''## 16. 全 N の一覧表（まとめ）

各 N の節で個別に導いた数え上げと計算結果を、最後に一覧する（§1 の方針どおり、ここで初めて表にする）。

| N | 頂点 | 波 M | 辺 | 対角線 | 面 | 位相の種類 | d² 角度の種類 | rank | 軸の形 |
|---|---|---|---|---|---|---|---|---|---|
'''+chr(10).join(rows)+'''

注：
- rank は全 N で最大（N−1）に達する。例外は N=3 のみ（rank 1）。
- 軸がすべて等長（丸い）なのは N=4, 5, 7 の 3 つだけ。N≥9 は奇数の設計でも等長にならない。
- d² が実数（実の擬ユークリッド空間に置ける）のは N=5 のみ（§4.6）。
- μ/r² = −(N−1) は N=4〜16 のすべてで成立（N=3 は −3/2）。
'''
s=s+chr(10).join(parts)+chr(10)+summary
open(P,'w').write(s); print('assembled:', len(s), 'chars')
