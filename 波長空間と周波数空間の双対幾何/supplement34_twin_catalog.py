#!/usr/bin/env python3
# S=18 双子系 (9,9): 2016 配置の関係クラスカタログと分離スペクトル
import numpy as np, itertools
from collections import defaultdict

# shell-9 cells: (1,1,1,1)型 16 + (2,1,0,0)型 48
cells=[]
for k in itertools.product(range(-2,3),repeat=4):
    if sum((abs(t)+0.5)**2 for t in k)==9.0: cells.append(k)
cells=[np.array(k) for k in cells]
print(f"shell-9 セル数 = {len(cells)}")

# B4 群 (符号付き置換): 4!×2^4=384
perms=list(itertools.permutations(range(4)))
signs=list(itertools.product([1,-1],repeat=4))
def b4_orbit_key(va,vb):
    best=None
    for p in perms:
        for sg in signs:
            a=tuple(sg[i]*va[p[i]] for i in range(4))
            b=tuple(sg[i]*vb[p[i]] for i in range(4))
            for x,y in ((a,b),(b,a)):  # 双子は区別不能 → swap も同一視
                key=(x,y)
                if best is None or key<best: best=key
    return best

typ=lambda v: tuple(sorted(np.abs(v)))
classes=defaultdict(int)
info={}
pairs=list(itertools.combinations(range(len(cells)),2))
print(f"双子配置数 = {len(pairs)} (C(64,2)=2016 期待)")
for i,j in pairs:
    va,vb=cells[i],cells[j]
    d=va-vb; ssum=va+vb
    key=b4_orbit_key(tuple(va),tuple(vb))
    classes[key]+=1
    if key not in info:
        info[key]=(typ(va),typ(vb),int(d@d),int(ssum@ssum),int(va@vb))

# クラス集約: (型a,型b,|Δ|²,|Σ|²,内積) で表示
rows=defaultdict(int)
for key,cnt in classes.items():
    rows[info[key]]+=cnt
print(f"B4×swap 軌道クラス数 = {len(classes)} / 不変量タプルの種類 = {len(rows)}")
print()
print("関係クラス表 [(型a, 型b), |Δ|², |Σ|², <a,b>, 配置数, 比率%]  (|Δ|² 昇順)")
tot=sum(rows.values())
for (ta,tb,d2,s2,ip),cnt in sorted(rows.items(), key=lambda kv:(kv[0][2],kv[0][0],kv[0][1],kv[0][4])):
    print(f"  {ta}x{tb}  |Δ|²={d2:2d}  |Σ|²={s2:2d}  <a,b>={ip:+d}  W={cnt:4d}  ({100*cnt/tot:.1f}%)")
print()
# 分離スペクトル (|Δ|² の周辺分布)
dd=defaultdict(int)
for (ta,tb,d2,s2,ip),cnt in rows.items(): dd[d2]+=cnt
print("分離スペクトル |Δ|² : 配置数 (比率%)")
for d2 in sorted(dd):
    print(f"  |Δ|²={d2:2d} : {dd[d2]:4d}  ({100*dd[d2]/tot:.1f}%)   |Δ|={d2**0.5:.3f}  縞間隔=1/|Δ|={1/d2**0.5:.3f}")
print()
print(f"|Δ|²=1 の配置数 = {dd.get(1,0)}  (隣接禁止の確認)")
# 記録識別性: 不変量タプル (型対なし、|Δ|²,|Σ|² のみ=干渉縞で読める量) でクラスが分離するか
rec=defaultdict(list)
for (ta,tb,d2,s2,ip),cnt in rows.items(): rec[(d2,s2)].append(((ta,tb,ip),cnt))
amb=[k for k,v in rec.items() if len(v)>1]
print(f"縞データ (|Δ|²,|Σ|²) だけで分離できないクラス組 = {len(amb)} 件")
for k in amb[:10]:
    print(f"  (|Δ|²,|Σ|²)={k}: ", [t[0] for t in rec[k]])
