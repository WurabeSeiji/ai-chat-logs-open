#!/usr/bin/env python3
# 安定種 s=1,3,5 の代表周波数の調査
# (a) 符号語読み: 自分の容器 (周期 2√s) での論理波(箱和集合)の厳密スペクトル
# (b) アルファベット読み: 親格子の占有セル k の線 (補遺17)
import numpy as np, itertools, math

def centers4(s):
    R=math.sqrt(s); K=int(math.floor(R))+1
    rng=np.arange(-K,K+1)
    g=np.stack(np.meshgrid(rng,rng,rng,rng,indexing='ij'),axis=-1).reshape(-1,4)
    return g[((np.abs(g)+0.5)**2).sum(axis=1)<=s+1e-9]

for s in [1,3,5]:
    R=math.sqrt(s); P=2*R
    cells=centers4(s).astype(np.float64); N=len(cells)
    nc=6
    ax=np.arange(-nc,nc+1)
    n4=np.array(list(itertools.product(ax,repeat=4)),dtype=np.float64)
    ph=cells@n4.T
    S=np.exp(-2j*np.pi/P*ph).sum(axis=0)
    def phi(n):
        out=np.ones_like(n); nz=n!=0
        out[nz]=np.sin(np.pi*n[nz]/P)/(np.pi*n[nz]/P); return out
    PHI=phi(n4[:,0])*phi(n4[:,1])*phi(n4[:,2])*phi(n4[:,3])
    c=S*PHI/(P**4)
    mag=np.abs(c)
    dc_mask=(np.abs(n4).sum(axis=1)==0)
    mag_ac=mag.copy(); mag_ac[dc_mask]=0
    # 振動数 (セル単位 cycles/unit): nu = n/P
    nu=np.linalg.norm(n4,axis=1)/P
    # エネルギー重心
    E=mag_ac**2
    centroid=(E*nu).sum()/E.sum()
    # 上位5線 (縮退をまとめる: |n| 成分のソート形でグループ化)
    order=np.argsort(-mag_ac)
    seen={}
    for i in order[:4000]:
        if mag_ac[i]<1e-12: break
        key=tuple(sorted(np.abs(n4[i]).astype(int)))
        if key not in seen:
            seen[key]=(mag_ac[i], nu[i])
        if len(seen)>=5: break
    print(f"--- s={s} (R'=√s={R:.4f}, P=2√s={P:.4f}, セル数={N}) ---")
    print(f"  符号語(自容器)スペクトル 上位線 [形 (|n1..n4|), 振幅, ν=|n|/P (cycles/cell)]:")
    for key,(m,v) in sorted(seen.items(), key=lambda kv:-kv[1][0]):
        print(f"    n形={key}  振幅={m:.5f}  ν={v:.4f}")
    print(f"  基本波 ν0=1/(2√s) = {1/(2*R):.4f}   AC エネルギー重心 <ν> = {centroid:.4f}")
    print(f"  恒等式チェック: ν0 × R' = {1/(2*R)*R:.4f} (=1/2 なら零点)")
print()
print("アルファベット読み (親格子の占有セル線, 補遺17):")
print("  s=1: k=0          -> 線=DC,        裸の ν_obs=0,      零点ドレッシング: 0+0+1=1 ✓")
print("  s=3: k=±e_i       -> 線ノルム=1,   ν_obs=1,           1+1+1=3 ✓")
print("  s=5: k=(±1,±1,0,0)-> 線ノルム=√2,  ν_obs=1.4142,      2+2+1=5 ✓")
