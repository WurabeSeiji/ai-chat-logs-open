#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺73: 補遺72 §4 引き継ぎの実行 — 独立再検証・Δ(13)・graded変種・完全結合振幅版・零チャネル検査
import itertools, time
from fractions import Fraction as F
from math import comb
from collections import Counter, defaultdict
src=open('supplement57_convention_dependence.py',encoding='utf-8').read()
exec(src.split('# 現行(字句順序)規約')[0])
CAP={1:1,3:8,5:24,7:40,9:64,11:96,13:96}
def M_channel(f,used):
    tot=1
    for m,k in Counter(f).items():
        free=CAP[m]-used.get(m,0)
        if free<k: return 0
        tot*=comb(free,k)
    return tot
def seq_twin(s):
    fins=all_finals(s)
    M1={f:M_channel(f,{}) for f in fins}; T1=sum(M1.values())
    P=defaultdict(lambda: F(0))
    for f1 in fins:
        if M1[f1]==0: continue
        used=Counter(f1)
        M2={f:M_channel(f,used) for f in fins}; T2=sum(M2.values())
        for f2 in fins:
            if M2[f2]==0: continue
            P[tuple(sorted([f1,f2]))]+=F(M1[f1],T1)*F(M2[f2],T2)
    return dict(P)
def batch_twin(s):
    fins=all_finals(s); W=defaultdict(int)
    for i,f1 in enumerate(fins):
        for f2 in fins[i:]:
            joint=Counter(f1)+Counter(f2); tot=1
            for m,k in joint.items():
                if CAP[m]<k: tot=0; break
                tot*=comb(CAP[m],k)
            if tot>0: W[tuple(sorted([f1,f2]))]+=tot
    T=sum(W.values()); return {k:F(v,T) for k,v in W.items()}
def derived_twin(s,Wd,graded=False):
    fins=all_finals(s); P={}; tot=0.0
    for i,f1 in enumerate(fins):
        for f2 in fins[i:]:
            joint=Counter(f1)+Counter(f2)
            if any(CAP[m]<k for m,k in joint.items()): continue
            g=1.0
            if graded:
                Mj=1
                for m,k in joint.items(): Mj*=comb(CAP[m],k)
                Mi=M_channel(f1,{})*M_channel(f2,{})
                g=Mj/Mi if Mi>0 else 0.0
            w=(2 if f1!=f2 else 1)*Wd.get(f1,0)*Wd.get(f2,0)*g
            if w>0: P[tuple(sorted([f1,f2]))]=w; tot+=w
    return {k:v/tot for k,v in P.items()}
def TV(P,Q):
    ks=set(P)|set(Q)
    return sum(abs(float(P.get(k,0))-float(Q.get(k,0))) for k in ks)/2
CH2={}
def ch2(vx,c,d):
    key=(vx,c,d)
    if key in CH2: return CH2[key]
    tot=0+0j
    if c==d:
        sh=SH[c]
        for ic in range(len(sh)):
            for idd in range(ic+1,len(sh)):
                e2=tphase(vx,sh[ic],sh[idd])
                if e2 is not None: tot+=e2
    else:
        for vc in SH[c]:
            for vd in SH[d]:
                e2=tphase(vx,vc,vd)
                if e2 is not None: tot+=e2
    CH2[key]=tot; return tot
def runW(vp,s):
    fins=all_finals(s); z={f:0+0j for f in fins}
    for f in fins:
        for ((a,b),x,(c,d),d1) in two_step_paths(s,f):
            for ia,va in enumerate(SH[a]):
                for ib,vb in enumerate(SH[b]):
                    if a==b and not ia<ib: continue
                    e1=tphase(vp,va,vb)
                    if e1 is None: continue
                    for vx in ([va,vb] if a==b else [va if x==a else vb]):
                        z[f]+=e1*ch2(vx,c,d)
    return {f:abs(z[f])**2 for f in fins}

print("==== V: claude.ai 結果の独立再検証 ====")
W9=runW((2,1,0,0),9)
print(f"s=9 W: { {f:round(w) for f,w in W9.items() if w>1e-9} }")
seq9=seq_twin(9); bat9=batch_twin(9); der9=derived_twin(9,W9)
X9=tuple(sorted([(1,3,5),(3,3,3)]))
PX_seq=float(seq9[X9]); PX_bat=float(bat9[X9])
print(f"逐次 P(X)={PX_seq:.5f} (396/403={396/403:.5f}), 一括={PX_bat:.5f} (60/61={60/61:.5f})")
print(f"導出 P(X)={der9[X9]:.5f} (0.98195), Δ(9)=TV={TV(der9,seq9):.5f} (0.00068)")
M531=M_channel((1,3,5),{}); M333=M_channel((3,3,3),{})
print(f"単一崩壊 P(531)={M531}/{M531+M333}={M531/(M531+M333):.4f} (24/31=0.7742)")
W11A=runW((2,1,1,0),11); W11B=runW((-2,1,1,0),11)
print(f"s=11 W(A): { {f:round(w) for f,w in W11A.items() if w>1e-9} }")
print(f"s=11 W(B): { {f:round(w) for f,w in W11B.items() if w>1e-9} }")
seq11=seq_twin(11); bat11=batch_twin(11)
dA=derived_twin(11,W11A); dB=derived_twin(11,W11B)
print(f"Δ(11) A={TV(dA,seq11):.5f} (0.24136), B={TV(dB,seq11):.5f} (0.43971), TV(seq,bat)={TV(seq11,bat11):.5f} (0.06929)")

print()
print("==== 1: Δ(13) 親型3種×セクター + 結合多重集合衝突 ====")
fins13=all_finals(13)
joint_map=defaultdict(list)
for i,f1 in enumerate(fins13):
    for f2 in fins13[i:]:
        key=tuple(sorted(Counter(f1)+Counter(f2)))
        joint_map[key].append((f1,f2))
coll=[v for v in joint_map.values() if len(v)>1]
print(f"s=13 結合多重集合の衝突: {len(coll)} 件 {coll[:3] if coll else ''}")
seq13=seq_twin(13); bat13=batch_twin(13)
reps={'(2,1,1,1)A':(2,1,1,1),'(2,1,1,1)B':(-2,1,1,1),
      '(2,2,0,0)opp':(-2,2,0,0),'(2,2,0,0)++':(2,2,0,0),'(2,2,0,0)--':(-2,-2,0,0),
      '(3,0,0,0)A':(3,0,0,0),'(3,0,0,0)B':(-3,0,0,0)}
t0=time.time()
W13={}
for name,vp in reps.items():
    W13[name]=runW(vp,13)
    d13=derived_twin(13,W13[name])
    d13g=derived_twin(13,W13[name],graded=True)
    print(f"  {name}: Δ(13)={TV(d13,seq13):.5f}, graded版={TV(d13g,seq13):.5f}, TV(二値,graded)={TV(d13,d13g):.5f} ({time.time()-t0:.0f}s)")
print(f"  参考 TV(seq,bat) s=13 = {TV(seq13,bat13):.5f}")

print()
print("==== 2a: graded 変種 (s=9, 11) ====")
for s,Wd,name in ((9,W9,'s=9 A'),(11,W11A,'s=11 A'),(11,W11B,'s=11 B')):
    sq=seq_twin(s)
    db=derived_twin(s,Wd); dg=derived_twin(s,Wd,graded=True)
    print(f"  {name}: Δ二値={TV(db,sq):.5f}, Δgraded={TV(dg,sq):.5f}, 規約幅 TV(二値,graded)={TV(db,dg):.5f}")

print()
print("==== 2b: 完全結合振幅版 (s=9, 双子コヒーレント+セル水準排他) ====")
def eta_paths_full(vp,parent,final):
    out=[]
    for ((a,b),x,(c,d),d1) in two_step_paths(parent,final):
        for ia,va in enumerate(SH[a]):
            for ib,vb in enumerate(SH[b]):
                if a==b and not ia<ib: continue
                e1=tphase(vp,va,vb)
                if e1 is None: continue
                for vx,kept in ([(va,vb),(vb,va)] if a==b else [((va if x==a else vb),(vb if x==a else va))]):
                    for ic,vc in enumerate(SH[c]):
                        for idd,vd in enumerate(SH[d]):
                            if c==d and not ic<idd: continue
                            e2=tphase(vx,vc,vd)
                            if e2 is None: continue
                            out.append((e1*e2, frozenset([kept,vc,vd]), (kept,vc,vd)))
    return out
vp1=(2,1,0,0)
for vp2 in [(1,2,0,0),(2,-1,0,0),(0,0,2,1)]:
    P1={f:eta_paths_full(vp1,9,f) for f in [(1,3,5),(3,3,3)]}
    P2={f:eta_paths_full(vp2,9,f) for f in [(1,3,5),(3,3,3)]}
    zconf=defaultdict(complex)
    for F1 in P1:
        for F2 in P2:
            for (eta1,fs1,_) in P1[F1]:
                if len(fs1)<3: continue
                for (eta2,fs2,_) in P2[F2]:
                    if len(fs2)<3: continue
                    if fs1 & fs2: continue   # セル水準排他
                    zconf[(frozenset([fs1,fs2]) if fs1!=fs2 else fs1, tuple(sorted([F1,F2])))]+=eta1*eta2
    Pjoint=defaultdict(float)
    for (cfg,ch),z in zconf.items(): Pjoint[ch]+=abs(z)**2
    tot=sum(Pjoint.values())
    Pjoint={k:v/tot for k,v in Pjoint.items()}
    XJ=tuple(sorted([(1,3,5),(3,3,3)]))
    PX=Pjoint.get(XJ,0)
    d9=derived_twin(9,W9)
    print(f"  vp2={vp2}: 結合振幅 P(X)={PX:.5f} vs D3因子化 {float(d9[XJ]):.5f}, 全結合分布 TV={TV(Pjoint,d9):.5f}")

print()
print("==== 3: 零チャネル両側一致 (s=9,11,13) ====")
for s,Wd_list in ((9,[W9]),(11,[W11A,W11B]),(13,list(W13.values()))):
    fins=all_finals(s)
    for f in fins:
        M=M_channel(f,{})
        Wmax=max(wd.get(f,0) for wd in Wd_list)
        if (M==0)!=(Wmax<1e-6):
            print(f"  s={s} {f}: M={M}, Wmax={Wmax:.1f} ← 片側のみ零!")
    zeroM=[f for f in fins if M_channel(f,{})==0]
    zeroW=[f for f in fins if max(wd.get(f,0) for wd in Wd_list)<1e-6]
    print(f"  s={s}: M=0 チャネル {zeroM} / W=0(全代表親) {zeroW} → 台の一致 {set(zeroM)==set(zeroW)}")
