#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺65 検証: 射影と保存則の可換性
#  C1: R² 加法保存 — 全管理経路で s_parent = s(kept)+s(c)+s(d) (記録時点)、中間段は厳密に ±1 借用で相殺
#  C2: Q(ε) 乗法保存 — 振幅を運ぶ経路(e1,e2≠None)上で ε(vp)=ε(kept)·ε(vc)·ε(vd) が成立するか
#  C3: 記録からの構成読み(matched filter) — ホログラムだけから占有セル集合を復元 → s_read の加法性
import numpy as np, itertools
src=open('supplement57_convention_dependence.py',encoding='utf-8').read()
exec(src.split('[1]')[0].split('print(')[0])
def sval(v): return sum((abs(t)+0.5)**2 for t in v)
def eps(v): return (-1)**sum(abs(t) for t in v)

# ---------- C1 + C2: 経路台帳 ----------
for s in (9,11,13):
    parents = SH[s] if s!=9 else [v for v in SH[9] if sorted(map(abs,v),reverse=True)==[2,1,0,0]]
    fins = all_finals(s)
    n_paths=0; r_ok=0; borrow_ok=0; eps_ok=0; eps_bad=0
    for vp in parents[:8]:   # 代表親8(全数は対称性で同型、B4同変)
        for f in fins:
            for ((a,b),x,(c,d),d1) in two_step_paths(s,f):
                for va in SH[a]:
                    for vb in SH[b]:
                        e1=tphase(vp,va,vb)
                        if e1 is None: continue
                        vx_shell = x; kept = vb if x==a else va
                        vdecay = va if x==a else vb
                        for vc in SH[c]:
                            for vd in SH[d]:
                                e2=tphase(vdecay,vc,vd)
                                if e2 is None: continue
                                n_paths+=1
                                # R²: 記録時点の加法保存
                                if abs(sval(vp)-(sval(kept)+sval(vc)+sval(vd)))<1e-9: r_ok+=1
                                # 中間段の借用: (a+b)-s と (c+d)-x が ±1 で相殺
                                b1=(a+b)-s; b2=(c+d)-x
                                if abs(b1)==1 and b1+b2==0: borrow_ok+=1
                                # ε 乗法保存
                                if eps(vp)==eps(kept)*eps(vc)*eps(vd): eps_ok+=1
                                else: eps_bad+=1
    print(f"s={s}: 寄与経路 {n_paths}")
    print(f"  C1a R²加法(記録時点): {r_ok}/{n_paths}")
    print(f"  C1b 中間借用 ±1 相殺: {borrow_ok}/{n_paths}")
    print(f"  C2  ε乗法保存:       {eps_ok}/{n_paths} (違反 {eps_bad})")

# ---------- C3: 記録からの構成読みと s_read 加法性 ----------
# 4D 8^4 格子、辞書: k>0→√2cos(2πkx), k<0→√2sin(2π|k|x), k=0→1 (正規直交)
Ng=8; g=np.arange(Ng)/Ng
def axwave(k):
    if k==0: return np.ones(Ng)
    return np.sqrt(2)*(np.cos(2*np.pi*abs(k)*g) if k>0 else np.sin(2*np.pi*abs(k)*g))
def cellwave(v):
    w=axwave(v[0]).reshape(-1,1,1,1)*axwave(v[1]).reshape(1,-1,1,1)
    return w*axwave(v[2]).reshape(1,1,-1,1)*axwave(v[3]).reshape(1,1,1,-1)
ALL=[v for r in (1,3,5,7,9) for v in SH[r]]
import random
random.seed(11)
trials=0; recover=0; s_add=0
vp=(2,1,0,0)
fins=all_finals(9)
cases=[]
for f in fins:
    for ((a,b),x,(c,d),d1) in two_step_paths(9,f):
        for va in SH[a][:4]:
            for vb in SH[b][:4]:
                if tphase(vp,va,vb) is None: continue
                kept = vb if x==a else va; vdec = va if x==a else vb
                for vc in SH[c][:4]:
                    for vd in SH[d][:4]:
                        if tphase(vdec,vc,vd) is None: continue
                        if len({kept,vc,vd})==3: cases.append((kept,vc,vd))
random.shuffle(cases); cases=cases[:60]
for (k1,k2,k3) in cases:
    trials+=1
    I=(1.0+cellwave(k1)+cellwave(k2)+cellwave(k3))**2
    # matched filter: 全625セルとの内積(転写項 2Φ が δ を与える; 交差項が汚染候補)
    scores={v: float(np.mean(I*cellwave(v))) for v in ALL}
    top=sorted(scores,key=lambda v:-abs(scores[v]))[:3]
    if set(top)=={k1,k2,k3}:
        recover+=1
        s_read=sum(sval(v) for v in top)
        if abs(s_read-9.0)<1e-9: s_add+=1
print(f"C3: 構成読み(60無作為終状態): 占有セル完全復元 {recover}/{trials}, s_read 加法(=9) {s_add}/{trials}")
print("C3 帰結: R は記録から読んだ構成の導出量 √s_read であり、時計レート ν_t=√s_read を同時に決める(t の従属性)")
